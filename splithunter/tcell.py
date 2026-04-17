#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""
TcellExTRECT-style T-cell fraction estimation from WES/WGS BAM files.

This module now exposes two code paths:

* ``legacy-whole-window``: the original Rust median log-ratio shortcut used by
  this repo before this patch.
* ``upstream-exon-smoothed``: an exon-aware port of the original
  McGranahanLab/TcellExTRECT workflow.  It restricts the signal to capture
  exons, applies the per-exon running median filter, drops failed exons, then
  estimates the focal dip from a smoothed log-ratio track with a QC score.

The exon-aware path becomes the default whenever exon intervals are available
for the analysed locus.  This matches the original TcellExTRECT semantics for
WES/TRA runs.  Whole-window mode remains available for loci without vetted
exon sets and for explicit ``--no-default-targets`` runs.
"""

import argparse
import json
import math
import sys

import numpy as np
import pandas as pd

from . import _core
from .loci import (
    HG38_TCELL_SEGMENTS,
    LocusSegments,
    default_exons_for,
    load_exons_bed,
)
from .utils import DefaultHelpParser


def _empty_result(**extra):
    result = {
        "baseline_median": float("nan"),
        "focal_median": float("nan"),
        "log2_ratio": float("nan"),
        "tcell_fraction": float("nan"),
        "tcell_fraction_lwr": float("nan"),
        "tcell_fraction_upr": float("nan"),
        "focal_positions": 0,
        "baseline_positions": 0,
        "qc_fit": float("nan"),
        "gc_corrected": False,
    }
    result.update(extra)
    return result


def _clamp01(x):
    if x != x:
        return x
    if x < 0.0:
        return 0.0
    if x > 1.0:
        return 1.0
    return x


def _coverage_df_from_pairs(cov):
    if not cov:
        return pd.DataFrame({"pos": [], "reads": []}, dtype=float)
    pos, reads = zip(*cov)
    return pd.DataFrame({"pos": np.asarray(pos, dtype=np.int64),
                         "reads": np.asarray(reads, dtype=float)})


def _interval_mask(pos, interval):
    start, end = interval
    return (pos >= start) & (pos <= end)


def _baseline_mask(pos, segs):
    return _interval_mask(pos, segs.left_baseline) | _interval_mask(pos, segs.right_baseline)


def _running_median_constant(values, k):
    values = np.asarray(values, dtype=float)
    n = values.size
    width = 2 * int(k) + 1
    if width > n:
        if n == 0:
            width = 1
        elif n % 2 == 0:
            width = n - 1
        else:
            width = n
    if width <= 1:
        return values.copy()

    half = width // 2
    out = np.empty(n, dtype=float)
    for i in range(n):
        lo = max(0, i - half)
        hi = min(n, i + half + 1)
        window = values[lo:hi]
        need = width - window.size
        if need > 0:
            pad = np.repeat(values[0] if lo == 0 else values[-1], need)
            window = np.concatenate([pad, window]) if lo == 0 else np.concatenate([window, pad])
        out[i] = float(np.median(window))
    return out


def _select_exon_coverage(coverage_df, exons, median_k, median_thresh):
    rows = []
    removed = []
    for exon_idx, (start, end) in enumerate(exons, 1):
        tmp = coverage_df[_interval_mask(coverage_df["pos"], (start, end))].copy()
        if tmp.empty:
            removed.append(exon_idx)
            continue
        tmp["reads"] = _running_median_constant(tmp["reads"].to_numpy(), median_k)
        exon_med = float(np.nanmedian(tmp["reads"].to_numpy()))
        if (not np.isfinite(exon_med)) or exon_med < median_thresh:
            removed.append(exon_idx)
            continue
        tmp["exon"] = exon_idx
        rows.append(tmp)

    if not rows:
        empty = coverage_df.iloc[0:0].copy()
        empty["exon"] = pd.Series(dtype=int)
        return empty, removed

    exon_df = pd.concat(rows, ignore_index=True)
    exon_df = exon_df.drop_duplicates(subset=["pos"], keep="first")
    return exon_df, removed


def _adaptive_loess_with_ci(x, y, grid, frac=0.2, min_points=100):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    grid = np.asarray(grid, dtype=float)
    if x.size == 0:
        nan = np.full(grid.size, np.nan)
        return nan, nan, nan

    k = min(x.size, max(int(math.ceil(frac * x.size)), int(min_points)))
    est = np.empty(grid.size, dtype=float)
    lower = np.empty(grid.size, dtype=float)
    upper = np.empty(grid.size, dtype=float)

    for i, g in enumerate(grid):
        dist = np.abs(x - g)
        h = np.partition(dist, k - 1)[k - 1]
        if h <= 0:
            nz = dist[dist > 0]
            h = float(nz.min()) if nz.size else 1.0
        u = dist / h
        w = np.where(u < 1.0, (1.0 - u ** 3) ** 3, 0.0)
        sw = float(w.sum())
        if sw <= 0.0:
            est[i] = np.nan
            lower[i] = np.nan
            upper[i] = np.nan
            continue

        z = x - g
        s0 = sw
        s1 = float((w * z).sum())
        s2 = float((w * z * z).sum())
        t0 = float((w * y).sum())
        t1 = float((w * z * y).sum())
        det = s0 * s2 - s1 * s1

        if abs(det) < 1e-12:
            fit0 = t0 / s0
            fitted = np.full_like(y, fit0)
        else:
            fit0 = (t0 * s2 - t1 * s1) / det
            fit1 = (t1 * s0 - t0 * s1) / det
            fitted = fit0 + fit1 * z

        resid_var = float((w * np.square(y - fitted)).sum() / sw)
        eff_n = sw * sw / max(float(np.square(w).sum()), 1e-12)
        se = math.sqrt(max(resid_var, 0.0) / max(eff_n - 2.0, 1.0))

        est[i] = fit0
        lower[i] = fit0 - 1.96 * se
        upper[i] = fit0 + 1.96 * se

    return est, lower, upper


def _qc_fit_value(smoothed):
    smoothed = np.asarray(smoothed, dtype=float)
    if smoothed.size < 3 or not np.isfinite(smoothed).any():
        return float("nan")
    smoothed = smoothed[np.isfinite(smoothed)]
    if smoothed.size < 3:
        return float("nan")

    maxima = []
    minima = []
    for i in range(1, smoothed.size - 1):
        left, cur, right = smoothed[i - 1], smoothed[i], smoothed[i + 1]
        if cur >= left and cur >= right:
            maxima.append(cur)
        if cur <= left and cur <= right:
            minima.append(cur)

    extrema = np.asarray(minima + maxima, dtype=float)
    if extrema.size == 0:
        return 1.0
    denom = float(np.max(np.abs(extrema)))
    if denom <= 0.0:
        return 1.0
    return float(np.mean(np.abs(extrema / denom)))


def _gc_fraction(seq):
    seq = seq.upper()
    valid = sum(base in "ACGT" for base in seq)
    if valid == 0:
        return float("nan")
    gc = sum(base in "GC" for base in seq)
    return gc / valid


def _compute_gc_corrected_ratio(logr_df, segs, exons, reference_fasta, sliding=1000):
    if not reference_fasta:
        return logr_df, False

    try:
        import pysam
    except Exception:
        return logr_df, False

    try:
        fasta = pysam.FastaFile(reference_fasta)
    except Exception:
        return logr_df, False

    try:
        exon_gc = {}
        for exon_idx, (start, end) in enumerate(exons, 1):
            seq = fasta.fetch(segs.chrom, start, end + 1)
            exon_gc[exon_idx] = _gc_fraction(seq)

        full_seq = fasta.fetch(segs.chrom, segs.full[0], segs.full[1] + 1)
    except Exception:
        return logr_df, False
    finally:
        fasta.close()

    full_len = len(full_seq)
    starts = np.arange(0, max(full_len - sliding, 0) + 1, sliding, dtype=int)
    if starts.size == 0:
        return logr_df, False

    gc_grid = pd.DataFrame({
        "pos": segs.full[0] + starts,
        "GC": np.asarray([_gc_fraction(full_seq[s:s + sliding]) for s in starts], dtype=float),
    })
    gc_grid = gc_grid[np.isfinite(gc_grid["GC"])].copy()
    if gc_grid.empty:
        return logr_df, False

    smooth_gc, _, _ = _adaptive_loess_with_ci(
        gc_grid["pos"].to_numpy(),
        gc_grid["GC"].to_numpy(),
        logr_df["pos"].to_numpy(),
        frac=0.15,
        min_points=25,
    )

    out = logr_df.copy()
    out["exon_gc"] = out["exon"].map(exon_gc).astype(float)
    out["smooth_gc"] = smooth_gc
    design = out[["exon_gc", "smooth_gc"]].copy()
    design["exon_gc2"] = np.square(design["exon_gc"])
    design["smooth_gc2"] = np.square(design["smooth_gc"])
    design.insert(0, "intercept", 1.0)
    mask = np.isfinite(design).all(axis=1) & np.isfinite(out["Ratio"])
    if mask.sum() < 5:
        return logr_df, False

    X = design.loc[mask].to_numpy(dtype=float)
    y = out.loc[mask, "Ratio"].to_numpy(dtype=float)
    coef, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    fitted = X @ coef
    out.loc[mask, "Ratio_gc_correct"] = y - fitted
    return out, True


def _estimate_upstream_from_coverage_df(
    coverage_df,
    segs,
    exons,
    min_cov=1,
    median_k=50,
    median_thresh=15,
    reference_fasta=None,
):
    exon_df, removed = _select_exon_coverage(coverage_df, exons, median_k, median_thresh)
    if exon_df.empty or len(removed) > 30:
        return _empty_result(
            coverage_mode="upstream-exon-smoothed",
            failure_reason="no_usable_exons",
            exons_total=len(exons),
            exons_removed=len(removed),
            exons_used=max(len(exons) - len(removed), 0),
        )

    exon_df = exon_df[exon_df["reads"] >= min_cov].copy()
    if exon_df.empty:
        return _empty_result(
            coverage_mode="upstream-exon-smoothed",
            failure_reason="no_positions_after_min_cov",
            exons_total=len(exons),
            exons_removed=len(removed),
            exons_used=max(len(exons) - len(removed), 0),
        )

    baseline_rows = exon_df[_baseline_mask(exon_df["pos"], segs)].copy()
    focal_rows = exon_df[_interval_mask(exon_df["pos"], segs.focal)].copy()
    if baseline_rows.empty:
        return _empty_result(
            coverage_mode="upstream-exon-smoothed",
            failure_reason="no_baseline_positions",
            exons_total=len(exons),
            exons_removed=len(removed),
            exons_used=max(len(exons) - len(removed), 0),
            focal_positions=int(focal_rows.shape[0]),
        )

    baseline_median = float(np.median(baseline_rows["reads"]))
    focal_median = (
        float(np.median(focal_rows["reads"])) if not focal_rows.empty else float("nan")
    )

    logr_df = exon_df[["pos", "reads", "exon"]].copy()
    logr_df["Ratio"] = np.log2(logr_df["reads"] / baseline_median)
    logr_df = logr_df[np.isfinite(logr_df["Ratio"])].copy()
    if logr_df.empty:
        return _empty_result(
            coverage_mode="upstream-exon-smoothed",
            failure_reason="no_finite_logr",
            exons_total=len(exons),
            exons_removed=len(removed),
            exons_used=max(len(exons) - len(removed), 0),
            baseline_median=baseline_median,
            focal_median=focal_median,
            focal_positions=int(focal_rows.shape[0]),
            baseline_positions=int(baseline_rows.shape[0]),
        )

    logr_df, gc_corrected = _compute_gc_corrected_ratio(
        logr_df, segs, exons, reference_fasta
    )
    ratio_col = "Ratio_gc_correct" if gc_corrected else "Ratio"
    ratio_vals = logr_df[ratio_col].to_numpy(dtype=float)

    baseline_logr = logr_df.loc[_baseline_mask(logr_df["pos"], segs), ratio_col].to_numpy(dtype=float)
    baseline_adjust = float(np.mean(baseline_logr))
    if baseline_logr.size > 1:
        baseline_ci = 1.96 * float(np.std(baseline_logr, ddof=1)) / math.sqrt(baseline_logr.size)
    else:
        baseline_ci = 0.0

    grid_full = np.arange(segs.full[0], segs.full[1] + 1, 100, dtype=float)
    smooth_full, _, smooth_upper = _adaptive_loess_with_ci(
        logr_df["pos"].to_numpy(),
        ratio_vals,
        grid_full,
        frac=0.2,
        min_points=100,
    )
    qc_fit = _qc_fit_value(smooth_full)

    focal_fit = (grid_full > segs.focal[0]) & (grid_full < segs.focal[1])
    if not focal_fit.any():
        focal_fit = np.abs(grid_full - np.mean(segs.focal)) == np.min(np.abs(grid_full - np.mean(segs.focal)))
    smooth_adjust = float(np.nanmean(smooth_full[focal_fit]))
    smooth_ci = float(np.nanmean(smooth_upper[focal_fit]) - np.nanmean(smooth_full[focal_fit]))

    adjust_value = max(baseline_adjust, smooth_adjust)
    norm_ci = smooth_ci if smooth_adjust > baseline_adjust else baseline_ci

    adjusted = ratio_vals - adjust_value
    grid_focal = np.arange(segs.focal[0], segs.focal[1] + 1, 100, dtype=float)
    focal_est, focal_lower, focal_upper = _adaptive_loess_with_ci(
        logr_df["pos"].to_numpy(),
        adjusted,
        grid_focal,
        frac=0.2,
        min_points=100,
    )

    log2_ratio = float(np.nanmean(focal_est))
    log2_ratio_min = float(np.nanmean(focal_lower)) - norm_ci
    log2_ratio_max = float(np.nanmean(focal_upper)) + norm_ci

    tcell_fraction = _clamp01(1.0 - 2.0 ** log2_ratio)
    tcell_fraction_upr = _clamp01(1.0 - 2.0 ** log2_ratio_min)
    tcell_fraction_lwr = _clamp01(1.0 - 2.0 ** log2_ratio_max)

    return {
        "baseline_median": baseline_median,
        "focal_median": focal_median,
        "log2_ratio": log2_ratio,
        "tcell_fraction": tcell_fraction,
        "tcell_fraction_lwr": tcell_fraction_lwr,
        "tcell_fraction_upr": tcell_fraction_upr,
        "focal_positions": int(focal_rows.shape[0]),
        "baseline_positions": int(baseline_rows.shape[0]),
        "qc_fit": qc_fit,
        "coverage_mode": "upstream-exon-smoothed",
        "gc_corrected": bool(gc_corrected),
        "exons_total": len(exons),
        "exons_removed": len(removed),
        "exons_used": max(len(exons) - len(removed), 0),
    }


def estimate_from_coverage(
    positions,
    depths,
    segs,
    target_intervals=None,
    min_cov=1,
    median_k=50,
    median_thresh=15,
    reference_fasta=None,
):
    coverage_df = pd.DataFrame({
        "pos": np.asarray(positions, dtype=np.int64),
        "reads": np.asarray(depths, dtype=float),
    })
    if target_intervals:
        return _estimate_upstream_from_coverage_df(
            coverage_df,
            segs,
            list(target_intervals),
            min_cov=min_cov,
            median_k=median_k,
            median_thresh=median_thresh,
            reference_fasta=reference_fasta,
        )

    result = dict(_core.fraction_from_coverage(
        list(map(int, positions)),
        list(map(int, depths)),
        segs.focal[0], segs.focal[1],
        segs.left_baseline[0], segs.left_baseline[1],
        segs.right_baseline[0], segs.right_baseline[1],
        min_cov,
    ))
    result["coverage_mode"] = "legacy-whole-window"
    result["gc_corrected"] = False
    result["qc_fit"] = float("nan")
    return result


def estimate(
    bam_path: str,
    segs: LocusSegments,
    min_mapq: int = 1,
    min_cov: int = 1,
    target_mask=None,
    reference_fasta=None,
    median_k: int = 50,
    median_thresh: int = 15,
):
    """
    Return the TcellExTRECT-style summary dict for one locus.

    When ``target_mask`` is provided it is interpreted as the exon/capture
    intervals that define the upstream TcellExTRECT signal substrate.  In that
    mode the estimator ignores intronic coverage, smooths depth within each
    exon, removes low-depth exons, and estimates the focal dip from the exon
    log-ratio track.  If ``target_mask`` is absent the function falls back to
    the repo's original whole-window Rust shortcut.
    """
    mask_list = list(target_mask) if target_mask else None
    if mask_list:
        cov = _core.region_coverage(
            bam_path,
            segs.chrom,
            segs.full[0],
            segs.full[1],
            min_mapq,
        )
        coverage_df = _coverage_df_from_pairs(cov)
        return _estimate_upstream_from_coverage_df(
            coverage_df,
            segs,
            mask_list,
            min_cov=min_cov,
            median_k=median_k,
            median_thresh=median_thresh,
            reference_fasta=reference_fasta,
        )

    result = _core.tcell_fraction(
        bam_path,
        segs.chrom,
        segs.focal[0],
        segs.focal[1],
        segs.left_baseline[0],
        segs.left_baseline[1],
        segs.right_baseline[0],
        segs.right_baseline[1],
        min_mapq,
        min_cov,
        [],
    )
    out = dict(result)
    out["coverage_mode"] = "legacy-whole-window"
    out["gc_corrected"] = False
    out["qc_fit"] = float("nan")
    return out


def _resolve_contig(bam_path, contig):
    import pysam

    with pysam.AlignmentFile(bam_path) as bf:
        refs = set(bf.references)
    if contig in refs:
        return contig
    alt = contig[3:] if contig.startswith("chr") else f"chr{contig}"
    if alt in refs:
        return alt
    raise KeyError(f"Contig {contig!r} (or {alt!r}) not present in BAM header")


def main(args=None):
    p = DefaultHelpParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("bam", help="Input BAM (indexed)")
    p.add_argument("--locus", default="TRA",
                   choices=sorted(HG38_TCELL_SEGMENTS.keys()),
                   help="V(D)J locus to analyse")
    p.add_argument("--min-mapq", type=int, default=1)
    p.add_argument("--min-cov", type=int, default=1)
    p.add_argument("--targets", default=None,
                   help="Capture-kit exon BED.  These intervals define the "
                        "upstream TcellExTRECT exon signal track.  Recommended "
                        "for WES BAMs.  TRA uses a packaged default exon set "
                        "if --targets is omitted.")
    p.add_argument("--no-default-targets", action="store_true",
                   help="Disable the packaged default exon set and use the "
                        "legacy whole-window estimator (appropriate for WGS "
                        "or for non-exome coverage tracks).")
    p.add_argument("--reference", default=None,
                   help="Indexed FASTA reference.  If provided together with "
                        "exon targets, enable the GC-corrected upstream path.")
    p.add_argument("--json", action="store_true",
                   help="Emit a JSON object instead of a human summary")
    parsed = p.parse_args(args)

    segs = HG38_TCELL_SEGMENTS[parsed.locus]
    # Accept BAMs indexed with or without the 'chr' prefix.
    try:
        contig = _resolve_contig(parsed.bam, segs.chrom)
    except Exception:
        contig = segs.chrom
    segs = segs._replace(chrom=contig)

    if parsed.targets:
        mask = load_exons_bed(parsed.targets, contig)
    elif parsed.no_default_targets:
        mask = None
    else:
        mask = default_exons_for(parsed.locus, contig)

    result = estimate(
        parsed.bam,
        segs,
        parsed.min_mapq,
        parsed.min_cov,
        mask,
        reference_fasta=parsed.reference,
    )
    result["locus"] = parsed.locus
    result["bam"] = parsed.bam
    result["mask_intervals"] = 0 if not mask else len(mask)

    if parsed.json:
        json.dump(result, sys.stdout, indent=2, sort_keys=True, default=float)
        sys.stdout.write("\n")
        return

    print(f"Locus          : {parsed.locus} ({segs.chrom})")
    print(f"Focal window   : {segs.focal[0]:,}-{segs.focal[1]:,}")
    print(f"Baseline L / R : {segs.left_baseline[0]:,}-{segs.left_baseline[1]:,} "
          f"/ {segs.right_baseline[0]:,}-{segs.right_baseline[1]:,}")
    if mask:
        print(f"Target exons   : {len(mask)} intervals "
              f"({result.get('coverage_mode', 'upstream-exon-smoothed')})")
    else:
        print(f"Target exons   : none (legacy whole-window coverage)")
    print(f"Depth (focal)  : median={result['focal_median']:.2f} "
          f"(n={result['focal_positions']})")
    print(f"Depth (base)   : median={result['baseline_median']:.2f} "
          f"(n={result['baseline_positions']})")
    print(f"log2 ratio     : {result['log2_ratio']:.4f}")
    print(f"T-cell fraction: {result['tcell_fraction']:.4f} "
          f"[{result['tcell_fraction_lwr']:.4f}, {result['tcell_fraction_upr']:.4f}]")
    if result.get("qc_fit") == result.get("qc_fit"):
        print(f"QC fit         : {result['qc_fit']:.4f}")


if __name__ == "__main__":
    main(sys.argv[1:])
