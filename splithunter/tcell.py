#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""
TcellExTRECT-style T-cell fraction estimation from WES/WGS BAM files.

Given a BAM and a TCRA (or any V(D)J) locus definition, compute the log2
coverage ratio of the focal V-J region against flanking baseline windows and
derive the T-cell (or B-cell) fraction as ``1 - 2^logR``.

The coverage and ratio computations run in Rust via ``splithunter._core``.
"""

import argparse
import json
import sys

from . import _core
from .loci import HG38_TCELL_SEGMENTS, HG38_LOCI_BY_NAME, LocusSegments
from .utils import DefaultHelpParser


def estimate(bam_path: str, segs: LocusSegments, min_mapq: int = 1, min_cov: int = 1):
    """
    Return the TcellExTRECT-style summary dict for one locus.
    """
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
    )
    return dict(result)


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

    result = estimate(parsed.bam, segs, parsed.min_mapq, parsed.min_cov)
    result["locus"] = parsed.locus
    result["bam"] = parsed.bam

    if parsed.json:
        json.dump(result, sys.stdout, indent=2, sort_keys=True, default=float)
        sys.stdout.write("\n")
        return

    print(f"Locus          : {parsed.locus} ({segs.chrom})")
    print(f"Focal window   : {segs.focal[0]:,}-{segs.focal[1]:,}")
    print(f"Baseline L / R : {segs.left_baseline[0]:,}-{segs.left_baseline[1]:,} "
          f"/ {segs.right_baseline[0]:,}-{segs.right_baseline[1]:,}")
    print(f"Depth (focal)  : median={result['focal_median']:.2f} "
          f"(n={result['focal_positions']})")
    print(f"Depth (base)   : median={result['baseline_median']:.2f} "
          f"(n={result['baseline_positions']})")
    print(f"log2 ratio     : {result['log2_ratio']:.4f}")
    print(f"T-cell fraction: {result['tcell_fraction']:.4f} "
          f"[{result['tcell_fraction_lwr']:.4f}, {result['tcell_fraction_upr']:.4f}]")


if __name__ == "__main__":
    main(sys.argv[1:])
