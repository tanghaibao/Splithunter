#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""
Compile splithunter calls into a tsv file from a set of JSON files.
"""

import argparse
import datetime as _dt
import html as _html
import json
import math
import os.path as op
import sys

from multiprocessing import Pool, cpu_count

import pandas as pd

from .loci import HG38_LOCI, HG38_TCELL_SEGMENTS
from .utils import DefaultHelpParser


def df_to_tsv(df, tsvfile, index=False):
    dd = ["SampleKey"]
    columns = dd + sorted([x for x in df.columns if x not in dd])

    df = df.reindex(columns=columns)
    df = df.sort_values("SampleKey")
    df.to_csv(tsvfile, sep='\t', index=index)
    print("TSV output written to `{}` (# samples={})".format(tsvfile, df.shape[0]),
          file=sys.stderr)


def json_to_df_worker(jsonfile):
    with open(jsonfile) as fp:
        return json.load(fp)


def json_to_df(jsonfiles, cpus):
    """
    Compile a number of json files into a DataFrame.
    """
    with Pool(processes=cpus) as pool:
        records = pool.map(json_to_df_worker, jsonfiles)
    return pd.DataFrame.from_records(records)


def get_midpoint(s):
    """
    >>> get_midpoint("1:58,460-58,557(+)")
    (58508, '+')
    """
    _, startendstrand = s.split(':')
    startend, strand = startendstrand.split('(')
    start, end = startend.replace(',', '').split('-')
    start, end = int(start), int(end)
    strand = strand.split(')')[0]
    return (start + end) // 2, strand


def parse_alignments(s):
    """Convert alignment strings to tuples for plotting.

    >>> parse_alignments("1:58,460-58,557(+)|1:135,876-135,929(-);1:58,460-58,557(+)|1:135,876-135,929(-);")
    [(58508, 135902, '+-'), (58508, 135902, '+-')]
    """
    res = []
    for event in s.split(';'):
        if event.strip() == "":
            continue
        left, right = event.split("|")
        leftmid, leftstrand = get_midpoint(left)
        rightmid, rightstrand = get_midpoint(right)
        res.append((leftmid, rightmid, leftstrand + rightstrand))
    return res


def slicing_filter(details, orientation, offset=0, boundary=800000):
    sliced = [x for x in details if ((x[0] - offset) > boundary)
                                 or ((x[1] - offset) > boundary)]
    sliced = [x for x in sliced if x[-1] in orientation]
    if offset == 0:  # SR
        sliced = [x for x in sliced if x[0] > x[1]]
    return len(sliced)


def filter_TRA(df, TRA_start=21621904):
    data = {}
    for _, row in df.iterrows():
        samplekey = row["SampleKey"]
        SP, SR = row["TRA.SP-DETAILS"], row["TRA.SR-DETAILS"]
        SP_total, SR_total = row["TRA.SP-TOTAL"], row["TRA.SR-TOTAL"]
        SP_new = SR_new = 0
        if isinstance(SP, str):
            SP_new = slicing_filter(parse_alignments(SP), ['-+'], offset=TRA_start)
        if isinstance(SR, str):
            SR_new = slicing_filter(parse_alignments(SR), ['++'])

        SR_PPM = SR_new * 1e6 / SR_total if SR_total else 0
        SP_PPM = SP_new * 1e6 / SP_total if SP_total else 0
        total = SR_total + SP_total * 2
        PPM = (SR_new + SP_new * 2) * 1e6 / total if total else 0
        data[samplekey] = [SR_PPM, SP_PPM, PPM]

    xf = pd.DataFrame.from_dict(data, orient="index")
    xf.columns = ["TRA.SR-PPM", "TRA.SP-PPM", "TRA.PPM"]
    return xf


# ---------- HTML report helpers ------------------------------------------

_HTML_CSS = """
body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
       margin: 2rem auto; max-width: 1080px; color: #222; line-height: 1.5;
       padding: 0 1rem; }
h1 { border-bottom: 2px solid #2c3e50; padding-bottom: .4rem; }
h2 { color: #2c3e50; margin-top: 2rem; border-bottom: 1px solid #eee;
     padding-bottom: .2rem; }
h3 { color: #34495e; margin-top: 1.4rem; }
table { border-collapse: collapse; width: 100%; margin-top: 1rem;
        font-size: 0.92rem; }
th, td { border: 1px solid #ddd; padding: .45rem .65rem; text-align: right; }
th { background: #f4f6f8; text-align: center; }
td.locus, td.label { text-align: left; font-weight: 600; }
td.mono { font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
          font-size: 0.82rem; text-align: left; }
tr.na td { color: #aaa; background: #fafafa; }
.meta { background: #f8f9fa; padding: 1rem 1.2rem;
        border-left: 4px solid #3498db; border-radius: 0 4px 4px 0; }
.meta dt { font-weight: 600; display: inline-block; min-width: 140px; }
.bar { background: #ecf0f1; height: 14px; border-radius: 7px;
       overflow: hidden; min-width: 120px; }
.bar > span { display: block; height: 100%;
              background: linear-gradient(90deg,#27ae60,#2980b9); }
.grid { display: grid;
        grid-template-columns: repeat(auto-fit, minmax(180px,1fr));
        gap: 0.9rem; margin-top: 1rem; }
.card { border: 1px solid #e1e4e8; border-radius: 6px; padding: .9rem 1rem;
        background: #fff; box-shadow: 0 1px 2px rgba(0,0,0,.04); }
.card .k { color: #7f8c8d; font-size: 0.78rem; text-transform: uppercase;
           letter-spacing: 0.04em; }
.card .v { font-size: 1.55rem; font-weight: 700; color: #2c3e50;
           margin-top: .25rem; }
.card .sub { color: #95a5a6; font-size: 0.82rem; }
.legend { font-size: 0.82rem; color: #555; margin-top: .3rem; }
.legend span { display: inline-block; width: 12px; height: 12px;
               border-radius: 2px; vertical-align: middle; margin-right: 4px; }
code { background: #f4f6f8; padding: .1rem .35rem; border-radius: 3px;
       font-size: 0.88em; }
"""


def _isnan(x):
    try:
        return math.isnan(float(x))
    except (TypeError, ValueError):
        return x is None


def _fmt(x, spec=".4f", na="NA"):
    if x is None or _isnan(x):
        return na
    try:
        return format(float(x), spec)
    except (TypeError, ValueError):
        return str(x)


def _bin_coverage(cov, start, end, nbins=600):
    binw = (end - start) / nbins
    sums = [0.0] * nbins
    counts = [0] * nbins
    for p, d in cov:
        i = int((p - start) / binw)
        if 0 <= i < nbins:
            sums[i] += d
            counts[i] += 1
    xs, ys = [], []
    for i in range(nbins):
        if counts[i]:
            xs.append(int(start + (i + 0.5) * binw))
            ys.append(sums[i] / counts[i])
    return xs, ys


def _coverage_svg(xs, ys, full, focal, left, right, focal_med, base_med):
    if not xs:
        return ""
    W, H, PL, PR, PT, PB = 1000, 280, 70, 30, 30, 50
    PW, PH = W - PL - PR, H - PT - PB
    xmin, xmax = full
    ymax = max(ys) * 1.05 or 1.0

    def sx(x): return PL + (x - xmin) / (xmax - xmin) * PW
    def sy(y): return PT + PH - (y / ymax) * PH

    pts = " ".join(f"{sx(x):.1f},{sy(y):.1f}" for x, y in zip(xs, ys))
    poly = (f"{sx(xs[0]):.1f},{sy(0):.1f} " + pts +
            f" {sx(xs[-1]):.1f},{sy(0):.1f}")

    def rect(x0, x1, color, opa, label):
        return (f'<rect x="{sx(x0):.1f}" y="{PT}" '
                f'width="{sx(x1)-sx(x0):.1f}" height="{PH}" '
                f'fill="{color}" fill-opacity="{opa}"/>'
                f'<text x="{(sx(x0)+sx(x1))/2:.1f}" y="{PT-8}" '
                f'text-anchor="middle" font-size="11" fill="{color}" '
                f'font-weight="600">{label}</text>')

    yticks = []
    step = max(1, int(round(ymax / 5 / 5)) * 5)
    v = 0
    while v <= ymax:
        y = sy(v)
        yticks.append(
            f'<line x1="{PL-5}" x2="{PL}" y1="{y:.1f}" y2="{y:.1f}" stroke="#666"/>'
            f'<text x="{PL-8}" y="{y+3:.1f}" text-anchor="end" font-size="10" '
            f'fill="#555">{v}</text>'
            f'<line x1="{PL}" x2="{W-PR}" y1="{y:.1f}" y2="{y:.1f}" stroke="#eee"/>')
        v += step

    xticks = []
    span = xmax - xmin
    tick = 1e5 if span < 2e6 else 5e5
    v = math.ceil(xmin / tick) * tick
    while v <= xmax:
        x = sx(v)
        xticks.append(
            f'<line x1="{x:.1f}" x2="{x:.1f}" y1="{PT+PH}" y2="{PT+PH+5}" stroke="#666"/>'
            f'<text x="{x:.1f}" y="{PT+PH+18}" text-anchor="middle" '
            f'font-size="9" fill="#555">{v/1e6:.2f}</text>')
        v += tick * 2

    med_lines = ""
    if focal_med is not None and not _isnan(focal_med):
        med_lines += (f'<line x1="{sx(focal[0]):.1f}" x2="{sx(focal[1]):.1f}" '
                      f'y1="{sy(focal_med):.1f}" y2="{sy(focal_med):.1f}" '
                      f'stroke="#c0392b" stroke-width="2.5" stroke-dasharray="6 3"/>')
    if base_med is not None and not _isnan(base_med):
        for w in (left, right):
            med_lines += (f'<line x1="{sx(w[0]):.1f}" x2="{sx(w[1]):.1f}" '
                          f'y1="{sy(base_med):.1f}" y2="{sy(base_med):.1f}" '
                          f'stroke="#2c3e50" stroke-width="2.5" stroke-dasharray="6 3"/>')

    return (f'<svg viewBox="0 0 {W} {H}" xmlns="http://www.w3.org/2000/svg" '
            f'style="width:100%;height:auto;">'
            f'<rect x="{PL}" y="{PT}" width="{PW}" height="{PH}" fill="#fafbfc" stroke="#ddd"/>'
            f'{rect(left[0], left[1], "#2c3e50", 0.10, "left baseline")}'
            f'{rect(focal[0], focal[1], "#c0392b", 0.13, "focal V-J")}'
            f'{rect(right[0], right[1], "#2c3e50", 0.10, "right baseline")}'
            f'{"".join(yticks)}{"".join(xticks)}'
            f'<polygon points="{poly}" fill="#3498db" fill-opacity="0.35"/>'
            f'<polyline points="{pts}" fill="none" stroke="#2980b9" stroke-width="1.2"/>'
            f'{med_lines}'
            f'<text x="{PL-50}" y="{PT+PH/2}" font-size="11" fill="#555" '
            f'transform="rotate(-90 {PL-50} {PT+PH/2})" text-anchor="middle">'
            f'depth (reads)</text>'
            f'<text x="{W/2}" y="{H-10}" font-size="11" fill="#555" '
            f'text-anchor="middle">position (Mb)</text>'
            f'</svg>')


def _arc_svg(events_sr, events_sp, full, focal, left, right):
    if not events_sr and not events_sp:
        return ""
    W, H, PL, PR, PT, PB = 1000, 220, 70, 30, 20, 50
    PW = W - PL - PR
    xmin, xmax = full
    base_y = H - PB

    def sx(x): return PL + (x - xmin) / (xmax - xmin) * PW

    def arc(x1, x2, color, sw, opa):
        a, b = sorted((sx(x1), sx(x2)))
        mid = (a + b) / 2
        h = min(160, (b - a) * 0.35 + 30)
        return (f'<path d="M {a:.1f} {base_y} Q {mid:.1f} {base_y-h:.1f} '
                f'{b:.1f} {base_y}" fill="none" stroke="{color}" '
                f'stroke-width="{sw}" stroke-opacity="{opa}"/>'
                f'<circle cx="{a:.1f}" cy="{base_y}" r="3" fill="{color}"/>'
                f'<circle cx="{b:.1f}" cy="{base_y}" r="3" fill="{color}"/>')

    arcs = []
    for ev in events_sp:
        arcs.append(arc(ev[0], ev[1], "#e67e22", 2.2, 0.9))
    for ev in events_sr:
        arcs.append(arc(ev[0], ev[1], "#9b59b6", 1.6, 0.75))

    xticks = []
    span = xmax - xmin
    tick = 1e5 if span < 2e6 else 5e5
    v = math.ceil(xmin / tick) * tick
    while v <= xmax:
        x = sx(v)
        xticks.append(
            f'<line x1="{x:.1f}" x2="{x:.1f}" y1="{base_y}" y2="{base_y+5}" stroke="#666"/>'
            f'<text x="{x:.1f}" y="{base_y+18}" text-anchor="middle" '
            f'font-size="9" fill="#555">{v/1e6:.2f}</text>')
        v += tick * 2

    def ar(x0, x1, color, opa, label):
        return (f'<rect x="{sx(x0):.1f}" y="{PT}" '
                f'width="{sx(x1)-sx(x0):.1f}" height="{base_y-PT}" '
                f'fill="{color}" fill-opacity="{opa}"/>'
                f'<text x="{(sx(x0)+sx(x1))/2:.1f}" y="{PT+12}" '
                f'text-anchor="middle" font-size="10" fill="{color}" '
                f'font-weight="600">{label}</text>')

    return (f'<svg viewBox="0 0 {W} {H}" xmlns="http://www.w3.org/2000/svg" '
            f'style="width:100%;height:auto;">'
            f'{ar(left[0], left[1], "#2c3e50", 0.08, "left")}'
            f'{ar(focal[0], focal[1], "#c0392b", 0.10, "V-J")}'
            f'{ar(right[0], right[1], "#2c3e50", 0.08, "right")}'
            f'<line x1="{PL}" x2="{W-PR}" y1="{base_y}" y2="{base_y}" '
            f'stroke="#333" stroke-width="1"/>'
            f'{"".join(arcs)}{"".join(xticks)}'
            f'<text x="{W/2}" y="{H-10}" font-size="11" fill="#555" '
            f'text-anchor="middle">arcs link breakpoint pairs '
            f'(orange = SP, purple = SR)</text></svg>')


def _try_coverage(bam, chrom, start, end):
    try:
        from . import _core
        return _core.region_coverage(bam, chrom, start, end, 1)
    except Exception:
        return None


def _resolve_contig(bam, chrom):
    try:
        import pysam
        with pysam.AlignmentFile(bam) as bf:
            refs = set(bf.references)
        if chrom in refs:
            return chrom
        alt = chrom[3:] if chrom.startswith("chr") else f"chr{chrom}"
        if alt in refs:
            return alt
    except Exception:
        pass
    return chrom


def _events_with_midpoints(detail_str):
    try:
        evs = parse_alignments(detail_str)
    except Exception:
        return []
    return [(a, b, s) for (a, b, s) in evs]


def _render_locus_block(record, locus_name, locus_obj, segs):
    bam = record.get("bam", "")
    pre = f"{locus_name}."
    sr_total = record.get(pre + "SR-TOTAL")
    sr_signal = record.get(pre + "SR-SIGNAL")
    sr_ppm = record.get(pre + "SR-PPM")
    sp_total = record.get(pre + "SP-TOTAL")
    sp_signal = record.get(pre + "SP-SIGNAL")
    sp_ppm = record.get(pre + "SP-PPM")
    sr_details = record.get(pre + "SR-DETAILS", "") or ""
    sp_details = record.get(pre + "SP-DETAILS", "") or ""

    has_signal = any(v is not None for v in (sr_total, sp_total))
    if not has_signal:
        return ""

    sr_events = _events_with_midpoints(sr_details)
    sp_events = _events_with_midpoints(sp_details)

    focal_med = record.get(pre + "TCELL-focal_median")
    base_med = record.get(pre + "TCELL-baseline_median")
    log2r = record.get(pre + "TCELL-log2_ratio")
    tfrac = record.get(pre + "TCELL-tcell_fraction")
    tlwr = record.get(pre + "TCELL-tcell_fraction_lwr")
    tupr = record.get(pre + "TCELL-tcell_fraction_upr")
    fpos = record.get(pre + "TCELL-focal_positions")
    bpos = record.get(pre + "TCELL-baseline_positions")

    parts = [f'<h2>Locus {locus_name} — {locus_obj.chrom}:'
             f'{locus_obj.start:,}-{locus_obj.end:,}</h2>']

    cards = [
        ("SR signal", _fmt(sr_signal, "g"),
         f"of {_fmt(sr_total, ',.0f')} reads ({_fmt(sr_ppm, '.2f')} PPM)"),
        ("SP signal", _fmt(sp_signal, "g"),
         f"of {_fmt(sp_total, ',.0f')} pairs ({_fmt(sp_ppm, '.2f')} PPM)"),
    ]
    if tfrac is not None:
        pct = "NA" if _isnan(tfrac) else f"{float(tfrac)*100:.2f}%"
        ci = ("NA" if _isnan(tlwr) or _isnan(tupr)
              else f"95% CI {float(tlwr)*100:.2f} – {float(tupr)*100:.2f}%")
        cards += [("T-cell fraction", pct, ci),
                  ("log₂ coverage ratio", _fmt(log2r), "focal vs baseline")]
    parts.append('<div class="grid">' + "".join(
        f'<div class="card"><div class="k">{k}</div>'
        f'<div class="v">{v}</div><div class="sub">{s}</div></div>'
        for k, v, s in cards) + '</div>')

    parts.append('<h3>SR / SP summary</h3>'
                 '<table><thead><tr><th>Class</th>'
                 '<th>Total reads / pairs</th><th>Signal events</th>'
                 '<th>PPM</th></tr></thead><tbody>'
                 f'<tr><td class="label">SR (split reads)</td>'
                 f'<td>{_fmt(sr_total, ",.0f")}</td>'
                 f'<td>{_fmt(sr_signal, "g")}</td>'
                 f'<td>{_fmt(sr_ppm, ".2f")}</td></tr>'
                 f'<tr><td class="label">SP (split pairs)</td>'
                 f'<td>{_fmt(sp_total, ",.0f")}</td>'
                 f'<td>{_fmt(sp_signal, "g")}</td>'
                 f'<td>{_fmt(sp_ppm, ".2f")}</td></tr>'
                 '</tbody></table>')

    if sr_events or sp_events:
        rows = []
        for i, (a, b, st) in enumerate(sp_events, 1):
            rows.append(f'<tr><td>{i}</td><td class="mono">SP</td>'
                        f'<td class="mono">{a:,}</td>'
                        f'<td class="mono">{b:,}</td><td>{_html.escape(st)}</td></tr>')
        for i, (a, b, st) in enumerate(sr_events, 1):
            rows.append(f'<tr><td>{i}</td><td class="mono">SR</td>'
                        f'<td class="mono">{a:,}</td>'
                        f'<td class="mono">{b:,}</td><td>{_html.escape(st)}</td></tr>')
        parts.append('<h3>Breakpoint detail</h3>'
                     '<table><thead><tr><th>#</th><th>Class</th>'
                     '<th>Left midpoint</th><th>Right midpoint</th>'
                     '<th>Strand</th></tr></thead><tbody>'
                     + "".join(rows) + '</tbody></table>')

    if segs is not None:
        parts.append('<h3>TcellExTRECT-style coverage</h3>')
        if tfrac is not None and not _isnan(tfrac):
            pct = f"{float(tfrac)*100:.2f}%"
            bar = (f'<div class="bar"><span style="width:'
                   f'{max(0, min(100, float(tfrac)*100)):.2f}%"></span></div>')
        else:
            pct, bar = "NA", "—"
        parts.append(
            '<table><thead><tr>'
            '<th>Focal med depth</th><th>Focal positions</th>'
            '<th>Baseline med depth</th><th>Baseline positions</th>'
            '<th>log₂ ratio</th><th>T-cell fraction</th>'
            '<th>95% CI</th><th>Visual</th></tr></thead><tbody><tr>'
            f'<td>{_fmt(focal_med, ".2f")}</td>'
            f'<td>{_fmt(fpos, ",.0f")}</td>'
            f'<td>{_fmt(base_med, ".2f")}</td>'
            f'<td>{_fmt(bpos, ",.0f")}</td>'
            f'<td>{_fmt(log2r)}</td><td>{pct}</td>'
            f'<td>{_fmt(tlwr, ".4f")} – {_fmt(tupr, ".4f")}</td>'
            f'<td>{bar}</td></tr></tbody></table>')

        if bam and op.exists(bam):
            chrom = _resolve_contig(bam, segs.chrom)
            cov = _try_coverage(bam, chrom, segs.full[0], segs.full[1])
            if cov:
                xs, ys = _bin_coverage(cov, segs.full[0], segs.full[1])
                parts.append('<h3>Coverage track</h3>')
                parts.append(_coverage_svg(
                    xs, ys, segs.full, segs.focal,
                    segs.left_baseline, segs.right_baseline,
                    focal_med, base_med))
                parts.append('<div class="legend">'
                             '<span style="background:#3498db;opacity:.5"></span>read depth · '
                             '<span style="background:#c0392b"></span>focal V-J · '
                             '<span style="background:#2c3e50"></span>baseline · '
                             'dashed = median</div>')

    if sr_events or sp_events:
        parts.append('<h3>Breakpoint arc plot</h3>')
        plot_full = segs.full if segs else (locus_obj.start, locus_obj.end)
        plot_focal = segs.focal if segs else plot_full
        plot_left = segs.left_baseline if segs else (plot_full[0], plot_full[0])
        plot_right = segs.right_baseline if segs else (plot_full[1], plot_full[1])
        parts.append(_arc_svg(sr_events, sp_events, plot_full,
                              plot_focal, plot_left, plot_right))
        parts.append('<div class="legend">'
                     '<span style="background:#e67e22"></span>SP (split pair) · '
                     '<span style="background:#9b59b6"></span>SR (split read)</div>')

    return "\n".join(parts)


def _render_html_report(records, html_path):
    now = _dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    body = ['<!DOCTYPE html><html lang="en"><head><meta charset="UTF-8">'
            '<title>Splithunter report</title>'
            f'<style>{_HTML_CSS}</style></head><body>'
            '<h1>Splithunter report</h1>'
            '<p>V(D)J split-read / split-pair detection plus '
            'TcellExTRECT-style coverage analysis. Engine: '
            '<code>splithunter._core</code> (Rust + BWA-MEM FFI).</p>'
            f'<div class="meta"><dl><dt>Generated</dt><dd>{now}</dd><br>'
            f'<dt>Samples</dt><dd>{len(records)}</dd></dl></div>']

    for rec in records:
        sk = rec.get("SampleKey", "?")
        bam = rec.get("bam", "")
        body.append(f'<h2 style="border-bottom:2px solid #2c3e50">'
                    f'Sample: {_html.escape(str(sk))}</h2>')
        body.append(f'<p>BAM: <code>{_html.escape(str(bam))}</code></p>')
        for locus in HG38_LOCI:
            block = _render_locus_block(
                rec, locus.name, locus, HG38_TCELL_SEGMENTS.get(locus.name))
            if block:
                body.append(block)

    body.append('</body></html>')
    with open(html_path, "w") as fw:
        fw.write("\n".join(body))
    print(f"HTML report written to `{html_path}`", file=sys.stderr)


def main(args=None):
    p = DefaultHelpParser(
        description=__doc__,
        prog=__file__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("files", nargs="*")
    p.add_argument('--tsv', default="out.tsv", help="Path to the tsv file")
    p.add_argument('--html', default=None,
                   help="If set, also emit an HTML report with SR/SP "
                        "breakpoint arcs and TcellExTRECT coverage tracks")
    p.add_argument('--cpus', default=cpu_count(), type=int,
                   help='Number of threads')
    args = p.parse_args(args)

    files = args.files
    tsvfile = args.tsv

    if files:
        nfiles = len(files)
        cpus = min(nfiles, args.cpus)
        print("Using {} cpus to parse {} JSON files".format(cpus, nfiles),
              file=sys.stderr)
        df = json_to_df(files, cpus)
        df_to_tsv(df, tsvfile)
    else:
        if op.exists(tsvfile):
            df = pd.read_csv(tsvfile, sep="\t")
        else:
            sys.exit(not p.print_help())

    if df.empty:
        print("Dataframe empty - check input files", file=sys.stderr)
        sys.exit(1)

    print("Filtering TRA locus for age prediction", file=sys.stderr)
    xf = filter_TRA(df)
    tra_tsv = tsvfile.rsplit(".", 1)[0] + ".TRA.tsv"
    xf = xf.sort_index()
    xf.to_csv(tra_tsv, sep='\t', index_label="SampleKey")
    print("TSV output written to `{}` (# samples={})".format(tra_tsv, xf.shape[0]),
          file=sys.stderr)

    if args.html:
        records = df.to_dict(orient="records")
        _render_html_report(records, args.html)


if __name__ == '__main__':
    main(sys.argv[1:])
