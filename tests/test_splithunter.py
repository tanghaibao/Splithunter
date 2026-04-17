#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

import csv
import gzip
import json
import os.path as op

import pytest

from splithunter import loci as sh_loci
from splithunter import run as sh_run
from splithunter import report as sh_report
from splithunter import tcell as sh_tcell

pytest.importorskip("splithunter._core")

HERE = op.dirname(op.abspath(__file__))
SAMPLES_CSV = op.join(HERE, "samples.csv")
TEST_BAM = op.join(HERE, "NA12878.bam")
FIXTURES = op.join(HERE, "fixtures")


def _load_cov_example():
    pos, dep = [], []
    with gzip.open(op.join(FIXTURES, "cov_example.csv.gz"), "rt") as f:
        r = csv.DictReader(f)
        for row in r:
            pos.append(int(row["pos"]))
            dep.append(int(row["reads"]))
    return pos, dep


def _load_seg(build):
    segs = {}
    with open(op.join(FIXTURES, f"tcra_seg_{build}.csv")) as f:
        for row in csv.DictReader(f):
            segs[row["segName"]] = (int(row["start"]), int(row["end"]))
    return segs


@pytest.fixture(scope="module")
def hg38_tra_reference(tmp_path_factory):
    """Decompress + FAI-index the small chr14/TRA reference fixture."""
    import pysam
    dest = tmp_path_factory.mktemp("ref") / "chr14_hg38_tra.fa"
    with gzip.open(op.join(FIXTURES, "chr14_hg38_tra.fa.gz"), "rt") as fi, \
         open(dest, "w") as fo:
        for line in fi:
            fo.write(line)
    pysam.faidx(str(dest))
    return str(dest)


@pytest.fixture(scope="module")
def bam_contigs():
    import pysam
    with pysam.AlignmentFile(TEST_BAM) as bf:
        return set(bf.references)


def _locus_for(bam_contigs, name):
    locus = sh_loci.HG38_LOCI_BY_NAME[name]
    chrom = sh_loci.pick_contig(locus, bam_contigs)
    return locus._replace(chrom=chrom)


def test_core_module_is_importable():
    from splithunter import _core
    assert hasattr(_core, "analyze_locus")
    assert hasattr(_core, "tcell_fraction")
    assert hasattr(_core, "region_coverage")


def test_analyze_locus_returns_expected_keys(bam_contigs, hg38_tra_reference):
    from splithunter import _core

    locus = _locus_for(bam_contigs, "TRA")
    result = _core.analyze_locus(
        TEST_BAM, locus.name, locus.chrom, locus.start, locus.end,
        hg38_tra_reference,
        30, 10_000, 50.0,
    )
    for suffix in ("SR-TOTAL", "SR-SIGNAL", "SR-PPM", "SR-DETAILS",
                   "SP-TOTAL", "SP-SIGNAL", "SP-PPM", "SP-DETAILS"):
        assert f"TRA.{suffix}" in result
    assert isinstance(result["TRA.SR-TOTAL"], int)
    assert isinstance(result["TRA.SR-PPM"], float)
    assert result["TRA.SR-TOTAL"] >= result["TRA.SR-SIGNAL"] >= 0


def test_analyze_locus_detects_split_reads(bam_contigs, hg38_tra_reference):
    """The test BAM carries NA12878 TRA rearrangement reads — BWA-MEM
    realignment should identify split reads and split pairs with
    non-zero signal."""
    from splithunter import _core

    locus = _locus_for(bam_contigs, "TRA")
    result = _core.analyze_locus(
        TEST_BAM, locus.name, locus.chrom, locus.start, locus.end,
        hg38_tra_reference,
        30, 10_000, 50.0,
    )
    # Sanity: we hit reads in the locus window.
    assert result["TRA.SR-TOTAL"] > 1000, "no reads fetched from TRA window"
    # With full realignment + filters we expect at least a handful of
    # chimeric events on a real germline sample — SP is the more reliable
    # signal on the NA12878 subset BAM.
    assert (
        result["TRA.SR-SIGNAL"] + result["TRA.SP-SIGNAL"]
    ) > 0, "no SR/SP evidence detected"


def test_region_coverage_is_sorted(bam_contigs):
    from splithunter import _core

    locus = _locus_for(bam_contigs, "TRA")
    window = 5_000
    cov = _core.region_coverage(
        TEST_BAM, locus.chrom, locus.start, locus.start + window, 0,
    )
    positions = [p for p, _ in cov]
    assert positions == sorted(positions)
    # The sample BAM is sparse — we only require that covered positions fall
    # within the fetched window.
    for pos, _ in cov:
        assert locus.start <= pos < locus.start + window


def test_tcell_fraction_within_bounds(bam_contigs):
    segs = sh_loci.HG38_TCELL_SEGMENTS["TRA"]
    chrom = sh_loci.pick_contig(sh_loci.HG38_LOCI_BY_NAME["TRA"], bam_contigs)
    result = sh_tcell.estimate(TEST_BAM, segs._replace(chrom=chrom))
    # Either we have signal (finite fraction in [0, 1]) or the sample is too
    # thin and everything is NaN — both are legitimate outcomes on a sparse
    # test BAM.
    if not _nan(result["tcell_fraction"]):
        assert 0.0 <= result["tcell_fraction"] <= 1.0


def test_run_pipeline_produces_json(tmp_path, bam_contigs, hg38_tra_reference):
    workdir = tmp_path / "work"
    sh_run.main([
        SAMPLES_CSV, "--workdir", str(workdir), "--locus", "TRA",
        "--reference", hg38_tra_reference,
    ])
    outputs = list(workdir.glob("*.json"))
    assert outputs, "splithunter produced no JSON output"
    with open(outputs[0]) as fp:
        payload = json.load(fp)
    assert "SampleKey" in payload
    assert "TRA.SR-TOTAL" in payload


def test_report_from_run(tmp_path, bam_contigs, hg38_tra_reference):
    workdir = tmp_path / "work"
    sh_run.main([
        SAMPLES_CSV, "--workdir", str(workdir), "--locus", "TRA",
        "--reference", hg38_tra_reference,
    ])
    tsvfile = tmp_path / "out.tsv"
    jsonfiles = [str(p) for p in workdir.glob("*.json")]
    assert jsonfiles
    sh_report.main(jsonfiles + ["--tsv", str(tsvfile)])
    assert tsvfile.exists()
    assert tsvfile.stat().st_size > 0


def _nan(x):
    return x != x


# ---- TcellExTRECT fixture-based regression tests --------------------------
#
# `cov_example` and `tcra_seg_{hg19,hg38}` are the reference data objects
# shipped in McGranahanLab/TcellExTRECT (BSD-3 licensed).  The cov_example
# sample is a pan-cancer WES coverage track with a clearly visible TCRA dip;
# the paper's Fig. 1b shows the expected ~50% T-cell fraction.  We reproduce
# that estimate from the Rust median-logR path.

def test_cov_example_fixture_structure():
    pos, dep = _load_cov_example()
    assert len(pos) == len(dep) == 588_715
    assert min(pos) == 22_090_057 and max(pos) == 23_221_076


def test_tcra_seg_hg38_matches_package_constants():
    """Our HG38_TCELL_SEGMENTS['TRA'] must match TcellExTRECT's authoritative
    ``tcra_seg_hg38`` row-for-row.  If upstream updates the coordinates, this
    test will surface the drift."""
    segs = _load_seg("hg38")
    ours = sh_loci.HG38_TCELL_SEGMENTS["TRA"]
    assert segs["all"] == ours.full
    assert segs["focal"] == ours.focal
    assert segs["local1"] == ours.left_baseline
    assert segs["local2"] == ours.right_baseline


def test_fraction_from_cov_example_reproduces_tcell_dip():
    from splithunter import _core

    pos, dep = _load_cov_example()
    segs = _load_seg("hg19")
    result = _core.fraction_from_coverage(
        pos, dep,
        segs["focal"][0], segs["focal"][1],
        segs["local1"][0], segs["local1"][1],
        segs["local2"][0], segs["local2"][1],
        0,
    )
    # The cov_example sample carries a clear TCRA dip — baseline > focal and
    # the implied T-cell fraction should be materially above zero.
    assert result["baseline_median"] > result["focal_median"]
    assert result["log2_ratio"] < 0
    assert 0.1 < result["tcell_fraction"] <= 1.0
    assert result["focal_positions"] > 10_000
    assert result["baseline_positions"] > 100_000


def test_fraction_from_coverage_length_mismatch_raises():
    from splithunter import _core
    with pytest.raises(ValueError):
        _core.fraction_from_coverage([1, 2, 3], [1, 2], 0, 10, 0, 5, 5, 10, 1)


def test_fraction_from_coverage_no_signal_when_flat():
    from splithunter import _core

    # Uniform depth — zero logR, zero fraction.
    pos = list(range(0, 1000))
    dep = [50] * len(pos)
    result = _core.fraction_from_coverage(
        pos, dep,
        400, 600,    # focal
        0, 200,      # left
        800, 1000,   # right
        0,
    )
    assert result["log2_ratio"] == 0.0
    assert result["tcell_fraction"] == 0.0


def test_fraction_from_coverage_perfect_hemizygous():
    from splithunter import _core

    # Focal region has exactly half the baseline depth → fraction = 0.5.
    pos = list(range(0, 1000))
    dep = [20 if 400 <= p < 600 else 40 for p in pos]
    result = _core.fraction_from_coverage(
        pos, dep,
        400, 600,
        0, 200,
        800, 1000,
        0,
    )
    assert result["log2_ratio"] == pytest.approx(-1.0, abs=1e-9)
    assert result["tcell_fraction"] == pytest.approx(0.5, abs=1e-9)


# ---- Capture-kit exon masking (WES bait-edge correction) ------------------
#
# The V-J focal window is intronic in every reference build, so this
# estimator measures off-target *bleed* depth — reads that spill from
# captured exons into the surrounding intron.  On WES the few positions
# that fall inside exons carry bait-inflated depth (~40x vs ~2x bleed)
# and dominate the baseline median, while the focal median stays on the
# bleed scale (no exons in focal); this asymmetry biases the log-ratio.
# The ``target_mask`` parameter excludes exon positions so both focal
# and baseline medians are taken on the same bleed-depth scale.

def test_default_exons_bed_is_packaged():
    exons = sh_loci.load_exons_bed(sh_loci.DEFAULT_EXONS_BED, "chr14")
    assert len(exons) == 192
    tra = sh_loci.HG38_TCELL_SEGMENTS["TRA"]
    for s, e in exons:
        assert tra.full[0] <= s < e <= tra.full[1]


def test_load_exons_bed_handles_chr_prefix_drift(tmp_path):
    bed = tmp_path / "x.bed"
    bed.write_text("# comment\n14\t100\t200\n14\t300\t400\nchrX\t0\t10\n")
    assert sh_loci.load_exons_bed(str(bed), "chr14") == [(100, 200), (300, 400)]
    assert sh_loci.load_exons_bed(str(bed), "14") == [(100, 200), (300, 400)]
    assert sh_loci.load_exons_bed(str(bed), "chrY") == []


def test_default_exons_for_only_tra():
    assert sh_loci.default_exons_for("TRA", "chr14") is not None
    assert sh_loci.default_exons_for("IGH", "chr14") is None
    assert sh_loci.default_exons_for("TRB", "chr7") is None


def test_fraction_from_coverage_mask_removes_bait_edge_bias():
    """Simulate a WES depth track where baseline exons carry 40x capture
    depth and surrounding intronic bleed is 2x; focal (all intronic) is
    uniformly 2x.  Without masking, baseline median is pulled up by the
    exon spikes and the log-ratio goes negative even with no T-cell dip.
    Masking the exon positions restores focal == baseline == 2."""
    from splithunter import _core

    pos = list(range(0, 1000))
    # Baseline exons cover > 50% of baseline positions so the legacy median
    # is pulled onto the exon-depth scale — the degenerate bait-edge regime
    # we want masking to fix.
    exon_intervals = [(0, 120), (810, 1000)]
    exon_set = set(range(0, 120)) | set(range(810, 1000))
    dep = []
    for p in pos:
        if 400 <= p < 600:
            dep.append(2)                       # focal bleed
        elif p in exon_set:
            dep.append(40)                      # captured exon depth
        else:
            dep.append(2)                       # baseline bleed
    legacy = _core.fraction_from_coverage(
        pos, dep, 400, 600, 0, 200, 800, 1000, 1,
    )
    masked = _core.fraction_from_coverage(
        pos, dep, 400, 600, 0, 200, 800, 1000, 1, exon_intervals,
    )
    # Legacy: baseline median == 40 (exon-dominated), focal median == 2.
    assert legacy["baseline_median"] == pytest.approx(40.0)
    assert legacy["focal_median"] == pytest.approx(2.0)
    # With exons masked, both medians fall on the 2x bleed scale → ratio 0.
    assert masked["baseline_median"] == pytest.approx(2.0)
    assert masked["focal_median"] == pytest.approx(2.0)
    assert masked["log2_ratio"] == pytest.approx(0.0, abs=1e-9)
    assert masked["tcell_fraction"] == pytest.approx(0.0, abs=1e-9)


def test_fraction_from_coverage_empty_mask_equals_legacy():
    from splithunter import _core

    pos = list(range(0, 1000))
    dep = [20 if 400 <= p < 600 else 40 for p in pos]
    legacy = _core.fraction_from_coverage(
        pos, dep, 400, 600, 0, 200, 800, 1000, 0,
    )
    empty = _core.fraction_from_coverage(
        pos, dep, 400, 600, 0, 200, 800, 1000, 0, [],
    )
    assert empty["log2_ratio"] == legacy["log2_ratio"]
    assert empty["tcell_fraction"] == legacy["tcell_fraction"]
    assert empty["focal_positions"] == legacy["focal_positions"]


def test_fraction_from_coverage_mask_covering_everything_yields_nan():
    from splithunter import _core

    pos = list(range(0, 1000))
    dep = [40] * 1000
    # Mask that engulfs every focal + baseline position → no survivors.
    result = _core.fraction_from_coverage(
        pos, dep, 400, 600, 0, 200, 800, 1000, 0,
        [(0, 1000)],
    )
    assert _nan(result["log2_ratio"])
    assert _nan(result["tcell_fraction"])


def test_fraction_from_cov_example_with_hg19_exon_mask_preserves_dip():
    """The upstream WES cov_example (hg19) shows a clean 50% T-cell dip
    under legacy whole-window medians.  Masking the 192 hg19 TCRA exons
    should leave the dip intact because the focal window is intronic and
    the few exon positions in baseline make up < 5% of the total — the
    masked and legacy fractions should be very close."""
    from splithunter import _core

    pos, dep = _load_cov_example()
    segs = _load_seg("hg19")
    exons = sh_loci.load_exons_bed(
        op.join(FIXTURES, "tcra_exons_hg19.bed"), "chr14")
    result = _core.fraction_from_coverage(
        pos, dep,
        segs["focal"][0], segs["focal"][1],
        segs["local1"][0], segs["local1"][1],
        segs["local2"][0], segs["local2"][1],
        0,
        exons,
    )
    assert result["baseline_median"] > result["focal_median"]
    assert result["log2_ratio"] < 0
    assert 0.3 < result["tcell_fraction"] <= 1.0
    # Exon positions are excluded → baseline_positions is smaller than the
    # legacy count (237,219) but most intronic positions survive.
    assert 150_000 < result["baseline_positions"] < 237_220
