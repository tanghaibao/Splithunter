"""
Immunoglobulin and T-cell receptor locus coordinates (hg38).

Two parallel maps are provided:

* ``HG38_LOCI`` - whole-locus windows used by the split-read hunter; define
  the region fetched from the BAM for each V(D)J locus.
* ``HG38_TCELL_SEGMENTS`` - per-locus focal V-J window plus two flanking
  baseline windows, following the TcellExTRECT (Bentham et al., Nature 2021)
  parameterisation extended heuristically to the IG loci.  The coordinates
  are 0-based half-open (BED-style).
"""

from collections import namedtuple

Locus = namedtuple("Locus", ["name", "chrom", "start", "end"])

# Focal + flanking windows for the coverage-dip fraction estimator.
LocusSegments = namedtuple(
    "LocusSegments",
    ["name", "chrom", "full", "focal", "left_baseline", "right_baseline"],
)


HG38_LOCI = [
    Locus("TRA", "chr14", 21621904, 22752132),
    Locus("TRB", "chr7", 142299011, 142813287),
    Locus("TRG", "chr7", 38240024, 38368055),
    Locus("IGH", "chr14", 105586437, 106879844),
    Locus("IGK", "chr2", 88857361, 90235368),
    Locus("IGL", "chr22", 22026076, 22922913),
]

HG38_LOCI_BY_NAME = {l.name: l for l in HG38_LOCI}


# Focal V-J regions (approx.) with flanking baselines drawn from the V-distal
# 5' end and the C-distal 3' end of each locus.  Values are hg38; TRA is
# taken from the TcellExTRECT paper, the rest are derived from IMGT / UCSC
# coordinates of the V-proximal J-distal boundaries.
HG38_TCELL_SEGMENTS = {
    # TRA coordinates match TcellExTRECT's authoritative `tcra_seg_hg38`
    # (data-raw in the McGranahanLab/TcellExTRECT repo).
    "TRA": LocusSegments(
        "TRA", "chr14",
        full=(21_621_904, 22_752_132),
        focal=(22_331_570, 22_411_598),
        left_baseline=(21_621_904, 21_830_070),
        right_baseline=(22_547_506, 22_752_132),
    ),
    "TRB": LocusSegments(
        "TRB", "chr7",
        full=(142_299_011, 142_813_287),
        focal=(142_299_011, 142_510_993),
        left_baseline=(142_299_011, 142_299_500),
        right_baseline=(142_520_000, 142_800_000),
    ),
    "TRG": LocusSegments(
        "TRG", "chr7",
        full=(38_240_024, 38_368_055),
        focal=(38_240_024, 38_368_055),
        left_baseline=(38_230_000, 38_240_000),
        right_baseline=(38_368_055, 38_378_000),
    ),
    "IGH": LocusSegments(
        "IGH", "chr14",
        full=(105_586_437, 106_879_844),
        focal=(105_860_000, 106_879_844),
        left_baseline=(105_586_437, 105_860_000),
        right_baseline=(106_879_844, 106_889_000),
    ),
    "IGK": LocusSegments(
        "IGK", "chr2",
        full=(88_857_361, 90_235_368),
        focal=(88_857_361, 89_330_000),
        left_baseline=(88_840_000, 88_857_361),
        right_baseline=(89_330_000, 90_235_368),
    ),
    "IGL": LocusSegments(
        "IGL", "chr22",
        full=(22_026_076, 22_922_913),
        focal=(22_026_076, 22_380_000),
        left_baseline=(22_020_000, 22_026_076),
        right_baseline=(22_380_000, 22_922_913),
    ),
}


def pick_contig(locus, bam_header_contigs):
    """
    Return the contig name matching the locus as it appears in the BAM header.
    Handles ``chr`` prefix drift between UCSC-style and Ensembl-style BAMs.
    """
    if locus.chrom in bam_header_contigs:
        return locus.chrom
    stripped = locus.chrom.removeprefix("chr")
    if stripped in bam_header_contigs:
        return stripped
    raise KeyError(
        f"Locus {locus.name} contig {locus.chrom!r} not found in BAM "
        "(neither with nor without 'chr')"
    )
