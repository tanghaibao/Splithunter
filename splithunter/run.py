#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""
Splithunter: parse target V(D)J regions from a BAM file, detect split reads
and split read pairs, and optionally compute a TcellExTRECT-style coverage
ratio.  Core heavy-lifting runs in Rust via ``splithunter._core``.
"""

import argparse
import json
import logging
import os
import os.path as op
import sys
import time

from datetime import timedelta
from multiprocessing import Pool, cpu_count

import pandas as pd

from . import _core
from .loci import (
    HG38_LOCI,
    HG38_LOCI_BY_NAME,
    HG38_TCELL_SEGMENTS,
    default_exons_for,
    load_exons_bed,
    pick_contig,
)
from .utils import DefaultHelpParser, get_abs_path, mkdir

logging.basicConfig()
logger = logging.getLogger(__name__)

PACKAGE_ROOT = op.dirname(op.abspath(__file__))
DATA_DIR = get_abs_path(op.join(PACKAGE_ROOT, "data"))
HLI_BAMS = op.join(DATA_DIR, "HLI_bams.csv.gz")
LOCI = tuple(l.name for l in HG38_LOCI)


def set_argparse():
    p = DefaultHelpParser(
        description=__doc__,
        prog="splithunter_run",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument('infile', nargs='?',
                   help="Input path (BAM, list of BAMs, or csv format)")
    p.add_argument('--locus', default="", choices=LOCI,
                   help="Compute only one locus (default: all)")
    p.add_argument('--cpus',
                   help='Number of CPUs to use', type=int, default=cpu_count())
    p.add_argument('--log', choices=("INFO", "DEBUG"), default="INFO",
                   help='Print debug logs, DEBUG=verbose')
    p.add_argument('--tcell-fraction', action='store_true',
                   help='Also compute TcellExTRECT-style coverage fraction')
    p.add_argument('--targets', default=None,
                   help='Capture-kit exon BED. Positions inside listed '
                        'intervals are EXCLUDED from the coverage median so '
                        'that bait-inflated depth does not skew the baseline '
                        'while the focal (intronic) V-J window is untouched. '
                        'Recommended for WES BAMs.  TRA falls back to a '
                        'packaged default exon set when --targets is omitted.')
    p.add_argument('--no-default-targets', action='store_true',
                   help='Disable the packaged default exon mask — compute '
                        'over full focal/baseline windows (WGS mode).')
    p.add_argument('--pad', type=int, default=30,
                   help='Significant match / clip padding')
    p.add_argument('--indel', type=int, default=10_000,
                   help='Minimum split distance (bp)')
    p.add_argument('--min-entropy', type=float, default=50.0,
                   help='Minimum trinucleotide entropy')
    p.add_argument('--reference', required=True,
                   help='Indexed FASTA reference (.fai alongside) — used to '
                        'build a per-locus in-memory BWA-MEM index for '
                        're-aligning soft-clipped reads')
    set_aws_opts(p)
    return p


def set_aws_opts(p):
    group = p.add_argument_group("AWS and Docker options")
    group.add_argument("--sample_id", help="Sample ID")
    group.add_argument("--workflow_execution_id", help="Workflow execution ID")
    group.add_argument("--input_bam_path", help="Input s3 path, override infile")
    group.add_argument("--output_path", help="Output s3 path")
    group.add_argument("--workdir", default=os.getcwd(), help="Specify work dir")


def bam_path(bam):
    """Return network URLs untouched; resolve local paths to absolute."""
    if bam.startswith(("s3://", "http://", "https://", "ftp://")):
        return bam
    return op.abspath(bam)


def _bam_contigs(path):
    import pysam
    with pysam.AlignmentFile(path) as bf:
        return set(bf.references)


def analyze_bam(samplekey, bam, args, loci):
    """Run the Rust core over each requested locus and return a result dict."""
    result = {"SampleKey": samplekey, "bam": bam}
    try:
        contigs = _bam_contigs(bam)
    except Exception as e:
        logger.error("Cannot open BAM %s: %s", bam, e)
        return None

    for locus in loci:
        try:
            chrom = pick_contig(locus, contigs)
        except KeyError as e:
            logger.warning("%s: %s", samplekey, e)
            continue
        logger.debug("Analyzing %s @ %s:%d-%d", locus.name, chrom, locus.start, locus.end)
        try:
            locus_result = _core.analyze_locus(
                bam, locus.name, chrom, locus.start, locus.end,
                args.reference,
                args.pad, args.indel, args.min_entropy,
            )
        except Exception as e:
            logger.error("%s %s failed: %s", samplekey, locus.name, e)
            return None
        result.update(locus_result)

        if args.tcell_fraction and locus.name in HG38_TCELL_SEGMENTS:
            seg = HG38_TCELL_SEGMENTS[locus.name]._replace(chrom=chrom)
            if args.targets:
                mask = load_exons_bed(args.targets, chrom)
            elif args.no_default_targets:
                mask = None
            else:
                mask = default_exons_for(locus.name, chrom)
            try:
                tc = _core.tcell_fraction(
                    bam, seg.chrom,
                    seg.focal[0], seg.focal[1],
                    seg.left_baseline[0], seg.left_baseline[1],
                    seg.right_baseline[0], seg.right_baseline[1],
                    1, 1,
                    mask,
                )
                for k, v in dict(tc).items():
                    result[f"{locus.name}.TCELL-{k}"] = v
                result[f"{locus.name}.TCELL-mask_intervals"] = (
                    0 if not mask else len(mask)
                )
            except Exception as e:
                logger.warning("%s %s fraction failed: %s", samplekey, locus.name, e)

    return result


def run(arg):
    samplekey, bam, args = arg
    loci = HG38_LOCI
    if args.locus:
        loci = [HG38_LOCI_BY_NAME[args.locus]]
    return analyze_bam(samplekey, bam, args, loci)


def to_json(results):
    if not results:
        return
    samplekey = results['SampleKey']
    bam = results['bam']
    logger.debug("Writing JSON for %s `%s`", samplekey, bam)
    jsonfile = ".".join((samplekey, "json"))
    with open(jsonfile, "w") as fw:
        json.dump(results, fw, sort_keys=True, indent=4,
                  separators=(',', ': '), default=float)
        fw.write("\n")


def get_HLI_bam(samplekey):
    """From @176449128, retrieve the S3 path of the BAM."""
    df = pd.read_csv(HLI_BAMS, index_col=0, dtype=str, header=None,
                     names=["SampleKey", "BAM"])
    return df.loc[samplekey]["BAM"]


def read_csv(csvfile, args):
    # Mode 0: HLI id, starting with @
    if csvfile[0] == '@':
        samplekey = csvfile[1:]
        bam = get_HLI_bam(samplekey)
        return [(samplekey, bam)]

    # Mode 1: single BAM/CRAM file
    if csvfile.endswith((".bam", ".cram")):
        bam = bam_path(csvfile)
        if args.workflow_execution_id and args.sample_id:
            samplekey = "_".join((args.workflow_execution_id, args.sample_id))
        else:
            samplekey = op.basename(bam).rsplit(".", 1)[0]
        return [(samplekey, bam)]

    contents = []
    with open(csvfile) as fp:
        header = fp.readline().strip()

        # Mode 2: just a list of BAM files
        if header.endswith(".bam") and header.count(",") == 0:
            fp.seek(0)
            for row in fp:
                bam = bam_path(row.strip())
                samplekey = op.basename(bam).rsplit(".", 1)[0]
                contents.append((samplekey, bam))
            return contents

        # Mode 3: CSV with samplekey,bam
        fp.seek(0)
        for row in fp:
            atoms = row.strip().split(",")
            if len(atoms) < 2:
                continue
            samplekey, bam = atoms[:2]
            bam = bam_path(bam)
            if bam.endswith(".bam"):
                contents.append((samplekey, bam))
    return contents


def main(args=None):
    p = set_argparse()
    args = p.parse_args(args)

    loglevel = getattr(logging, args.log.upper(), logging.INFO)
    logger.setLevel(loglevel)
    logger.debug("Commandline Arguments: %s", vars(args))

    start = time.time()
    workdir = args.workdir
    cwd = os.getcwd()
    if workdir != cwd:
        mkdir(workdir, logger=logger)

    infile = args.input_bam_path or args.infile
    if not infile:
        sys.exit(not p.print_help())

    samples = read_csv(infile, args)
    logger.debug("Total samples: %d", len(samples))

    task_args = []
    os.chdir(workdir)
    for samplekey, bam in samples:
        jsonfile = ".".join((samplekey, "json"))
        if op.exists(jsonfile):
            continue
        task_args.append((samplekey, bam, args))

    cpus = min(args.cpus, len(task_args)) if task_args else 0
    if cpus == 0:
        logger.debug("All jobs already completed.")
    else:
        logger.debug("Starting %d threads for %d jobs.", cpus, len(task_args))
        if cpus == 1:
            for ta in task_args:
                to_json(run(ta))
        else:
            with Pool(processes=cpus) as pool:
                for results in pool.imap(run, task_args):
                    to_json(results)

    print(f"Elapsed time={timedelta(seconds=time.time() - start)}",
          file=sys.stderr)
    os.chdir(cwd)


if __name__ == '__main__':
    main(sys.argv[1:])
