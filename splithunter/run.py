#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""
HLI sjTREC caller: parse the target regions from the input BAM file, extract
split reads and split read pairs to test whether they are associated with aging.
"""

import argparse
import json
import logging
import os
import os.path as op
import shutil
import subprocess
import sys
import time

from datetime import timedelta
from multiprocessing import Pool, cpu_count

import pandas as pd

from .utils import DefaultHelpParser, get_abs_path, mkdir, which

logging.basicConfig()
logger = logging.getLogger(__name__)

PACKAGE_ROOT = op.dirname(op.abspath(__file__))
REPO_ROOT = op.dirname(PACKAGE_ROOT)
SRC_DIR = op.join(REPO_ROOT, "src")
DATA_DIR = get_abs_path(op.join(SRC_DIR, "data"))
HLI_BAMS = op.join(DATA_DIR, "HLI_bams.csv.gz")
LOCI = ("TRA", "TRB", "TRG", "IGH", "IGK", "IGL")

exec_path = which("Splithunter", paths=[SRC_DIR])


def set_argparse():
    p = DefaultHelpParser(
        description=__doc__,
        prog=__file__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument('infile', nargs='?',
                   help="Input path (BAM, list of BAMs, or csv format)")
    p.add_argument('--locus', default="", choices=LOCI,
                   help="Compute only one locus")
    p.add_argument('--cpus',
                   help='Number of CPUs to use', type=int, default=cpu_count())
    p.add_argument('--log', choices=("INFO", "DEBUG"), default="INFO",
                   help='Print debug logs, DEBUG=verbose')
    set_aws_opts(p)
    return p


def set_aws_opts(p):
    group = p.add_argument_group("AWS and Docker options")
    # https://github.com/hlids/infrastructure/wiki/Docker-calling-convention
    group.add_argument("--sample_id", help="Sample ID")
    group.add_argument("--workflow_execution_id", help="Workflow execution ID")
    group.add_argument("--input_bam_path", help="Input s3 path, override infile")
    group.add_argument("--output_path", help="Output s3 path")
    group.add_argument("--workdir", default=os.getcwd(), help="Specify work dir")


def bam_path(bam):
    """
    Return network URLs untouched; resolve local paths to absolute.
    """
    if bam.startswith(("s3://", "http://", "https://", "ftp://")):
        return bam
    return op.abspath(bam)


def check_bam(bam):
    bbam = op.basename(bam)
    baifile = bbam + ".bai"
    altbaifile = bbam.rsplit(".", 1)[0] + ".bai"

    for bai in (baifile, altbaifile):
        if op.exists(bai):
            logger.debug("Remove index `{}`".format(bai))
            os.remove(bai)

    logger.debug("Working on `{}`".format(bam))
    return bam


def cleanup(cwd, samplekey):
    os.chdir(cwd)
    shutil.rmtree(samplekey, ignore_errors=True)


def run(arg):
    """
    :param arg: (samplekey, bam, args)
    :return: dict of calls or None on failure
    """
    samplekey, bam, args = arg
    if exec_path is None or not op.exists(exec_path):
        logger.error("Splithunter binary not found. Please check or recompile.")
        sys.exit(1)

    cwd = os.getcwd()
    mkdir(samplekey)
    os.chdir(samplekey)

    res = {'SampleKey': samplekey, 'bam': bam}
    if check_bam(bam) is None:
        cleanup(cwd, samplekey)
        return res

    cmd = [exec_path, bam, "-s", samplekey]
    if args.locus:
        cmd += ["-l", args.locus]
    try:
        print(" ".join(cmd), file=sys.stderr)
        subprocess.run(cmd, check=True)
        jsonfile = "{}.json".format(samplekey)
        with open(jsonfile) as fp:
            res = json.load(fp)
    except (subprocess.CalledProcessError, OSError, ValueError) as e:
        logger.error("Exception on `{}` {} ({})".format(bam, samplekey, e))
        cleanup(cwd, samplekey)
        return None

    cleanup(cwd, samplekey)
    return res


def to_json(results):
    if not results:
        return

    samplekey = results['SampleKey']
    bam = results['bam']
    logger.debug("Writing JSON for {} `{}`".format(samplekey, bam))

    jsonfile = ".".join((samplekey, "json"))
    with open(jsonfile, "w") as fw:
        json.dump(results, fw, sort_keys=True, indent=4, separators=(',', ': '))
        fw.write("\n")


def get_HLI_bam(samplekey):
    """
    From @176449128, retrieve the S3 path of the BAM
    """
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
    logger.debug('Commandline Arguments:{}'.format(vars(args)))

    start = time.time()
    workdir = args.workdir
    cwd = os.getcwd()

    if workdir != cwd:
        mkdir(workdir, logger=logger)

    infile = args.input_bam_path or args.infile
    if not infile:
        sys.exit(not p.print_help())

    samples = read_csv(infile, args)
    logger.debug("Total samples: {}".format(len(samples)))

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
        logger.debug("Starting {} threads for {} jobs.".format(cpus, len(task_args)))

        if cpus == 1:
            for ta in task_args:
                to_json(run(ta))
        else:
            with Pool(processes=cpus) as pool:
                for results in pool.imap(run, task_args):
                    to_json(results)

    print("Elapsed time={}".format(timedelta(seconds=time.time() - start)),
          file=sys.stderr)
    os.chdir(cwd)


if __name__ == '__main__':
    main(sys.argv[1:])
