#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
HLI sjTREC caller: parse the target regions from the input BAM file, extract
split reads and split read pairs to test whether they are associated with aging.
"""

import argparse
import shutil
import os
import os.path as op
import json
import sys
import time
import logging

import pandas as pd

from .utils import DefaultHelpParser, mkdir, get_abs_path, which
from datetime import timedelta
from multiprocessing import Pool, cpu_count

logging.basicConfig()
logger = logging.getLogger(__name__)

datadir = get_abs_path(op.join(op.dirname(__file__) + "/..", 'data'))
datafile = lambda x: op.join(datadir, x)
HLI_BAMS = datafile("HLI_bams.csv.gz")
LOCI = "TRA|TRB|TRG|IGH|IGK|IGL"
exec_path = which("Splithunter", paths=[op.dirname(__file__) + "/.."])


def set_argparse():
    p = DefaultHelpParser(description=__doc__, prog=__file__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('infile', nargs='?',
                    help="Input path (BAM, list of BAMs, or csv format)")
    p.add_argument('--locus', default="", choices=LOCI.split("|"),
                    help="Compute only one locus, must be one of {}".format(LOCI))
    p.add_argument('--cpus',
                    help='Number of CPUs to use', type=int, default=cpu_count())
    p.add_argument('--log', choices=("INFO", "DEBUG"), default="INFO",
                    help='Print debug logs, DEBUG=verbose')
    set_aws_opts(p)
    return p


def set_aws_opts(p):
    group = p.add_argument_group("AWS and Docker options")
    p.add_argument_group(group)
    # https://github.com/hlids/infrastructure/wiki/Docker-calling-convention
    group.add_argument("--sample_id", help="Sample ID")
    group.add_argument("--workflow_execution_id", help="Workflow execution ID")
    group.add_argument("--input_bam_path", help="Input s3 path, override infile")
    group.add_argument("--output_path", help="Output s3 path")
    group.add_argument("--workdir", default=os.getcwd(), help="Specify work dir")


def bam_path(bam):
    '''
    Is the file network-based?
    '''
    if bam.startswith("s3://") or bam.startswith("http://") or \
       bam.startswith("ftp://") or bam.startswith("https://"):
        return bam
    else:
        return op.abspath(bam)


def check_bam(bam):
    # Check indices - remove if found, otherwise distinct bams may try to use
    # the previous index, leading to an error
    bbam = op.basename(bam)
    baifile = bbam + ".bai"
    altbaifile = bbam.rsplit(".", 1)[0] + ".bai"

    for bai in (baifile, altbaifile):
        if op.exists(bai):
            logger.debug("Remove index `{}`".format(bai))
            os.remove(bai)

    # Does the file exist?
    logger.debug("Working on `{}`".format(bam))

    return bam


def cleanup(cwd, samplekey):
    """
    Change back to the parent folder and remove the samplekey folder after done
    """
    os.chdir(cwd)
    shutil.rmtree(samplekey)


def run(arg):
    '''
    :param bam: path to bam file
    :return: dict of calls
    '''
    samplekey, bam, args = arg
    if not op.exists(exec_path):
        logging.error("{} not found. Please check or recompile."\
                    .format(exec_path))
        sys.exit(1)

    cwd = os.getcwd()
    mkdir(samplekey)
    os.chdir(samplekey)

    res = { 'SampleKey': samplekey, 'bam': bam }
    failed = False
    if check_bam(bam) is None:
        cleanup(cwd, samplekey)
        return res

    cmd = "{} {} -s {}".format(exec_path, bam, samplekey)
    if args.locus:
        cmd += " -l {}".format(args.locus)
    try:
        print >> sys.stderr, cmd
        os.system(cmd)
        jsonfile = "{}.json".format(samplekey)

        # Read JSON from program output
        fp = open(jsonfile)
        res = json.load(fp)
        fp.close()
    except Exception as e:
        logger.error("Exception on `{}` {} ({})".format(bam, samplekey, e))
        failed = True

    cleanup(cwd, samplekey)
    return None if failed else res


def to_json(results):
    if results is None:
        return

    samplekey = results['SampleKey']
    if not results:
        logger.debug("No calls are found for {} `{}`".format(samplekey, bam))
        return

    jsonfile = ".".join((samplekey, "json"))
    js = json.dumps(results, sort_keys=True, indent=4, separators=(',', ': '))

    # Write JSON with Python
    fw = open(jsonfile, "w")
    print >> fw, js
    fw.close()


def get_HLI_bam(samplekey):
    """
    From @176449128, retrieve the S3 path of the BAM
    """
    df = pd.read_csv(HLI_BAMS, index_col=0, dtype=str, header=None,
                     names=["SampleKey", "BAM"])
    return df.ix[samplekey]["BAM"]


def read_csv(csvfile, args):
    # Mode 0: See if this is just a HLI id, starting with @
    if csvfile[0] == '@':
        samplekey = csvfile[1:]
        bam = get_HLI_bam(samplekey)
        return [(samplekey, bam)]

    # Mode 1: See if this is just a BAM file
    if csvfile.endswith(".bam") or csvfile.endswith(".cram"):
        bam = csvfile
        bam = bam_path(bam)
        if args.workflow_execution_id and args.sample_id:
            samplekey = "_".join((args.workflow_execution_id, args.sample_id))
        else:
            samplekey = op.basename(bam).rsplit(".", 1)[0]
        return [(samplekey, bam)]

    fp = open(csvfile)
    # Mode 2: See if the file contains JUST list of BAM files
    header = fp.next().strip()
    contents = []
    if header.endswith(".bam") and header.count(",") == 0:
        fp.seek(0)
        for row in fp:
            bam = row.strip()
            bam = bam_path(bam)
            samplekey = op.basename(bam).rsplit(".", 1)[0]
            contents.append((samplekey, bam))
        return contents

    # Mode 3: Continue reading, this is a CSV file
    fp.seek(0)
    for row in fp:
        atoms = row.strip().split(",")
        samplekey, bam = atoms[:2]
        bam = bam_path(bam)
        if bam.endswith(".bam"):
            contents.append((samplekey, bam))
    return contents


def main(args):
    p = set_argparse()
    args = p.parse_args(args)

    loglevel = getattr(logging, args.log.upper(), "INFO")
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

    # Parallel processing
    for i, (samplekey, bam) in enumerate(samples):
        jsonfile = ".".join((samplekey, "json"))
        if op.exists(jsonfile):
            continue
        task_args.append((samplekey, bam, args))

    cpus = min(args.cpus, len(task_args))
    if cpus == 0:
        logger.debug("All jobs already completed.")

    else:
        logger.debug("Starting {} threads for {} jobs.".format(cpus, len(task_args)))

        if cpus == 1:  # Serial
            for ta in task_args:
                results = run(ta)
                to_json(results)
        else:
            p = Pool(processes=cpus)
            for results in p.imap(run, task_args):
                to_json(results)

    print >> sys.stderr, "Elapsed time={}"\
            .format(timedelta(seconds=time.time() - start))
    os.chdir(cwd)


if __name__ == '__main__':
    main(sys.argv[1:])
