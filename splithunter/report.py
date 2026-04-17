#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""
Compile splithunter calls into a tsv file from a set of JSON files.
"""

import argparse
import json
import os.path as op
import sys

from multiprocessing import Pool, cpu_count

import pandas as pd

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


def main(args):
    p = DefaultHelpParser(
        description=__doc__,
        prog=__file__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("files", nargs="*")
    p.add_argument('--tsv', default="out.tsv", help="Path to the tsv file")
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


if __name__ == '__main__':
    main(sys.argv[1:])
