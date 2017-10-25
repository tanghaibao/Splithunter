#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Compile splithunter calls into a tsv file from a set of JSON files.
"""

import argparse
import json
import sys
import os.path as op
import pandas as pd

from multiprocessing import Pool, cpu_count

from .utils import DefaultHelpParser


def df_to_tsv(df, tsvfile, index=False):
    dd = ["SampleKey"]
    columns = dd + sorted([x for x in df.columns if (x not in dd)])

    df = df.reindex_axis(columns, axis='columns')
    df.sort_values("SampleKey")
    df.to_csv(tsvfile, sep='\t', index=index)
    print >> sys.stderr, "TSV output written to `{}` (# samples={})"\
                .format(tsvfile, df.shape[0])


def json_to_df_worker(jsonfile):
    return json.load(open(jsonfile))


def json_to_df(jsonfiles, tsvfile, cpus):
    """
    Compile a number of json files into tsv file for easier manipulation.
    """
    df = pd.DataFrame()
    p = Pool(processes=cpus)
    results = []
    r = p.map_async(json_to_df_worker, jsonfiles,
                    callback=results.append)
    r.wait()

    for res in results:
        df = df.append(res, ignore_index=True)
    return df


def get_midpoint(s):
    '''
    >>> get_midpoint("1:58,460-58,557(+)")
    (58508, '+')
    '''
    chr, startendstrand = s.split(':')
    startend, strand = startendstrand.split('(')
    start, end = startend.replace(',', '').split('-')
    start, end = int(start), int(end)
    strand = strand.split(')')[0]
    return (start + end) / 2, strand


def parse_alignments(s):
    '''Convert strings like `1:58,460-58,557(+)|1:135,876-135,929(-);` to tuples for plotting.

    >>> parse_alignments("1:58,460-58,557(+)|1:135,876-135,929(-);1:58,460-58,557(+)|1:135,876-135,929(-);")
    [(58508, 135902, '+-'), (58508, 135902, '+-')]
    '''
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
    # Test slicing rule >800Kb
    sliced = [x for x in details if ((x[0] - offset) > boundary) \
                                 or ((x[1] - offset) > boundary)]
    sliced = [x for x in sliced if x[-1] in orientation]
    if offset == 0: # SR
        sliced = [x for x in sliced if x[0] > x[1]]

    return len(sliced)


def filter_TRA(df, TRA_start=21621904):
    data = {}
    for i, row in df.iterrows():
        samplekey = row["SampleKey"]
        SP, SR = row["TRA.SP-DETAILS"], row["TRA.SR-DETAILS"]
        SP_total, SR_total = row["TRA.SP-TOTAL"], row["TRA.SR-TOTAL"]
        SP_new = SR_new = 0
        if isinstance(SP, basestring):
            SP_details = parse_alignments(SP)
            SP_new = slicing_filter(SP_details, ['-+'], offset=TRA_start)
        if isinstance(SR, basestring):
            SR_details = parse_alignments(SR)
            SR_new = slicing_filter(SR_details, ['++'])

        data[samplekey] = [SR_new * 1e6 / SR_total, SP_new * 1e6 / SP_total,
                           (SR_new + SP_new * 2) * 1e6 / (SR_total + SP_total * 2)]

    xf = pd.DataFrame.from_dict(data, orient="index")
    xf.columns = ["TRA.SR-PPM", "TRA.SP-PPM", "TRA.PPM"]
    return xf


def main(args):
    p = DefaultHelpParser(description=__doc__, prog=__file__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("files", nargs="*")
    p.add_argument('--tsv', default="out.tsv",
                   help="Path to the tsv file")
    p.add_argument('--cpus', default=cpu_count(),
                   help='Number of threads')
    args = p.parse_args(args)

    files = args.files
    tsvfile = args.tsv

    if files:
        nfiles = len(files)
        cpus = min(nfiles, args.cpus)
        suffix = "JSON"
        print >> sys.stderr, "Using {} cpus to parse {} {} files"\
                .format(cpus, nfiles, suffix)
        df = json_to_df(files, tsvfile, cpus)
        df_to_tsv(df, tsvfile)
    else:
        if op.exists(tsvfile):
            df = pd.read_csv(tsvfile, sep="\t")
        else:
            sys.exit(not p.print_help())

    if df.empty:
        print >> sys.stderr, "Dataframe empty - check input files"
        sys.exit(1)

    # Filter TRA locus for age prediction
    print >> sys.stderr, "Filtering TRA locus for age prediction"
    print df
    xf = filter_TRA(df)
    print xf
    tsvfile = tsvfile.rsplit(".", 1)[0] + ".TRA.tsv"
    xf = xf.sort_index()
    xf.to_csv(tsvfile, sep='\t', index_label="SampleKey")
    print >> sys.stderr, "TSV output written to `{}` (# samples={})"\
                .format(tsvfile, xf.shape[0])


if __name__ == '__main__':
    main(sys.argv[1:])
