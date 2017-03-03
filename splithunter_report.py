#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Compile splithunter calls into a tsv file from a set of JSON files.
"""

import argparse
import json
import sys
import pandas as pd

from multiprocessing import Pool, cpu_count

from pysrc.utils import DefaultHelpParser


def df_to_tsv(df, tsvfile):
    dd = ["samplekey"]
    columns = dd + sorted([x for x in df.columns if (x not in dd)])

    df = df.reindex_axis(columns, axis='columns')
    df.sort_values("samplekey")
    df.to_csv(tsvfile, sep='\t', index=False)
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


def main():
    p = DefaultHelpParser(description=__doc__, prog=__file__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("files", nargs="*")
    p.add_argument('--tsv', default="out.tsv",
                   help="Path to the tsv file")
    p.add_argument('--cpus', default=cpu_count(),
                   help='Number of threads')
    args = p.parse_args()

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
        df = pd.read_csv(tsvfile, sep="\t")

    if df.empty:
        print >> sys.stderr, "Dataframe empty - check input files"
        sys.exit(1)

if __name__ == '__main__':
    main()
