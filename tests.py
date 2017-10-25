#!/usr/bin/env python
# -*- coding: UTF-8 -*-


import pytest


def test_splithunter():
    """ Run splithunter.py on sample CSV file
    """
    from splithunter.run import main
    main(["tests/samples.csv", "--workdir", "work", "--locus", "TRA"])


def test_splithunter_report():
    """ Run splithunter_report.py on JSON file
    """
    from splithunter.report import main
    main(["work/NA12878.json", "--tsv", "tests.splithunter.tsv"])
