#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

import json
import os
import os.path as op
import shutil

import pytest

from splithunter import run as sh_run
from splithunter import report as sh_report

HERE = op.dirname(op.abspath(__file__))
SAMPLES_CSV = op.join(HERE, "samples.csv")


def _binary_available():
    return sh_run.exec_path is not None and op.exists(sh_run.exec_path)


needs_binary = pytest.mark.skipif(
    not _binary_available(),
    reason="Splithunter C++ binary not built; run `cd src && make`",
)


@needs_binary
def test_splithunter_produces_json(tmp_path):
    workdir = tmp_path / "work"
    sh_run.main([SAMPLES_CSV, "--workdir", str(workdir), "--locus", "TRA"])

    outputs = list(workdir.glob("*.json"))
    assert outputs, "splithunter produced no JSON output"

    with open(outputs[0]) as fp:
        payload = json.load(fp)
    assert "SampleKey" in payload
    assert "TRA.SR-TOTAL" in payload


@needs_binary
def test_splithunter_report_produces_tsv(tmp_path):
    workdir = tmp_path / "work"
    sh_run.main([SAMPLES_CSV, "--workdir", str(workdir), "--locus", "TRA"])

    tsvfile = tmp_path / "out.tsv"
    jsonfiles = [str(p) for p in workdir.glob("*.json")]
    assert jsonfiles, "no JSON inputs for report"
    sh_report.main(jsonfiles + ["--tsv", str(tsvfile)])

    assert tsvfile.exists()
    assert tsvfile.stat().st_size > 0
