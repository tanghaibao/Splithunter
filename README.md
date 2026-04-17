# Split read hunter

[![CI](https://github.com/tanghaibao/Splithunter/actions/workflows/ci.yml/badge.svg)](https://github.com/tanghaibao/Splithunter/actions/workflows/ci.yml)

|||
|---|---|
| Author | Haibao Tang ([tanghaibao](http://github.com/tanghaibao)) |
| Email | <tanghaibao@gmail.com> |
| License | BSD-3-Clause |

## Description

Detect V(D)J split reads, split read pairs, and TcellExTRECT-style coverage
dips directly from a BAM file.  The performance-sensitive core is implemented
in Rust (via `rust-htslib`) and exposed to Python as a PyO3 extension.

## Installation

Requires Rust (stable) and Python 3.9+.  Build via [maturin](https://www.maturin.rs):

```bash
pip install maturin
maturin develop --release          # editable install for development
# or
pip install .                      # regular install (also invokes maturin)
```

## Usage

### Split-read / split-pair detection

```bash
splithunter_run HLI_bams.csv --workdir hli --locus TRA
splithunter_report hli/*.json --tsv hli.splithunter.tsv
```

Add `--tcell-fraction` to also emit TcellExTRECT-style coverage ratios per
locus in each JSON.

### Standalone T-cell fraction (TcellExTRECT port)

```bash
splithunter_tcell path/to/sample.bam --locus TRA
# or:
splithunter_tcell path/to/sample.bam --locus TRA --json
```

This computes the median log2 ratio between the focal V-J window and the
flanking baseline windows and derives the fraction as `1 - 2^logR` (clamped to
[0, 1]).

## Development

```bash
pip install maturin "pytest" pysam pandas
maturin develop --release
pytest tests/
cd rust && cargo test --release
```
