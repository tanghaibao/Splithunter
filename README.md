# Split read hunter

[![CI](https://github.com/tanghaibao/Splithunter/actions/workflows/ci.yml/badge.svg)](https://github.com/tanghaibao/Splithunter/actions/workflows/ci.yml)

|||
|---|---|
| Author | Haibao Tang ([tanghaibao](http://github.com/tanghaibao)) |
| Email | <tanghaibao@gmail.com> |
| License | BSD-3-Clause |

## Description

Identify split reads and read pairs in a particular region.

## Installation

```bash
git clone --recursive https://github.com/tanghaibao/splithunter
cd splithunter
cd src && make -j && cd ..
pip install .
```

Requires Python 3.8+.

## Usage

- Run batch jobs listed in ``HLI_bams.csv``:

```bash
splithunter_run HLI_bams.csv --workdir hli --locus TRA
```

- Collect results into a TSV file:

```bash
splithunter_report hli/*.json --tsv hli.splithunter.tsv
```

## Build custom database

The default regions include `TRA`, `TRB`, `TRG`, `IGH`, `IGK`, `IGL`. To build
BWA indices for other regions, use ``BuildDB``:

```bash
src/BuildDB src/data/TR_IG.bed -r ~/projects/ref/hg38.upper.fa
```

## Development

```bash
pip install -e ".[test]"
pytest tests/
```
