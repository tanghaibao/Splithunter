# Split read hunter

[![Travis-CI](https://travis-ci.org/tanghaibao/Splithunter.svg?branch=master)](https://travis-ci.org/tanghaibao/Splithunter)

|||
|---|---|
| Author | Haibao Tang ([tanghaibao](http://github.com/tanghaibao)) |
| Email | <htang@humanlongevity.com> |
| License | Restricted |

## Description

Identify split reads and read pairs in a particular region.

## Installation

- Checkout a copy of Splithunter and install:

```bash
git clone --recursive https://github.com/tanghaibao/splithunter
cd src && make && cd ..
```

## Usage

- Run batch jobs in ``HLI_bams.csv``:

```bash
python splithunter_run.py HLI_bams.csv --workdir hli --locus TRA
```

- Collect results into a TSV file:

```bash
python splithunter_report.py hli/*.json --tsv hli.splithunter.tsv
```

## Build custom database

The default regions includes `TRA`, `TRB`, `TRG`, `IGH`, `IGK`, `IGL`. If
other regions are desired, use ``BuildDB`` to build BWA indices by specifying
the reference genome (e.g. ``hg38.upper.fa``).

- You can build your custom BWA indices using:

```bash
BuildDB src/data/TR_IG.bed -r ~/projects/ref/hg38.upper.fa
```
