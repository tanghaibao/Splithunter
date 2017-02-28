# sjTREC caller

|||
|---|---|
| Author | Haibao Tang ([tanghaibao](http://github.com/tanghaibao)) |
| Email | <htang@humanlongevity.com> |
| License | Restricted |

## Description

Identify split reads in a particular region.

## Installation

- Checkout a copy of Splithunter and install:

  ```bash
  git clone --recursive https://github.com/tanghaibao/splithunter
  make
  ```

For accessing BAMs that are located on S3, please refer to
`Dockerfiles/sjtrec.dockerfile` for installation of SAMTOOLS/pysam with S3
support.
