#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Copyright (c) 2015-2017 Human Longevity Inc.

Author: Haibao Tang <htang@humanlongevity.com>
License: Non-Commercial Use Only. For details, see `LICENSE` file

HLI TRA deletion caller: parse the TRA regions from the input BAM file, extract
deletion events.
"""

import sys
from splithunter.run import main


if __name__ == '__main__':
    main(sys.argv[1:])
