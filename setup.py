#!/usr/bin/env python

import os.path as op

from setuptools import setup
from setup_helper import SetupHelper


name = "splithunter"
classifiers = [
    'Development Status :: 4 - Beta',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: BSD License',
    'Programming Language :: Python',
    'Programming Language :: Python :: 2',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
]
with open('requirements.txt') as f:
    required = f.read().splitlines()

# Use the helper
h = SetupHelper(initfile=op.join(name, "__init__.py"), readmefile="README.md")
h.check_version(name, majorv=2, minorv=7)

setup(
      name=name,
      version=h.version,
      author=h.author,
      author_email=h.email,
      license=h.license,
      long_description=h.long_description,
      packages=[name],
      package_dir={name : name},
      include_package_data=True,
      package_data={name: ["data/*.*"]},
      scripts=["splithunter.py", "splithunter_report.py"],
      classifiers=classifiers,
      zip_safe=False,
      url='https://github.com/tanghaibao/splithunter',
      description='Split read hunter',
      install_requires=required + ['pysam']
)
