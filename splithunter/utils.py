#!/usr/bin/env python3
# -*- coding: UTF-8 -*-


import argparse
import os
import os.path as op
import shutil
import sys


class DefaultHelpParser(argparse.ArgumentParser):

    def error(self, message):
        sys.stderr.write('error: {}\n\n'.format(message))
        sys.exit(not self.print_help())


def get_abs_path(path):
    return op.realpath(path)


def mkdir(dirname, overwrite=False, logger=None):
    """
    Wraps around os.makedirs(), but checks for existence first.
    """
    if op.isdir(dirname):
        if overwrite:
            shutil.rmtree(dirname)
            os.makedirs(dirname)
            if logger:
                logger.debug("Overwrite folder `{}`.".format(dirname))
        else:
            return False
    else:
        os.makedirs(dirname, exist_ok=True)
        if logger:
            logger.debug("`{}` not found. Creating new.".format(dirname))

    return True


def is_exe(fpath):
    """ Determines if a file is executable
    """
    return op.isfile(fpath) and os.access(fpath, os.X_OK)


def which(program, paths=None):
    """
    Emulates the unix which command.
    """
    paths = paths or []
    fpath, _ = op.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep) + list(paths):
            exe_file = op.join(path, program)
            if is_exe(exe_file):
                return op.abspath(exe_file)

    return None
