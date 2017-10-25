#!/usr/bin/env python
# -*- coding: UTF-8 -*-


import argparse
import os
import os.path as op
import shutil
import sys
import logging

from subprocess import PIPE, call


class InputParams:
    '''
    Encapsulates all input parameters to the BamParser class
    '''
    KWARGS_LOG = 'log'

    def __init__(self, bam, **kwargs):
        self.bam = bam
        self.kwargs = kwargs

    def getLogLevel(self, defaultLevel='INFO'):
        '''
        :return: log level specified by 'log' in kwargs, default level if level is not specified or invalid
        '''
        levelName = self.kwargs.get(InputParams.KWARGS_LOG, defaultLevel)
        numericLevel = getattr(logging, levelName.upper(), defaultLevel)
        return numericLevel


class DefaultHelpParser(argparse.ArgumentParser):

    def error(self, message):
        sys.stderr.write('error: {}\n\n'.format(message))
        sys.exit(not self.print_help())


def get_abs_path(link_name):
    source = link_name
    if op.islink(source):
        source = os.readlink(source)
    else:
        source = op.basename(source)

    link_dir = op.dirname(link_name)
    source = op.normpath(op.join(link_dir, source))
    source = op.abspath(source)
    if source == link_name:
        return source
    else:
        return get_abs_path(source)


def mkdir(dirname, overwrite=False, logger=None):
    """
    Wraps around os.mkdir(), but checks for existence first.
    """
    if op.isdir(dirname):
        if overwrite:
            shutil.rmtree(dirname)
            os.mkdir(dirname)
            if logger:
                logger.debug("Overwrite folder `{}`.".format(dirname))
        else:
            return False  # Nothing is changed
    else:
        try:
            os.mkdir(dirname)
        except:
            os.makedirs(dirname)

        if logger:
            logger.debug("`{}` not found. Creating new.".format(dirname))

    return True


def Popen(cmd, stdin=None, stdout=PIPE, debug=False, shell="/bin/bash"):
    """
    Capture the cmd stdout output to a file handle.
    """
    from subprocess import Popen as P
    if debug:
        logging.debug(cmd)
    # See: <https://blog.nelhage.com/2010/02/a-very-subtle-bug/>
    proc = P(cmd, bufsize=1, stdin=stdin, stdout=stdout, \
             shell=True, executable=shell)
    return proc


def popen(cmd, debug=True, shell="/bin/bash"):
    return Popen(cmd, debug=debug, shell=shell).stdout


def s3ify(address):
    if not address.startswith("s3://"):
        address = "s3://" + address.lstrip("/")
    return address


def ls_s3(s3_store_obj_name):
    s3_store_obj_name = s3ify(s3_store_obj_name)
    cmd = "aws s3 ls {0}/".format(s3_store_obj_name)
    contents = []
    for row in popen(cmd):
        contents.append(row.split()[-1])
    return contents


def push_to_s3(s3_store, obj_name):
    cmd = "sync" if op.isdir(obj_name) else "cp"
    s3address = "{0}/{1}".format(s3_store, obj_name)
    s3address = s3ify(s3address)
    cmd = "aws s3 {0} {1} {2} --sse".format(cmd, obj_name, s3address)
    sh(cmd)
    return s3address


def sh(cmd, infile=None, outfile=None, errfile=None,
        append=False, background=False, threaded=None, logger=None,
        shell="/bin/bash"):
    """
    simple wrapper for system calls
    """
    if not cmd:
        return 1

    if infile:
        cat = "cat"
        if infile.endswith(".gz"):
            cat = "zcat"
        cmd = "{0} {1} |".format(cat, infile) + cmd
    if outfile and outfile != "stdout":
        if outfile.endswith(".gz"):
            cmd += " | gzip"
        tag = ">"
        if append:
            tag = ">>"
        cmd += " {0}{1}".format(tag, outfile)
    if errfile:
        if errfile == outfile:
            errfile = "&1"
        cmd += " 2>{0}".format(errfile)
    if background:
        cmd += " &"

    if logger:
        logger.debug(cmd)

    return call(cmd, shell=True, executable=shell)


def is_exe(fpath):
    """ Determines if a file is executable
    """
    return op.isfile(fpath) and os.access(fpath, os.X_OK)


def which(program, paths=[]):
    """
    Emulates the unix which command.
    >>> which("cat")
    "/bin/cat"
    >>> which("nosuchprogram")
    """
    fpath, fname = op.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep) + paths:
            exe_file = op.join(path, program)
            if is_exe(exe_file):
                return op.abspath(exe_file)

    return None
