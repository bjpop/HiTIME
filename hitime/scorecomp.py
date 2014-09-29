#!/bin/env python

import sys
import os.path
import math
import argparse
import logging
import time

import numpy as np
from rtree import index

# add dir of this (and following) file to path
sys.path.append(os.path.realpath(__file__))
import md_filter
import md_io

#default full width at half max height of signal for retention time. In number of scans.
DEFAULTRTGAP = 3.0  # sec

# default mz range for finding local maxima
DEFAULTMZGAP = 0.5

# command line argument parser
parser = argparse.ArgumentParser(description='Generate list of doublet locations from score data')
parser.add_argument('--rtGap',
                    help='Retention Time gap size for finding local maxima (sec)',
                    default=DEFAULTRTGAP,
                    metavar='W',
                    action='store',
                    type=float)
parser.add_argument('inputFile',
                    help='score data input data file',
                    type=str)
parser.add_argument('outputFile',
                    help='file name to save text data to',
                    type=str)
parser.add_argument('--logFile',
                    metavar='FILENAME',
                    type=str,
                    help='log progress in FILENAME')
parser.add_argument('--mzGap',
                    type=float,
                    default=DEFAULTMZGAP,
                    help='m/z gap size for finding local maxima (m/z)')
parser.add_argument('--removeLow',
                    type=int,
                    help='Remove score values below the given signal level')
parser.add_argument('--outDir',
                    metavar='DIRECTORY',
                    type=str,
                    required=False,
                    help='save output in DIRECTORY, if it does not exist it will be created')


def main():
    options = parser.parse_args()

    if options.logFile:
        logging.basicConfig(filename=options.logFile + '-' + str(RANK),
                            level=logging.DEBUG,
                            filemode='w',
                            format='%(asctime)s %(message)s',
                            datefmt='%m/%d/%Y %H:%M:%S')
    logging.info('program started')
    logging.info('command line: {}'.format(' '.join(sys.argv)))

    data_out = open(options.outputFile, "w")
    Drt = options.rtGap/2.0
    Dmz = options.mzGap/2.0
    idx = index.Index()
    logging.info('region: +/- {} sec, +/- {} mz'.format(Drt,Dmz))
    count = 0
    with open(options.inputFile, "r") as data_in:
        for datum in data_in:
#            rt, mz, amp, score = (float(v) for v in datum.split(","))
            rt, mz, amp, score = [float(v) for v in datum.split(",")][:4]
            coord = (rt-Drt, mz-Dmz, rt+Drt, mz+Dmz)
            if idx.count(coord) == 0:
                idx.insert(count, coord)
                print >> data_out, "{}, {}, {}, {}".format(rt, mz, amp, score)
                count += 1
    print "found {} regions".format(count)

    logging.debug('main mem: {}'.format(md_io.memory_usage_resource()))
    logging.info('program completed')


if __name__ == '__main__':

    main()
