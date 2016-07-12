#!/bin/env python

import sys
import os.path
import math
import argparse
import logging
import time
import numpy as np
import pyopenms

# add dir of this (and following) file to path
sys.path.append(os.path.realpath(__file__))
import md_filter
import md_io

# default difference in mass of isotopes
DEFAULTMZDELTA = 6.0201
# default m/z tolerance in parts per million
DEFAULTPPM = 4.0
# fwmh ppm
DEFAULTFWHM = 150.0
#DEFAULTMZSIGMA = 3.0
DEFAULTMZSIGMA = 1.5
# default ration of isotopes
DEFAULTINTENSITYRATIO = 1.0
#default full width at half max height of signal for retention time. In number of scans.
DEFAULTRTWIDTH = 17.0  # scans
#DEFAULTRTSIGMA = 3.0
DEFAULTRTSIGMA = 1.5

# min number of samples in score regions
DEFAULTMINSAMPLE = DEFAULTRTWIDTH * DEFAULTRTSIGMA/2.355

# command line argument parser
parser = argparse.ArgumentParser(description='Filter mass spec data for isotope doublets.')
parser.add_argument('--format',
                    help='file format used for input mass spec data, options are: mzml; mzdata',
                    default='mzml',
                    type=str)
parser.add_argument('--intensityRatio',
                    help='ratio of intensities for a doublet (isotope amount/parent amount)',
                    default=DEFAULTINTENSITYRATIO,
                    metavar='R',
                    action='store',
                    type=float)
parser.add_argument('--rtWidth',
                    help='Retention Time full width at half maximum in number of scans',
                    default=DEFAULTRTWIDTH,
                    metavar='W',
                    action='store',
                    type=float)
parser.add_argument('--rtSigma',
                    help='Boundary for retention time width in standard deviations',
                    default=DEFAULTRTSIGMA,
                    action='store',
                    type=float)
parser.add_argument('--ppm',
                    help='m/z tolerance in parts per million',
                    default=DEFAULTPPM,
                    metavar='P',
                    action='store',
                    type=int)
parser.add_argument('--mzWidth',
                    help='m/z full width at half maximum in parts per million',
                    default=DEFAULTFWHM,
                    metavar='F',
                    action='store',
                    type=int)
parser.add_argument('--mzSigma',
                    help='Boundary for mz window in standard deviations',
                    default=DEFAULTMZSIGMA,
                    action='store',
                    type=float)
parser.add_argument('inputFile',
                    help='mass spec input data file',
                    type=str)
parser.add_argument('outputFile',
                    help='file name to save text data to',
                    type=str)
parser.add_argument('--logFile',
                    metavar='FILENAME',
                    type=str,
                    help='log progress in FILENAME')
parser.add_argument('--mzDelta',
                    metavar='D',
                    type=float,
                    default=DEFAULTMZDELTA,
                    help='m/z difference for doublets')
parser.add_argument('--removeLow',
                    type=int,
                    help='Remove intensity values below the given signal level')
parser.add_argument('--outDir',
                    metavar='DIRECTORY',
                    type=str,
#                    required=True,
                    required=False,
                    help='save output in DIRECTORY, if it does not exist it will be created')
parser.add_argument('--noScore',
                    action='store_true',
                    help='process without scoring.  Use for data exploration')
parser.add_argument('--minSample',
                    help='minimum number of data points required in each sample region',
                    default=DEFAULTMINSAMPLE,
                    action='store',
                    type=float)


def main(MPI=None):

    # try parallel
    COMM = None
    RANK = 0
    SIZE = 1
    if MPI:
        COMM = MPI.COMM_WORLD
        RANK = COMM.Get_rank()
        SIZE = COMM.Get_size()
        status = MPI.Status()

    if RANK == 0:
        options = parser.parse_args()
    else:
        options = []

    if SIZE > 1:
        options = COMM.bcast(options)

    if options.logFile:
        logging.basicConfig(filename=options.logFile + '-' + str(RANK),
                            level=logging.DEBUG,
                            filemode='w',
                            format='%(asctime)s %(message)s',
                            datefmt='%m/%d/%Y %H:%M:%S')
    logging.info('program started')
    logging.info('command line: {}'.format(' '.join(sys.argv)))

    if RANK == 0:
        if options.format == 'mzml':
            reader = md_io.MZMLtoSpectrum
        elif options.format == 'mzdata':
            # TODO: convert to MPI
            reader = md_io.parseMZDATA
        else:
            exit("Unknown mass spec data format: {}".format(options.format))


    # distribute to worker pool
    #           design
    # master            | workers
    # ---------  --------------------------------
    # while workers     | while !end
    #                   |   send result (or None)
    #   recv result     |
    #   store result    |
    #   get sender      |
    #   read data chunk |
    #   send to sender  |
    #                   |   recv data (or None)
    #                   |   do work (or end)
    #
    #
    # main split b/w master worker

    # read the input data file and extract useful contents
    if RANK == 0:
        #data_out = open(options.outputFile, "w")
        output_file = pyopenms.MzMLFile()
        output_experiment = pyopenms.MSExperiment()
        half_window = int(math.ceil(options.rtSigma * options.rtWidth / 2.355))
        logging.debug('half RT window {}'.format(half_window))
        if SIZE > 1:
            done = 1
        else:
            done = 0  # only needed for sequential
        raw_data = None  # only needed for sequential
        if not options.noScore:
            nextWindow = md_io.nextWindow(reader, options, half_window)
        else:
            nextWindow = reader(options)
        while done < SIZE:
            if SIZE > 1:
                raw_data, scores = COMM.recv(source=MPI.ANY_SOURCE, status=status)
                source = status.Get_source()
            if raw_data is not None:
                #md_io.writeResults(data_out, raw_data, scores)
                md_io.writeResults(output_experiment, raw_data, scores)

            ## Read data chunk
            try:
                spectra = nextWindow.next()
            except StopIteration:
                spectra = None

            # try removing low values
            if options.removeLow > 0 and spectra is not None:
                # use the specified low signal
                spectra = md_filter.removeLowSignal(spectra, options.removeLow)

            if SIZE > 1:
                COMM.send(spectra, dest=source)
            elif spectra is not None:    # do work sequentially
                ## do work
                if not options.noScore:
                    scores = md_filter.scoreSpectra(spectra, options)
                    raw_data = spectra[len(spectra)//2]
                else:
                    scores = None
                    raw_data = spectra
            if spectra is None:
                done += 1   # can only ever close each worker once
        output_file.store(options.outputFile, output_experiment)
    else:  # Worker
        scores = None
        raw = None
        # for stats
        send_time = []
        recv_time = []
        work_time = []
        in_mem = md_io.memory_usage_resource()
        while True:
            t1 = time.time()
            COMM.send((raw, scores))
            t2 = time.time()
            send_time.append(t2 - t1)
            raw = None
            scores = None
            logging.debug('rank {}, mem start, end: {:.1f} {:.1f}'.format(RANK, in_mem, md_io.memory_usage_resource()))
            in_mem = md_io.memory_usage_resource()
            t1 = time.time()
            spectra = COMM.recv()
            t2 = time.time()
            recv_time.append(t2 - t1)
            if spectra is not None:
                ## do work
                try:
                    t1 = time.time()
                    if not options.noScore:
                        scores = md_filter.scoreSpectra(spectra, options)
                        raw = spectra[len(spectra)//2]
                    else:
                        scores = None
                        raw = spectra
                    t2 = time.time()
                    work_time.append(t2 - t1)
                    spectra = None
                except MemoryError:
                    logging.debug('rank {} Memory Error'.format(RANK))
            else:
                break
        logging.info('rank {}, count {}'.format(RANK, len(work_time)))
        logging.info('rank {}, stats min 10% 25% 50% 75% 90% max'.format(RANK))
        limits = np.percentile(send_time, [0,10,25,50,75,90,100])
        logging.info('rank {}, send {}'.format(RANK, ', '.join(['{:.2f}'.format(i) for i in limits])))
        limits = np.percentile(recv_time, [0,10,25,50,75,90,100])
        logging.info('rank {}, recv {}'.format(RANK, ', '.join(['{:.2f}'.format(i) for i in limits])))
        limits = np.percentile(work_time, [0,10,25,50,75,90,100])
        logging.info('rank {}, work {}'.format(RANK, ', '.join(['{:.2f}'.format(i) for i in limits])))

    logging.debug('main mem: {}'.format(md_io.memory_usage_resource()))
    logging.info('program completed')


if __name__ == '__main__':

    # try parallel
    try:
        from mpi4py import MPI
    except ImportError:
        MPI = None

    main(MPI)
