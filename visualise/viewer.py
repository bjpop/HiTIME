#!/bin/env python

import sys
import os.path
import resource
import math
import argparse
import logging
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import SpanSelector

# add dir of this (and following) file to path
sys.path.append(os.path.realpath(__file__))
import md_filter
import md_io


# default difference in mass of isotopes
defaultMzDelta = 6.0201
# default m/z tolerance in parts per million
defaultPPM = 4.0
# fwmh ppm
#defaultFWHM = 12.0
#defaultFWHM = 1000.0
defaultFWHM = 120.0
# default ration of isotopes
defaultIntensityRatio = 1.0
#default full width at half max height of signal for retention time. In seconds.
defaultRTwidth = 9.0

# command line argument parser
parser = argparse.ArgumentParser(description='View mass-spec files and isotope doublets.')
parser.add_argument('--format',
                    help='file format used for input mass spec data, options are: mzml; mzdata; score',
                    default='score',
                    type=str)
parser.add_argument('inputFile',
                    help='mass spec input data file in mzdata (XML) format OR score data from hitime',
                    type=str)


def ReadScores(options):
    # input line of format: id, rt, mz, amp, score
    ids = []
    rts = []
    mzs = []
    amps = []
    scores = []
    with open(options.inputFile, 'r') as fin:
        for line in fin:
            # clean data to be space separated
            vals = line.strip().replace(',', ' ').split()
            try:
                ids.append(int(vals[0]))
                rts.append(float(vals[1]))
                mzs.append(float(vals[2]))
                amps.append(float(vals[3]))
                scores.append(float(vals[4]))
            except (ValueError, IndexError):
                    continue
    if ids:
        return np.array(ids), np.array(rts), np.array(mzs), np.array(amps), np.array(scores)
    return None, None, None, None, None


def DoImage(ids, times, mzs, amps, scores):


    # axis sharing for auto scaling
    fig = plt.figure()
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222, sharex=ax1, sharey=ax1)
    ax3 = fig.add_subplot(223, sharex=ax1, sharey=ax1)
    ax3 = fig.add_subplot(224, sharex=ax1, sharey=ax1)


    plt.figure()
    plt.subplot(4,1,1)
    rt_times = np.unique(times)
    # colapse mz data to give TIC
    TIC, _, _ = np.histogram2d(times, mzs, bins=(len(rt_times), 1), weights=amps)
    plt.plot(rt_times, TIC)
    TIC, _, _ = np.histogram2d(times, mzs, bins=(len(rt_times), 1), weights=scores)
    plt.plot(rt_times, TIC)

    plt.subplot(4,1,2)
    rt_view = 700.0
    rtwidth = 9
    scan_inx = (np.where((times > rt_view - rtwidth) & (times < rt_view + rtwidth)))
#    X, _, _ = np.histogram2d(mzs[scan_inx], times[scan_inx], bins=(len(scan_inx), 1), weights=amps[scan_inx])

#    plt.plot(mzs[scan_inx], X)
    plt.scatter(mzs[scan_inx], amps[scan_inx], marker='x', s=4, c=times[scan_inx], cmap=plt.cm.Reds )
    plt.scatter(mzs[scan_inx], scores[scan_inx], marker='x', s=4, c=times[scan_inx], cmap=plt.cm.Reds )

#    plt.subplot(2,1,2, projection='3d')

#    plt.scatter(times[scan_inx], mzs[scan_inx], c=amps[scan_inx], m='s', cmap=plt.cm.Reds, s=4)


    plt.show()

    ## expts
#    heatmap, xedges, yedges = np.histogram2d(times, mzs, bins=(len(rt_times), 1000), weights=amps)
#    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

#    plt.imshow(heatmap, extent=extent, cmap=plt.cm.Reds, aspect='auto')
#    plt.colorbar()
#    plt.scatter(x=times, y=mzs, c=scores, s=10, cmap=plt.cm.OrRd, alpha=0.9, linewidths=(0,))
#    plt.scatter(x=times, y=mzs, c=amps, s=4, cmap=plt.cm.Reds, marker='s', alpha=0.5, linewidths=(0,))


def main():

    options = parser.parse_args()

    # read the input data file and extract useful contents
    if options.format == 'mzml':
        numSpectra, delta_time, medianSignal = md_io.MZMLstats(options)
        spectra = md_io.MZMLtoSpectrum(options, delta_time)
    elif options.format == 'mzdata':
        numSpectra, delta_time, medianSignal = md_io.MZDATAstats(options)
        # TODO: MZDATA handling
        spectra = md_io.MZDATAtoSpectrum(options, delta_time)
    elif options.format == 'score':
        ids, times, mzs, amps, scores = ReadScores(options)
        numSpectra = len(ids)
    else:
        exit("Unknown mass spec data format: {}".format(options.format))

    if numSpectra == 0:
        exit("Zero spectra read from mz data file, did you specify the wrong input format?")


    DoImage(ids, times, mzs, amps, scores)


if __name__ == '__main__':

    main ()
