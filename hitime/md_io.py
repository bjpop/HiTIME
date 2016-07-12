#!/bin/env python

#from lxml import etree
import sys
import resource
import base64
import struct
import numpy as np
from itertools import *
import math
import csv
import logging
import os
import os.path
#import pymzml
from collections import deque
import resource
import pyopenms
import copy

# add dir of this (and following) file to path
sys.path.append(os.path.realpath(__file__))
import md_filter


# helper funtion for memory profiling
def memory_usage_resource():
    import resource
    rusage_denom = 1024
    if sys.platform == 'darwin':
        # ... it seems that in OSX the output is different units ...
        rusage_denom = rusage_denom * rusage_denom
    mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / float(rusage_denom)
    return mem

endianMap = { 'little': '<', 'big': '>' }
precisionMap = { '64' : 'd', '32': 'f' }


# convert the base64 encoded data into an array of floating point values
# given the endianness and precision of the raw data
def interpretBinary(data, endian, precision):
    precisionChar = precisionMap[precision]
    endianChar = endianMap[endian]
    decoded = base64.b64decode(data)
    count = len(decoded) / struct.calcsize(endianChar + precisionChar)
    return struct.unpack(endianChar + precisionChar * count, decoded[0:len(decoded)])

# find the binary data, endianness and precision of a single spectrum element
def getMZDATASpectrum(spectrum, tag):
    for child in spectrum:
        if child.tag == tag:
            for binaryChild in child:
                if binaryChild.tag == 'data':
                    endian = binaryChild.get('endian')
                    precision = binaryChild.get('precision')
                    binary = binaryChild.text
                    return interpretBinary(binary, endian, precision)

# find the timestamp of a spectrum, returns None if no timestamp was found.
def getSpectrumTime(spectrum):
    desc = spectrum.find('spectrumDesc')
    settings = desc.find('spectrumSettings')
    instrument = settings.find('spectrumInstrument')
    for param in instrument.iter('cvParam'):
        if param.get('name') == 'TimeInMinutes':
            return param.get('value')
    return None


# an encapsulation of a single Spectrum element containing its
# identity, decoded mz data and decoded intensity data
class Spectrum(object):
    def __init__(self, id, time, mzs, intensities, original):
        self.time = time
        self.mzs = mzs
        self.intensities = intensities
        self.id = int(id)
        self.original = original


'''
def parseMZDATA(options):
    filename = options.inputFile
    result = []
    # parse the XML document
    tree = etree.parse(filename)
    # get the root element
    root = tree.getroot()
    # iterate over the spectrum elements
    for spectrum in root.iter('spectrum'):
        # get the mz data for the spectrum
        mzData = getMZDATASpectrum(spectrum, 'mzArrayBinary')
        # get the intensity data for the spectrum
        intData = getMZDATASpectrum(spectrum, 'intenArrayBinary')
        time = getSpectrumTime(spectrum)
        result.append(Spectrum(spectrum.get('id'), time, mzData, intData))
    return result
'''

def writeResults(output_experiment, spectrum, scores):
    new_spectrum = copy.deepcopy(spectrum.original)
    scores = scores[:,0]
    #print(type(scores))
    #print(scores)
    #print(len(scores))

    new_spectrum.set_peaks((spectrum.mzs, scores))
    #old_peaks = new_spectrum.get_peaks()
    #old_mzs = old_peaks[0]
    #old_ints = old_peaks[1]
    #print(type(old_mzs))
    #print(type(old_ints))
    #print(old_mzs)
    #print(old_ints)
    #print(len(old_mzs))
    #exit(1) 
    #new_spectrum.set_peaks(new_spectrum.get_peaks())
    output_experiment.addSpectrum(new_spectrum)

'''
def writeResults(stream, spectrum, scores=None):
    if scores is not None:
        rt = spectrum.time
        for mz, amp, val in zip(spectrum.mzs, spectrum.intensities, scores):
#            if val > 0.0:
#                print >> stream, '{}, {}, {}, {}'.format(rt, mz, amp, val)
            if val[0] > 0.0:
                print >> stream, '{}, {}, {}, {}'.format(rt, mz, amp, ', '.join([str(v) for v in val]))
    else:
        rt = spectrum.time
        for mz, amp in zip(spectrum.mzs, spectrum.intensities):
            print >> stream, '{}, {}, {}'.format(rt, mz, amp)
'''

'''
def MZMLtoSpectrum(options):
    filename = options.inputFile
    delta_time = 0
    time_prev = 0
    points = 0
    mean = 0
    time = 0
    msrun = pymzml.run.Reader(filename)
    for n,spectrum in enumerate(msrun):
        mzData = np.array(spectrum.mz, dtype="float32")
        intData = np.array(spectrum.i, dtype="uint64")
        points += len(intData)
        mean += sum(intData)
        try:
            time = spectrum['MS:1000016']
            delta_time += (time - time_prev - delta_time)/(n+1)  # incremental update to mean delta_time
            time_prev = time
        except KeyError:
            time_prev = time
            if delta_time > 0:
                time += delta_time
            else:
                time += 1.0
        yield Spectrum(n, time, mzData, intData)

    if points > 0:
        mean /= float(points)
    else:
        exit("Zero spectra read from mz data file, did you specify the wrong input format?")
    logging.info('mzdata input file parsed, {0} ({1}) spectra (data points) read in'.format(n+1, points))
    logging.info('time delta: %g, mean signal: %g' % (delta_time, mean))
'''


def MZMLtoSpectrum(options):
    filename = options.inputFile
    delta_time = 0
    time_prev = 0
    points = 0
    mean = 0
    time = 0
    mzml_file = pyopenms.MzMLFile()
    experiment = pyopenms.MSExperiment()
    mzml_file.load(filename, experiment)
    for n,spectrum in enumerate(experiment):
        (mzData, intData) = spectrum.get_peaks()
        points += len(intData)
        mean += sum(intData)
        try:
            time = spectrum.getRT()
            delta_time += (time - time_prev - delta_time)/(n+1)  # incremental update to mean delta_time
            time_prev = time
        except KeyError:
            time_prev = time
            if delta_time > 0:
                time += delta_time
            else:
                time += 1.0
        yield Spectrum(n, time, mzData, intData, spectrum)

    if points > 0:
        mean /= float(points)
    else:
        exit("Zero spectra read from mz data file, did you specify the wrong input format?")
    logging.info('mzdata input file parsed, {0} ({1}) spectra (data points) read in'.format(n+1, points))
    logging.info('time delta: %g, mean signal: %g' % (delta_time, mean))


'''
def MZMLtoSpectrum(options):
    filename = options.inputFile
    delta_time = 0
    time_prev = 0
    points = 0
    mean = 0
    time = 0
    msrun = pymzml.run.Reader(filename)
    for n,spectrum in enumerate(msrun):
        mzData = np.array(spectrum.mz, dtype="float32")
        intData = np.array(spectrum.i, dtype="uint64")
        points += len(intData)
        mean += sum(intData)
        try:
            time = spectrum['MS:1000016']
            delta_time += (time - time_prev - delta_time)/(n+1)  # incremental update to mean delta_time
            time_prev = time
        except KeyError:
            time_prev = time
            if delta_time > 0:
                time += delta_time
            else:
                time += 1.0
        yield Spectrum(n, time, mzData, intData)

    if points > 0:
        mean /= float(points)
    else:
        exit("Zero spectra read from mz data file, did you specify the wrong input format?")
    logging.info('mzdata input file parsed, {0} ({1}) spectra (data points) read in'.format(n+1, points))
    logging.info('time delta: %g, mean signal: %g' % (delta_time, mean))
'''


def nextWindow(reader, options, half_window):
    '''
    Use iterators to serve up data when needed
    '''
    pad_front = repeat(Spectrum(0, 0.0, [0.0], [0.0], None), half_window + 1)  # extra one at start that gets ignored
    pad_back = repeat(Spectrum(0, 0.0, [0.0], [0.0], None), half_window)
    items = chain(pad_front, reader(options), pad_back)
    # 1st window
    data = list(islice(items, 0, 2 * half_window + 1 ))
    data_deque = deque(data)
    # rest
    for i, scan in enumerate(items):
        data_deque.popleft()
        data_deque.append(scan)
        yield list(data_deque)
