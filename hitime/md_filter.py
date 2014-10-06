#!/bin/env python

import resource
import numpy as np
import math
import logging
#import numexpr as ne
#from bx.intervals.intersection import IntervalTree

#from memory_profiler import profile

# For testing
def math_err(val):
    return np.nan_to_num(val)
#    return val

ROOT2PI = np.sqrt(2.0*np.pi)
#@profile
def scoreSpectra(spectra, options):
    '''
        Calculate correlation for each point using all points in region of hi/lo and rt window
        centred at point and point + delta mz.
        Calculate total correlation based on estimates in each region rescaled by proportion of data
        Use gaussian as ideal shape and rescale hi region by ratio of isotopes
        Gaussian extent is bounds of region.
        Define
        SSX = (x - mean(x))^2
        SSY = (f - mean(f))^2
        SXY = (x - mean(x))(f - mean(f))

        Correlation and linear regression from:
        http://en.wikipedia.org/wiki/Correlation_and_dependence
        http://en.wikipedia.org/wiki/Simple_linear_regression

        Tried parallel stats, but too slow AND memory intensive
        Parallel stats from:
        http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
    '''
    if spectra is None:
        return None

    rt_sigma = options.rtWidth / 2.355     # rt sigma from fwhm num scans
    mz_ppm_sigma = options.mzWidth/2.355e6    # half width ppm to sigma ppm
    mz_delta = options.mzDelta
    iso_ratio = options.intensityRatio
    minSample = options.minSample

    rt_len = len(spectra)
    mid_win = rt_len//2
    # mz values in centre RT of window
    # all regions are centred on these points
    mz_mu_vect = spectra[mid_win].mzs

    # create list of all start and stop thresholds
    # use 2 * sigma bounds
    lo_tol = 1.0 - options.mzSigma * mz_ppm_sigma
    hi_tol = 1.0 + options.mzSigma * mz_ppm_sigma

    points_lo_vect = np.vstack((mz_mu_vect * lo_tol, mz_mu_vect * hi_tol)).ravel(order='F')
    points_hi_vect = np.vstack(((mz_mu_vect + mz_delta) * lo_tol, (mz_mu_vect + mz_delta) * hi_tol)).ravel(order='F')

    data_lo = [[] for _ in range(mz_mu_vect.shape[0])]
    data_hi = [[] for _ in range(mz_mu_vect.shape[0])]
    shape_lo = [[] for _ in range(mz_mu_vect.shape[0])]
    shape_hi = [[] for _ in range(mz_mu_vect.shape[0])]
    len_lo = np.zeros(mz_mu_vect.shape)
    len_hi = np.zeros(mz_mu_vect.shape)

    # shape in rt dir
    pts = np.true_divide(np.arange(rt_len)-mid_win, float(rt_sigma))
    rt_shape = np.true_divide(np.exp( -0.5*pts*pts ), (rt_sigma * ROOT2PI))

    # for each scan row in rt window
    for rowi in range(rt_len):

#        logging.debug('window line: %d' % rowi)

        # mzs and amp per scan
        row_amp = spectra[rowi].intensities
        row_mzs = spectra[rowi].mzs

        # determine window bounds
        bounds_lo = np.searchsorted(row_mzs, points_lo_vect)  # find start and stop indices. assumes row_mzs sorted
        bounds_lo.shape = bounds_lo.shape[0]/2, 2

        bounds_hi = np.searchsorted(row_mzs, points_hi_vect)
        bounds_hi.shape = bounds_hi.shape[0]/2, 2

        # Numpy doesn't really vectorise array of irregular sized arrays
        # so must loop over all mz (masked array too memory intensive)
        # for each mz in scan row

        # row weight
        rt_lo = rt_shape[rowi]
#        rt_hi = iso_ratio * rt_lo
        rt_hi = rt_lo  # weight by iso_ratio later
        for mzi in range(mz_mu_vect.shape[0]):
            l1,r1 = bounds_lo[mzi]
            l2,r2 = bounds_hi[mzi]

            centre = mz_mu_vect[mzi]
            sigma = centre * mz_ppm_sigma

            if r1 > l1:
                mzs = (np.array(row_mzs[l1:r1]) - centre)/float(sigma)
                amp_lo = np.array(row_amp[l1:r1])
                shape_add = rt_lo * np.exp( -0.5 * mzs * mzs) / ( sigma * ROOT2PI )
                len_lo[mzi] += r1 - l1
            else:
                # all other terms at default (0 or 1 for count)
                amp_lo = np.zeros(1)
                shape_add = np.array([rt_lo/(sigma * ROOT2PI)])

            data_lo[mzi].extend(amp_lo)
            shape_lo[mzi].extend(shape_add)

            centre += mz_delta
            sigma = centre * mz_ppm_sigma

            if r2 > l2:
                mzs = (np.array(row_mzs[l2:r2]) - centre)/float(sigma)
                amp_hi = np.array(row_amp[l2:r2])
                shape_add = rt_hi * np.exp( -0.5 * mzs * mzs) / ( sigma * ROOT2PI )
                len_hi[mzi] += r2 - l2
            else:
                # all other terms at default (0 or 1 for count)
                amp_hi = np.zeros(1)
                shape_add = np.array([rt_hi/(sigma * ROOT2PI)])

            data_hi[mzi].extend(amp_hi)
            shape_hi[mzi].extend(shape_add)

    # convert to numpy arrays
    data_lo = [ np.array(v) if l >= minSample else np.zeros(1) for v,l in zip(data_lo, len_lo) ]
    data_hi = [ np.array(v) if l >= minSample else np.zeros(1) for v,l in zip(data_hi, len_hi) ]
    shape_lo = [ np.array(v) if l >= minSample else np.zeros(1) for v,l in zip(shape_lo, len_lo) ]
    shape_hi = [ iso_ratio * np.array(v) if l >= minSample else np.zeros(1) for v,l in zip(shape_hi, len_hi) ]
#    data_lo = [ np.array(v) for v in data_lo ]
#    data_hi = [ np.array(v) for v in data_hi ]
#    shape_lo = [ np.array(v) for v in shape_lo ]
#    shape_hi = [ iso_ratio * np.array(v) for v in shape_hi ]

    # data and shape for component regions
    dataAB = [ np.hstack((len(hi)*lo, len(lo)*hi)) for lo, hi in zip(data_lo, data_hi) ]
    shapeAB = [ np.hstack((len(hi)*lo, len(lo)*hi)) for lo, hi in zip(shape_lo, shape_hi) ]
    shapeA0 = [ np.hstack((len(hi)*lo, np.zeros(len(hi)))) for lo, hi in zip(shape_lo, shape_hi) ]
    shapeB0 = [ np.hstack((np.zeros(len(lo)), len(lo)*hi)) for lo, hi in zip(shape_lo, shape_hi) ]
    shape1r = [ np.hstack((len(hi)*np.ones(len(lo)), iso_ratio*len(lo)*np.ones(len(hi)))) for lo, hi in zip(shape_lo, shape_hi) ]
    nAB = np.array([ len(v) for v in dataAB ])

    # attempt to free memory
    # let python know these arn't needed any more
    del data_lo
    del data_hi
    del shape_lo
    del shape_hi

    # centre
    dataAB = [ v - np.mean(v) for v in dataAB ]
    shapeAB = [ v - np.mean(v) for v in shapeAB ]
    shapeA0 = [ v - np.mean(v) for v in shapeA0 ]
    shapeB0 = [ v - np.mean(v) for v in shapeB0 ]
    shape1r = [ v - np.mean(v) for v in shape1r ]

    # square
    data2AB = [ v*v for v in dataAB ]
    shape2AB = [ v*v for v in shapeAB ]
    shape2A0 = [ v*v for v in shapeA0 ]
    shape2B0 = [ v*v for v in shapeB0 ]
    shape21r = [ v*v for v in shape1r ]

    # sum of squares
    SSY = np.array([np.sum(v) for v in data2AB])
    SSXAB = np.array([np.sum(v) for v in shape2AB])
    SSXA0 = np.array([np.sum(v) for v in shape2A0])
    SSXB0 = np.array([np.sum(v) for v in shape2B0])
    SSX1r = np.array([np.sum(v) for v in shape21r])

    # mem hint
    del data2AB
    del shape2AB
    del shape2A0
    del shape2B0
    del shape21r

    # data * shape for each type of correlation
    datashape = [ v*u for v,u in zip(dataAB, shapeAB) ]
    SXYAB = np.array([np.sum(v) for v in datashape])
    datashape = [ v*u for v,u in zip(dataAB, shapeA0) ]
    SXYA0 = np.array([np.sum(v) for v in datashape])
    datashape = [ v*u for v,u in zip(dataAB, shapeB0) ]
    SXYB0 = np.array([np.sum(v) for v in datashape])
    datashape = [ v*u for v,u in zip(dataAB, shape1r) ]
    SXY1r = np.array([np.sum(v) for v in datashape])
    datashape = [ v*u for v,u in zip(shapeAB, shapeA0) ]
    SXYABA0 = np.array([np.sum(v) for v in datashape])
    datashape = [ v*u for v,u in zip(shapeAB, shapeB0) ]
    SXYABB0 = np.array([np.sum(v) for v in datashape])
    datashape = [ v*u for v,u in zip(shapeAB, shape1r) ]
    SXYAB1r = np.array([np.sum(v) for v in datashape])

    del dataAB
    del shapeAB
    del shapeA0
    del shapeB0
    del shape1r

    # full correlation
    correlAB = math_err(np.true_divide(SXYAB, np.sqrt(SSXAB*SSY))).clip(min=0)

    # Region lo
    correlA0 = math_err(np.true_divide(SXYA0, np.sqrt(SSXA0*SSY))).clip(min=0)

    # Region hi
    correlB0 = math_err(np.true_divide(SXYB0, np.sqrt(SSXB0*SSY))).clip(min=0)

    # no Region
    correl1r = math_err(np.true_divide(SXY1r, np.sqrt(SSX1r*SSY))).clip(min=0)

    # Correl between shapes AB, A0
    correlABA0 = math_err(np.true_divide(SXYABA0, np.sqrt(SSXAB*SSXA0))).clip(min=0)

    # Correl between shapes AB, B0
    correlABB0 = math_err(np.true_divide(SXYABB0, np.sqrt(SSXAB*SSXB0))).clip(min=0)

    # Correl between shapes AB, 1r
    correlAB1r = math_err(np.true_divide(SXYAB1r, np.sqrt(SSXAB*SSX1r))).clip(min=0)

    # Use Steiger (a.k.a Meng) Z-test
    rm2ABA0 = 0.5*(correlAB*correlAB + correlA0*correlA0)
    rm2ABB0 = 0.5*(correlAB*correlAB + correlB0*correlB0)
    rm2AB1r = 0.5*(correlAB*correlAB + correl1r*correl1r)
    fABA0 = math_err(np.true_divide((1.0 - correlABA0), 2.0*(1.0-rm2ABA0))).clip(max=1.0)
    fABB0 = math_err(np.true_divide((1.0 - correlABB0), 2.0*(1.0-rm2ABB0))).clip(max=1.0)
    fAB1r = math_err(np.true_divide((1.0 - correlAB1r), 2.0*(1.0-rm2AB1r))).clip(max=1.0)

    hABA0 = math_err(np.true_divide((1.0 - fABA0*rm2ABA0), (1.0 - rm2ABA0)))
    hABB0 = math_err(np.true_divide((1.0 - fABB0*rm2ABB0), (1.0 - rm2ABB0)))
    hAB1r = math_err(np.true_divide((1.0 - fAB1r*rm2AB1r), (1.0 - rm2AB1r)))

    zAB = np.arctanh(correlAB)
    zA0 = np.arctanh(correlA0)
    zB0 = np.arctanh(correlB0)
    z1r = np.arctanh(correl1r)

    sqrtn = np.sqrt(nAB-3.0)

    # Wrong
#    zABA0 = math_err(np.true_divide(np.fabs(zAB-zA0)*sqrtn, np.sqrt(2.0*(1.0-correlABA0)*hABA0)))
#    zABB0 = math_err(np.true_divide(np.fabs(zAB-zB0)*sqrtn, np.sqrt(2.0*(1.0-correlABB0)*hABB0)))
    # Right??
    zABA0 = math_err(np.true_divide((zAB-zA0)*sqrtn, np.sqrt(2.0*(1.0-correlABA0)*hABA0)))
    zABB0 = math_err(np.true_divide((zAB-zB0)*sqrtn, np.sqrt(2.0*(1.0-correlABB0)*hABB0)))
    zAB1r = math_err(np.true_divide((zAB-z1r)*sqrtn, np.sqrt(2.0*(1.0-correlAB1r)*hAB1r)))

#    score = np.minimum(zABA0, zABB0)

# this looks worse
#    Zscore = np.minimum(zABA0, zABB0)
    minScore = np.minimum(zABA0, zABB0, zAB1r).clip(min=0.0)
    score = np.vstack((minScore, correlAB, correlA0, correlB0, correl1r)).T

    return score


def removeLowSignal(spectra, cutoff):
    '''Remove data with intensity below cutoff.'''
    for spectrum in spectra:
        amps = np.array(spectrum.intensities)
        mzs = np.array(spectrum.mzs)
        indices = np.where(amps >= cutoff)  # NOTE! >= cutoff
        spectrum.intensities = amps[indices]
        spectrum.mzs = mzs[indices]
    return spectra


# parts per million tolerence
# individual mz, or average mz
def mzTol(mz, mz_accuracy):
    return mz*mz_accuracy/10**6


def median(xs):
    noZeros = filter(lambda x: x != 0, xs)
    size = len(noZeros)
    if size > 0:
        return sorted(noZeros)[size/2]
    else:
        return None
