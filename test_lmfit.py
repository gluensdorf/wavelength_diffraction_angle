# !/usr/bin/env python3
# -*- coding: utf-8 -*-

# link fuer einen vergleich von algorithmen?
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2631518/

import numpy as np
import scipy as sp
from scipy.signal import savgol_filter, find_peaks, find_peaks_cwt

"""
about savgol_filter:
take the signal data and apply savitzky-golay filter on it
second parameter is yet arbitrarily chosen - has to be odd (to have a center), smooths values which are left and right of center
  concerning padding or what happens at the edges of the signal:
      "When the 'interp' mode is selected (the default), no extension is used. Instead, a degree /polyorder/ polynomial
       is fit to the last /window_length/ values of the edges, and this polynomial is used to evaluate the last /window_length // 2/
       output values."
third parameter is the order of the polynomial used to fit the samples - 2 should be sufficient

is it sufficient to come up with a kernel-width (second parameter) by testing/experiments?
could be sufficient if the length of data (number of measured wavelengths/"pixel") from a
measurement won't change.

tested different kernel-widths, first correct maxima found with 43 using 13-N-0Grad2.SSM as input data
"""


def get_top_peaks(spectrum, kernel_width, polynomial_order):
    """
    Returns list of found top peaks from the given spectrum. Uses savgol_filter (Savitzky-Golay filter) to smooth the spectrum.

    :param spectrum: list of measured spectrum with xy-pairs (wavelength, amplitude)
    :returns: list of top peaks found
    """

    smoothed_signal = savgol_filter(spectrum, kernel_width, polynomial_order)
    peaks, _ = find_peaks(smoothed_signal[1])
    sort_peaks = np.argsort(smoothed_signal[1][peaks])

    print("peaks: ", peaks)
    # top_three_peaks -- from smoothed signal take the peaks and sort them by 'y' and take the last 3 (the biggest three)
    x_top_three_peaks = smoothed_signal[0][peaks][sort_peaks][-3:]
    y_top_three_peaks = smoothed_signal[1][peaks][sort_peaks][-3:]
    return np.array([x_top_three_peaks, y_top_three_peaks])


def get_top_peaks_cwt(spectrum, widths):
    """
    Returns list of found top peaks from the given spectrum. Uses savgol_filter (Savitzky-Golay filter) to smooth the spectrum.

    :param spectrum: list of measured spectrum with xy-pairs (wavelength, amplitude)
    :param widths: sequence of widths required by find_peaks_cwt()
    :returns: list of top peaks found
    """

    # smoothed_signal = savgol_filter(spectrum, kernel_width, polynomial_order)
    # peaks, _ = find_peaks(smoothed_signal[1])
    peak_indices = find_peaks_cwt(spectrum[1], widths)
    sort_peaks = np.argsort(spectrum[1][peak_indices])

    # top_three_peaks -- from smoothed signal take the peaks and sort them by 'y' and take the last 3 (the biggest three)
    x_top_three_peaks = spectrum[0][peak_indices][sort_peaks][-3:]
    y_top_three_peaks = spectrum[1][peak_indices][sort_peaks][-3:]
    return np.array([x_top_three_peaks, y_top_three_peaks])


"""
next step:
extract hill around the found peaks
fit a polynomial function through the hill
extract coefficients of polynomial function
get first order deviation of polynomial function (?)
"""


def get_n_highest_peaks(smoothed_signal, n, peaks):
    """
    Returns the n-highest peaks of the smoothed_signal.

    :param smoothed_signal: 2D numpy array containing a list of x-values and a list of y-values
    :param n: amount of peaks to find
    :param peaks: list of x-values of peaks found in the smoothed_signal

    :returns: 2D-list containing a list of x-values of n-highest peaks
                             and a list of y-values of n-highest peaks
    """
    # sort peaks by their amplitude
    sort_peaks = np.argsort(smoothed_signal[1][peaks])

    # top_three_peaks -- from smoothed signal take the peaks and sort them by 'y' and take the last 3 (the biggest three)
    n_highest_peaks = []

    if n > len(peaks):
        n = len(peaks)

    for i in reversed(range(0, n)):
        wavelength = np.flip(smoothed_signal[0][peaks][sort_peaks], 0)[i]
        amplitude = np.flip(smoothed_signal[1][peaks][sort_peaks], 0)[i]
        peak_index = list(smoothed_signal[0]).index(wavelength)
        n_highest_peaks.append({
            'peak_index': peak_index,
            'wavelength': wavelength,
            'amplitude': amplitude
        })
    return n_highest_peaks


def isolate_hills(spectrum, peaks, distance=0):
    """
    Returns a list of valley/minima pairs which enclose a peak/maxima in peaks.
    Pairs are stored in a dictionary {left: valley/minimum, right: valley/minimum}.

    :param peaks: list of peaks in a smoothed spectrum
    :param spectrum: 2D numpy array of a smoothed spectrum
    # :param distance: default=0, optional, minimal distance required between a peak and left/right valley/minimum
    :returns: dictionary of valley/minima pairs for each peak in peaks
    """

    # smooth_x = spectrum[0]
    smooth_y = spectrum[1]
    order_value = 30
    # peaks = sp.signal.argrelmax(smooth_y, order=order_value)
    valleys = list(sp.signal.argrelmin(smooth_y, order=order_value)[0])
    n_highest_peaks = get_n_highest_peaks(spectrum, 3, peaks)

    hills = []
    for element in n_highest_peaks:
        valleys.append(element['peak_index'])
        valleys.sort()
        valley_lhs = valleys[valleys.index(element['peak_index']) - 1]
        valley_rhs = valleys[valleys.index(element['peak_index']) + 1]
        hills.append({
            'peak_index': element['peak_index'],
            'left': valley_lhs,
            'right': valley_rhs
        })
    return hills


def select_valleys(smoothed_signal, hills):
    """
    Returns a tuple with a list of lhs- and rhs_valleys (with x,y) and a list
    of xy-pairs of all elements between a lhs- and rhs_valley pair.

    :param smoothed_signal: 2D numpy array of a smoothed signal
    :param hills: dictionary of valley/minima surrounding a peaks
    :returns: tuple of lhs_valleys, rhs_valleys and hills_samplepoints
    """

    lhs_valleys = [[], []]
    rhs_valleys = [[], []]
    hills_samplepoints = []
    for element in hills:
        lhs_valleys[0].append(smoothed_signal[0][element['left']])
        lhs_valleys[1].append(smoothed_signal[1][element['left']])
        rhs_valleys[0].append(smoothed_signal[0][element['right']])
        rhs_valleys[1].append(smoothed_signal[1][element['right']])
        hills_samplepoints.append([
            smoothed_signal[0][element['left']: element['right']],
            smoothed_signal[1][element['left']: element['right']]
        ])
    return lhs_valleys, rhs_valleys, hills_samplepoints
