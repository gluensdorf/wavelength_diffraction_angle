# !/usr/bin/env python3
# -*- coding: utf-8 -*-

# link fuer einen vergleich von algorithmen?
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2631518/

import math as math
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import glob
import peakutils as pu
from pprint import pprint
from scipy.signal import savgol_filter
from scipy.signal import find_peaks, find_peaks_cwt

dataset_paths = [
  'data/07_Serie_Prototyp/13-N-0Grad2.SSM',
  'data/07_Serie_Prototyp/13-N-hinten2,5Grad.SSM',
  'data/07_Serie_Prototyp/13-N-links2,5Grad.SSM',
  'data/07_Serie_Prototyp/13-N-rechts2,5Grad.SSM',
  'data/07_Serie_Prototyp/13-N-vorn2,5Grad.SSM'
]

dataset_paths = glob.glob('/home/darlokh/Documents/hiwi/code/data/Messwerte_clean/20181130_12_00*')

pass
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

  :param smoothed_signal: 2D-list containing a list of x-values and a list of y-values 
  :param n: amount of peaks to find
  :param peaks: list of x-values of peaks found in the smoothed_signal

  :returns: 2D-list containing a list of x-values of n-highest peaks
                           and a list of y-values of n-highest peaks
  """
  # sort peaks by their amplitude
  sort_peaks = np.argsort(smoothed_signal[1][peaks])

  # top_three_peaks -- from smoothed signal take the peaks and sort them by 'y' and take the last 3 (the biggest three)
  n_highest_peaks = []
  for index in reversed(range(0, n)):
    wavelength = smoothed_signal[0][peaks][sort_peaks][-index]
    amplitude = smoothed_signal[1][peaks][sort_peaks][-index]
    peak_index = list(smoothed_signal[0]).index(wavelength)
    n_highest_peaks.append({
        'peak_index' : peak_index,
        'wavelength' : wavelength,
        'amplitude' : amplitude
      })

  return n_highest_peaks

def isolate_hills(spectrum, peaks, distance=0):
  """
  Returns a list of valley/minima pairs which enclose a peak/maxima in peaks. 
  The pairs are stored in a dictionary {left: valley/minimum, right: valley/minimum}.

  # :param peaks: list of peaks in a smoothed spectrum
  :param spectrum: a smoothed spectrum
  :param distance: default=0, optional, minimal distance required between a peak and left/right valley/minimum
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
  

for path in dataset_paths:
  # skiprows set to 0 because data was cleaned
  data = np.loadtxt(path, skiprows=0, unpack=True)
  smoothed_signal = savgol_filter(data, 301, 2) #<<< useless, TODO delete - is used to plot the smoothed signal

  x = data[0]
  y = data[1]
  smooth_y = smoothed_signal[1]
  order_value = 20
  peaks, _ = sp.signal.find_peaks(smooth_y, distance=order_value, height=0.0025)#order=order_value)
  # valleys = sp.signal.argrelmin(smooth_y, order=order_value)

  distance = 0
  hills = isolate_hills(smoothed_signal, peaks, distance)

  left_valleys = [[], []]
  right_valleys = [[], []]
  for element in hills:
    left_valleys[0].append(smoothed_signal[0][element['left']])
    left_valleys[1].append(smoothed_signal[1][element['left']])
    right_valleys[0].append(smoothed_signal[0][element['right']])
    right_valleys[1].append(smoothed_signal[1][element['right']])
  
  hills_samplepoints = []
  for element in hills:
    hills_samplepoints.append([
      smoothed_signal[0][ element['left']: element['right']],
      smoothed_signal[1][ element['left']: element['right']]
      ])
  
  polys = []
  x_range = []
  for samplepoints in hills_samplepoints:
    polys.append(np.poly1d(np.polyfit(samplepoints[0], samplepoints[1], 2)))
    x_range.append(
      np.linspace(
        samplepoints[0][0], 
        samplepoints[0][-1], 
        500))

  for poly in polys:
    p_der = np.polyder(poly)
    print(p_der.roots)

  plt.plot(smoothed_signal[0], smoothed_signal[1], "r-")
  plt.plot(smoothed_signal[0][peaks], smoothed_signal[1][peaks], "gx")
  plt.plot(left_valleys[0], left_valleys[1], "go")
  plt.plot(right_valleys[0], right_valleys[1], "ro")
  plt.plot(x_range[0], polys[0](x_range[0]), "b-")
  plt.plot(x_range[1], polys[1](x_range[1]), "b-")
  plt.plot(x_range[2], polys[2](x_range[2]), "b-")

  plt.show()