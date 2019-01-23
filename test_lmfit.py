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

def get_valleys_around_peak(peak):
  """
  Returns a list of valley/minima pairs which enclose a peak/maxima in peaks. 
  The pairs are stored in a dictionary {lefthandside: valley/minimum, righthandside: valley/minimum}.

  :param peaks: list of peaks in a smoothed spectrum
  :param spectrum: a smoothed spectrum
  :returns: list of valley/minima pairs for each peak in peaks
  """

  return
  

for path in dataset_paths:
  # skiprows set to 0 because data was cleaned
  data = np.loadtxt(path, skiprows=0, unpack=True)
  smoothed_signal = savgol_filter(data, 301, 2) #<<< useless, TODO delete - is used to plot the smoothed signal

  x = data[0]
  y = data[1]
  smooth_x = smoothed_signal[0]
  smooth_y = smoothed_signal[1]

  peaks, properties = find_peaks(smooth_y, distance=150)#, height=0.02)
  # peaks = find_peaks_cwt(smooth_y, np.arange(50, 75))
  print(peaks)
  print(smooth_y[peaks])
  # prominences, left_bases, right_bases = sp.signal.peak_prominences(smooth_y, peaks, wlen=75)
  # peak_width_full = sp.signal.peak_widths(smooth_y, peaks, rel_height=1)
  # print("prominences: ", prominences)
  # print("left_bases: ", left_bases)
  # print("right_bases: ", right_bases)

  plt.plot(x-1, y, 'b-')
  plt.plot(smoothed_signal[0], smoothed_signal[1], "r-")
  plt.plot(smoothed_signal[0][peaks], smoothed_signal[1][peaks], "yx")

  valleys, properties = find_peaks(-smooth_y, distance=150)#, height=0.02)
  plt.plot(x-1, y, 'b-')
  plt.plot(smoothed_signal[0], smoothed_signal[1], "y-")
  plt.plot(smoothed_signal[0][valleys], smoothed_signal[1][valleys], "rx")
  print(valleys)
  print(smooth_y[valleys])
  # plt.plot(smoothed_signal[0][valleys], smoothed_signal[1][valleys], "yo")
  # plt.plot(smoothed_signal[0][right_bases], smoothed_signal[1][right_bases], "go")
  # plt.plot(smoothed_signal[0][left_bases], smoothed_signal[1][left_bases], "ro")
  # plt.hlines(*peak_width_full[1:], color="lime")

  plt.show()
  """
  widths = np.arange(1, 45)
  smoothed_signal = savgol_filter(data, 45, 2) #<<< useless, TODO delete - is used to plot the smoothed signal
  top_three_peaks = get_top_peaks_cwt(smoothed_signal, widths)
  # print("top_peaks: ", top_peaks)
  """

  """
  # former kernel_width was 45
  top_three_peaks = get_top_peaks(data, 301, 2)
  smoothed_signal = savgol_filter(data, 301, 2) #<<< useless, TODO delete - is used to plot the smoothed signal
  
  # top_three_peaks_int64 = np.array(top_three_peaks, dtype='int64')
  # smoothed_signal.astype(int64)
  # print("smoothed_signal[0]: ", smoothed_signal[0])
  # print("top_three_peaks_int64[0]: ", top_three_peaks_int64[0])
  # print(smoothed_signal[0][592])
  # prominences = sp.signal.peak_prominences(smoothed_signal[0], top_three_peaks_int64[0], wlen=25)
  # print(prominences)

  # print(top_three_peaks[0])

  plt.plot(x-1, y, 'b-')
  plt.plot(smoothed_signal[0], smoothed_signal[1], "r-")
  plt.plot(top_three_peaks[0], top_three_peaks[1], "yx")
  plt.show()
  """
