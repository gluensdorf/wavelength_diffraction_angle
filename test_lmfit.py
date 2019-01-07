# !/usr/bin/env python3
# -*- coding: utf-8 -*-

import math as math
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy.signal import savgol_filter
from scipy.signal import find_peaks

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

    # top_three_peaks -- from smoothed signal take the peaks and sort them by 'y' and take the last 3 (the biggest three)
    x_top_three_peaks = smoothed_signal[0][peaks][sort_peaks][-3:]
    y_top_three_peaks = smoothed_signal[1][peaks][sort_peaks][-3:]
    return (x_top_three_peaks, y_top_three_peaks)

dataset_paths = [
  'data/07_Serie_Prototyp/13-N-0Grad2.SSM',
  'data/07_Serie_Prototyp/13-N-hinten2,5Grad.SSM',
  'data/07_Serie_Prototyp/13-N-links2,5Grad.SSM',
  'data/07_Serie_Prototyp/13-N-rechts2,5Grad.SSM',
  'data/07_Serie_Prototyp/13-N-vorn2,5Grad.SSM'
]

for path in dataset_paths:
  data = np.loadtxt(path, skiprows=2, unpack=True)
  x = data[0]
  y = data[1]

  top_three_peaks = get_top_peaks(data, 45, 2)
  smoothed_signal = savgol_filter(data, 45, 2) #<<< useless, TODO delete 
  print(top_three_peaks[0])

  plt.plot(x-1, y, 'bo')
  plt.plot(smoothed_signal[0], smoothed_signal[1], 'r-')
  plt.plot(top_three_peaks[0], top_three_peaks[1], "yx")
  plt.show()
