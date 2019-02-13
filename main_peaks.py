# !/usr/bin/env python3
# -*- coding: utf-8 -*-

# link fuer einen vergleich von algorithmen?
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2631518/

import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import glob
from pprint import pprint
from scipy.signal import savgol_filter

import test_lmfit as lmfit

if __name__ == "__main__":
  """
  dataset_paths = [
    'data/07_Serie_Prototyp/13-N-0Grad2.SSM',
    'data/07_Serie_Prototyp/13-N-hinten2,5Grad.SSM',
    'data/07_Serie_Prototyp/13-N-links2,5Grad.SSM',
    'data/07_Serie_Prototyp/13-N-rechts2,5Grad.SSM',
    'data/07_Serie_Prototyp/13-N-vorn2,5Grad.SSM'
  ]
  """

  # return list of paths of all files in folder 
  dataset_paths = glob.glob('/home/darlokh/Documents/hiwi/code/data/Messwerte_clean/20181130_12_00*')
  for path in dataset_paths:
    # skiprows set to 0 because data was cleaned
    data = np.loadtxt(path, skiprows=0, unpack=True)
    smoothed_signal = savgol_filter(data, 301, 2)

    x = data[0] # wavelength
    y = data[1] # amplitude
    smooth_y = smoothed_signal[1]
    order_value = 20
    peaks, _ = sp.signal.find_peaks(smooth_y, distance=order_value, height=0.0025)
    # valleys = sp.signal.argrelmin(smooth_y, order=order_value)

    distance = 0
    hills = lmfit.isolate_hills(smoothed_signal, peaks, distance)

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

    fig = plt.figure()

    plt.plot(smoothed_signal[0], smoothed_signal[1], "r-")
    plt.plot(smoothed_signal[0][peaks], smoothed_signal[1][peaks], "gx")
    plt.plot(left_valleys[0], left_valleys[1], "go")
    plt.plot(right_valleys[0], right_valleys[1], "ro")
    plt.plot(x_range[0], polys[0](x_range[0]), "b-")
    plt.plot(x_range[1], polys[1](x_range[1]), "b-")
    plt.plot(x_range[2], polys[2](x_range[2]), "b-")

    toolbar_elements = fig.canvas.manager.toolbar.winfo_children()
    pprint(toolbar_elements)

    plt.show()