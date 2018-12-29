# !/usr/bin/env python3
# -*- coding: utf-8 -*-

import math as math
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy.signal import savgol_filter
from scipy.signal import find_peaks

"""
x = np.array([
    696.15911, 689.68461, 683.19140, 676.67985, 670.15033,
    663.60325, 657.03900, 650.45797, 643.86058, 637.24725,
    630.61840, 623.97447, 617.31589, 610.64312, 603.95661,
    597.25682, 590.54423, 583.81931, 577.08256, 570.33446,
    563.57552 #, 689.68461, 683.18245, 676.66290, 670.12634,
    # 663.57317, 657.00380, 650.41865, 643.81813, 637.20267,
    # 630.57270, 623.92867, 617.27102, 610.60021, 603.91669,
    # 597.22094, 590.51342, 583.79463, 577.06505, 570.32518,
    # 563.57552, 556.81658 
])
y = np.array([
    696.40479, 696.40479, 696.40479, 696.40479, 696.40479,
    696.40479, 696.40479, 696.40479, 696.40479, 696.40479,
    696.40479, 696.40479, 696.40479, 696.40479, 696.40479,
    696.40479, 696.40479, 696.40479, 696.40479, 696.40479,
    696.40479 #, 689.88437, 689.88437, 689.88437, 689.88437,
    # 689.88437, 689.88437, 689.88437, 689.88437, 689.88437,
    # 689.88437, 689.88437, 689.88437, 689.88437, 689.88437,
    # 689.88437, 689.88437, 689.88437, 689.88437, 689.88437,
    # 689.88437, 689.88437
])
"""

data = np.loadtxt('data/07_Serie_Prototyp/13-N-0Grad2.SSM', skiprows=2, unpack=True)
# data = np.loadtxt('data/07_Serie_Prototyp/13-N-hinten2,5Grad.SSM', skiprows=2, unpack=True)
# data = np.loadtxt('data/07_Serie_Prototyp/13-N-links2,5Grad.SSM', skiprows=2, unpack=True)
# data = np.loadtxt('data/07_Serie_Prototyp/13-N-rechts2,5Grad.SSM', skiprows=2, unpack=True)
# data = np.loadtxt('data/07_Serie_Prototyp/13-N-vorn2,5Grad.SSM', skiprows=2, unpack=True)
x = data[0]
y = data[1]

# about savgol_filter:
# take the signal data and apply savitzky-golay filter on it
# second parameter is yet arbitrarily chosen - has to be odd (to have a center), smooths values which are left and right of center 
#   concerning padding or what happens at the edges of the signal:
#       "When the 'interp' mode is selected (the default), no extension is used. Instead, a degree /polyorder/ polynomial
#        is fit to the last /window_length/ values of the edges, and this polynomial is used to evaluate the last /window_length // 2/
#        output values."
# third parameter is the order of the polynomial used to fit the samples - 2 should be sufficient
result = savgol_filter(data, 51, 2)
peaks, _ = find_peaks(result[1])
sort_peaks = np.argsort(result[1][peaks])

# top_three_peaks -- from smoothed signal take the peaks and sort them by 'y' and take the last 3 (the biggest three)
x_top_three_peaks = result[0][peaks][sort_peaks][-3:]
y_top_three_peaks = result[1][peaks][sort_peaks][-3:]

print(x_top_three_peaks)

plt.plot(x-1, y, 'bo')
plt.plot(result[0], result[1], 'r-')
plt.plot(x_top_three_peaks, y_top_three_peaks, "yx")
plt.show()
