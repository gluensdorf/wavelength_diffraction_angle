# !/usr/bin/env python3
# -*- coding: utf-8 -*-

# TODO: der link sollte hier spaeter raus
# link fuer einen vergleich von algorithmen?
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2631518/

import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import glob
from pprint import pprint
from scipy.signal import savgol_filter

import plotter as plotter
import test_lmfit as lmfit

if __name__ == "__main__":
    # return list of paths of all files in folder
    paths_dataset = glob.glob(
        '/home/darlokh/Documents/hiwi/code/data/Messwerte_clean/20181130_12_00*')

    graphs = []
    for path in paths_dataset:
        # skiprows set to 0 because data was cleaned
        data = np.loadtxt(path, skiprows=0, unpack=True)
        smoothed_signal = savgol_filter(data, 301, 2)

        x = data[0]  # wavelength
        y = data[1]  # amplitude
        smooth_y = smoothed_signal[1]
        order_value = 20
        peaks, _ = sp.signal.find_peaks(
            smooth_y, distance=order_value, height=0.0025)
        # valleys = sp.signal.argrelmin(smooth_y, order=order_value)

        distance = 0
        hills = lmfit.isolate_hills(smoothed_signal, peaks, distance)

        # outsource loop to a function
        lhs_valleys = [[], []]
        rhs_valleys = [[], []]
        for element in hills:
            lhs_valleys[0].append(smoothed_signal[0][element['left']])
            lhs_valleys[1].append(smoothed_signal[1][element['left']])
            rhs_valleys[0].append(smoothed_signal[0][element['right']])
            rhs_valleys[1].append(smoothed_signal[1][element['right']])

        # outsource loop to a function
        # TODO: ist in hills_samplepoints dasselbe wie in lhs/rhs_valleys?
        #       wenn ja koennen die valleys weg
        hills_samplepoints = []
        for element in hills:
            hills_samplepoints.append([
                smoothed_signal[0][element['left']: element['right']],
                smoothed_signal[1][element['left']: element['right']]
            ])

        # outsource loop to a function
        polys = []
        x_range = []
        for samplepoints in hills_samplepoints:
            polys.append(
                np.poly1d(
                    np.polyfit(samplepoints[0], samplepoints[1], 2)))
            x_range.append(
                np.linspace(
                    samplepoints[0][0], samplepoints[0][-1], 500))

        # outsource loop to a function
        # per poly calculate derivative and print zero-crossings
        for poly in polys:
            p_der = np.polyder(poly)
            # print(p_der.roots)

        graphs.append([smoothed_signal, x_range, peaks,
                       lhs_valleys, rhs_valleys, polys])
    plotter = plotter.plotter(graphs)
    plotter.plot(mode='hills')