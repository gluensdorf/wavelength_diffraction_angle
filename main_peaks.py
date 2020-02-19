# !/usr/bin/env python3
# -*- coding: utf-8 -*-

# TODO: der link sollte hier spaeter raus
# link fuer einen vergleich von algorithmen?
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2631518/

import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import glob
import sys
from pprint import pprint as pp
from scipy.signal import savgol_filter

import plotter as plotter
import test_lmfit as lmfit


def write_into_file(path_to_file, result):
    with open(path_to_file, "a") as myfile:
        for line in result:
            myfile.write(line)
    myfile.close()


def prepare_results(graphs, paths_dataset):
    width = 14
    precision = 10
    # header = f"  id | poly-id | peak(x)        [nm] | filename \n"
    # breakline = f"------------------------------------------------\n"
    header = f" 1. peak [nm] | 2. peak [nm] | 3. peak [nm] | filename \n"
    breakline = f"------------------------------------------------\n"
    result = []
    result.append(header)
    result.append(breakline)
    # iterate over every measurement (a graph)
    for graph_idx, graph in enumerate(graphs):
        poly_peaks = []
        # iterate over every created polynom of a graph
        for idx, poly in enumerate(graph[-1]):
            x_lhs_valleys = graphs[graph_idx][3][0][idx]
            x_rhs_valleys = graphs[graph_idx][4][0][idx]
            derivate = np.polyder(poly).roots

            r_derivate = derivate[derivate.imag == 0].real
            crossings = [x for x in r_derivate if x >=
                         x_lhs_valleys and x <= x_rhs_valleys]
            poly_peak = None
            for _, val in enumerate(crossings):
                if poly(val) == max(poly(crossings)):
                    poly_peak = val
                    poly_peaks.append(poly_peak)
            # pp(f"{graph_idx} maximum at: {poly_peak}")

            # line = f"{graph_idx:4} | {idx:^6} |{poly_peak:^{width}.{precision}}|"\
            #       f" {paths_dataset[graph_idx].split('/')[-1]}\n"
            # result.append(line)
        poly_peaks.sort()
        if len(poly_peaks) == 3:
            line = f"{poly_peaks[0]:^{width}.{precision}}|"\
                   f"{poly_peaks[1]:^{width}.{precision}}|"\
                   f"{poly_peaks[2]:^{width}.{precision}}|"
        else:
            line = f"{'NONE':^{width}.{precision}}|"\
                   f"{poly_peaks[0]:^{width}.{precision}}|"\
                   f"{poly_peaks[1]:^{width}.{precision}}|"
        line += f" {paths_dataset[graph_idx].split('/')[-1]}\n"
        result.append(line)
        # result.append(breakline)
    return result


if __name__ == "__main__":
    # return list of paths of all files in folder
    paths_dataset = glob.glob(
        '/home/darlokh/Documents/hiwi/code/data/Messwerte_clean/2018*')

    graphs = []
    for idx, path in enumerate(paths_dataset):
        if 0 == idx % 50:
            print(idx)
        # skiprows set to 0 because data was cleaned
        data = np.loadtxt(path, skiprows=0, unpack=True)
        smoothed_signal = savgol_filter(data, 301, 2)  # numpy array

        x = data[0]  # wavelength
        y = data[1]  # amplitude
        smooth_y = smoothed_signal[1]
        order_value = 50
        peaks, _ = sp.signal.find_peaks(
            smooth_y, distance=order_value, height=0.003)
        # valleys = sp.signal.argrelmin(smooth_y, order=order_value)

        distance = 0
        hills = lmfit.isolate_hills(smoothed_signal, peaks, distance)

        lhs_valleys, rhs_valleys, hills_samplepoints = lmfit.select_valleys(
            smoothed_signal, hills)

        # outsource loop to a function
        polys = []
        x_range = []
        for samplepoints in hills_samplepoints:
            polys.append(
                np.poly1d(
                    np.polyfit(x=samplepoints[0], y=samplepoints[1], deg=5)))
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
    write_into_file('3-auswertung_2018-11-30.txt', prepare_results(graphs, paths_dataset))
