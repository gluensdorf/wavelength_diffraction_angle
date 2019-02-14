# !/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Button
from pprint import pprint as pp


class plotter:
    def __init__(self, data):
        self.data = data
        self.index = 0
        self.data_sample = data[self.index]

    def __execute_mode(self):
        if self.mode == 'hills':
            self.update_hills()

    def __check_mode(self):
        if self.mode == 'hills':
            self.prepare_hills()
        else:
            raise Exception('mode is missing')

    def __next(self, event):
        self.index += 1
        self.index = self.index % len(self.data)
        self.data_sample = self.data[self.index]
        self.__execute_mode()
        self.__print_poly_maximum()
        plt.draw()

    def __prev(self, event):
        self.index -= 1
        self.index = self.index % len(self.data)
        self.data_sample = self.data[self.index]
        self.__execute_mode()
        self.__print_poly_maximum()
        plt.draw()

    def __on_key(self, event):
        if event.key == 'left':
            self.__prev(event)
        elif event.key == 'right':
            self.__next(event)

    def __print_poly_maximum(self):
        for idx, poly in enumerate(self.data[self.index][-1]):
            x_lhs_valleys = self.data[self.index][3][0][idx]
            x_rhs_valleys = self.data[self.index][4][0][idx]
            derivate = np.polyder(poly).roots

            r_derivate = derivate[derivate.imag == 0].real
            crossings = [x for x in r_derivate if x >=
                         x_lhs_valleys and x <= x_rhs_valleys]
            # pp(f"{self.index} : {crossings}")
            poly_peak = None
            for idx, val in enumerate(crossings):
                if poly(val) == max(poly(crossings)):
                    poly_peak = val
            pp(f"{self.index} maximum at: {poly_peak}")

            # # derivate = poly.deriv().r
            # r_derivate = derivate[derivate.imag==0].real
            # test = poly.deriv(2)(r_derivate)
            # maximum = max(poly(r_derivate))
            # print(f"{self.index} : {r_derivate[test<0]}")
            # print(f"maximum : {maximum}")
            # print(f"{self.index:>4} : {max(p_der.roots)}")
            # print(f"{self.index:>4} : {max(poly(p_der.roots))}")

    def plot(self, mode):
        self.mode = mode
        self.figure = plt.figure()
        self.figure.subplots_adjust(bottom=0.2)
        self.ax = self.figure.add_subplot(111)
        self.__check_mode()

        axprev = plt.axes([0.7, 0.05, 0.1, 0.075])
        axnext = plt.axes([0.81, 0.05, 0.1, 0.075])
        bnext = Button(axnext, 'Next')
        bnext.on_clicked(self.__next)
        bprev = Button(axprev, 'Previous')
        bprev.on_clicked(self.__prev)
        self.figure.canvas.mpl_connect('key_press_event', self.__on_key)

        plt.show()

    def update_hills(self):
        spectrum = self.data_sample[0]
        x_range = self.data_sample[1]
        peaks = self.data_sample[2]
        lhs_valleys = self.data_sample[3]
        rhs_valleys = self.data_sample[4]
        polys = self.data_sample[5]

        self.line1.set_xdata(spectrum[0])
        self.line1.set_ydata(spectrum[1])
        self.line2.set_xdata(spectrum[0][peaks])
        self.line2.set_ydata(spectrum[1][peaks])
        self.line3.set_xdata(lhs_valleys[0])
        self.line3.set_ydata(lhs_valleys[1])
        self.line4.set_xdata(rhs_valleys[0])
        self.line4.set_ydata(rhs_valleys[1])
        self.line5.set_xdata(x_range[0])
        self.line5.set_ydata(polys[0](x_range[0]))
        self.line6.set_xdata(x_range[1])
        self.line6.set_ydata(polys[1](x_range[1]))

        if len(x_range) >= 3:
            self.line7.set_xdata(x_range[2])
            self.line7.set_ydata(polys[2](x_range[2]))
        else:  # case if only two peaks had been found
            self.line7.set_xdata(x_range[1])
            self.line7.set_ydata(polys[1](x_range[1]))

    def prepare_hills(self):
        """
        Plots spectrum, fitted polynomial, valleys and peaks

        # :param spectrum: a list of [x,y] values
        # :param x_range: list of artificial samplepoints between element of 
        #                 lhs_valleys and corresponding element of rhs_valleys
        # :param peaks: list of x-values of peaks in spectrum
        # :param lhs_valleys: list of [x,y] values of a valley, order corresponds 
        #                     to peaks ordering
        # :param rhs_valleys: list of [x,y] values of a valley, order corresponds 
        #                     to peaks ordering
        # :param polys: list of polynomials, polynomials are fitted using the 
        #             samplepoints between element of lhs_valleys and 
        #             corresponding element of rhs_valleys
        """
        spectrum = self.data_sample[0]
        x_range = self.data_sample[1]
        peaks = self.data_sample[2]
        lhs_valleys = self.data_sample[3]
        rhs_valleys = self.data_sample[4]
        polys = self.data_sample[5]

        self.line1, = self.ax.plot(spectrum[0], spectrum[1], "r-")
        self.line2, = self.ax.plot(
            spectrum[0][peaks], spectrum[1][peaks], "gx")
        self.line3, = self.ax.plot(lhs_valleys[0], lhs_valleys[1], "go")
        self.line4, = self.ax.plot(rhs_valleys[0], rhs_valleys[1], "ro")
        self.line5, = self.ax.plot(x_range[0], polys[0](x_range[0]), "b-")
        self.line6, = self.ax.plot(x_range[1], polys[1](x_range[1]), "b-")

        if len(x_range) >= 3:
            self.line7, = self.ax.plot(x_range[2], polys[2](x_range[2]), "b-")
        else:  # case if only two peaks had been found
            self.line7, = self.ax.plot(x_range[1], polys[1](x_range[1]), "b-")
