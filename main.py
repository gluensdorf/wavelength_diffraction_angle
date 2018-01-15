#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import math as math
import numpy as np

class foobar:
    def __init__(self, g, m, hz, lambda_min, beta_Mb):
        """
        g is the raster constant
        m is the order)
        hz 'Abstand opt. Gitter - geom. Muster bzw. Pendel'
        lambda_min is wavelength minimum
        beta_Mb is the 'Messbereich'
        """
        self.g = g
        self.m = m
        self.h_z = h_z
        self.lambda_min = lambda_min
        self.beta_Mb = beta_Mb
        self.alpha_g = 2*beta_Mb # +/- 10°
        self.delta_lambda = (2 * g)/m * np.sin(alpha_g)
        self.lambda_max = delta_lambda + lambda_min
        self.alpha_e = np.arcsin(-(m/(2 * g))*(lambda_min + lambda_max))

        self.omega
        self.phi
        self.tau
        self.alpha_m
        self.A
        self.B
        self.C
        pass

    def calc_diffraction_angle(self):
        # we need to specify 'interval'
        # we propably need a loop for the interval - need more information on how it will be calculated
        interval = [self.lambda_min, self.lambda_max]
        self.alpha_m = np.arcsin((self.m * self.interval / self.g) + np.sin(self.alpha_e))
        # does alpha_m include the diag and vert version of it?

    def calc_x_0(self):
        a = 'what am I? how big am I?' # what value does 'a' have? 
        alpha_m_vert = self.alpha_m[0]
        alpha_m_diag = self.alpha_m[1]

        if alpha_m_vert == self.beta_Mb: 
            self.omega = -self.alpha_m_diag
            self.phi = 0
            self.tau = (-self.alpha_m_diag / self.beta_Mb) 
        else:
            self.omega = a * np.arctan(np.sqrt(math.pow(np.tan(
                    self.beta_Mb - alpha_m_vert), 2) + math.pow(np.tan(
                        alpha_m_vert - alpha_m_diag - self.beta_Mb), 2) 
                    )
                ),
            self.phi = np.arccos((a * np.tan(
                    alpha_m_vert - alpha_m_diag - self.beta_Mb)) / 
                    (np.sqrt(math.pow(np.tan(self.beta_Mb - alpa_m_vert)) + 
                        math.pow(np.tan(
                            alpha_m_vert - alpha_m_diag - self.beta_Mb), 2)
                        )
                    )
                ),
            self.tau = ((alpha_m_vert - alpha_m_diag) / self.beta_Mb) - 1

    def calc_jokabi_matrix(self):
        self.A = self.beta_Mb * np.sin(self.omega) * np.cos(self.beta_Mb * self.tau) * np.array(
            [
                np.sin(self.phi) * np.tan(self.beta_Mb * self.tau) * np.tan(self.beta_Mb) + np.cos(self.phi),
                -(np.sin(self.phi) + np.cos(self.phi)),
                np.cos(self.phi)*(1 - np.tan(self.beta_Mb * self.tau) * np.tan(self.beta_Mb))
            ]
        )
        self.B = self.beta_Mb * np.sin(self.omega) * np.cos(self.beta_Mb * self.tau) * np.array(
            [
                np.cos(self.phi) * np.tan(self.beta_Mb * self.tau) * np.tan(self.beta_Mb) - np.sin(self.phi) - np.tan(self.omega) * np.tan(self.alpha_m[0]),
                np.sin(self.phi) - np.cos(self.phi) + self.tan(self.omega) * self.tan(self.alpha_m[1]),
                np.tan(self.beta_Mb * self.tau) * np.tan(self.beta_Mb) * (np.tan(self.omega) * np.tan(self.alpha_m[1] + np.sin(self.phi)) - 
                    np.tan(self.omega) * np.tan(self.alpha_m[0] - np.sin(self.phi))
            ]
        self.C = math.pow(np.sin(self.omega), 2) * np.array(
            [
                np.cos(self.omega) + np.sin(self.phi) * self.tan(self.alpha_m[0]),
                -(np.cos(self.omega) + np.sin(self.phi) * np.tan(self.alpha_m[1])),
                np.cos(self.phi) * (np.tan(self.alpha_m[1]) - np.tan(self.alpha_m[0]))
            ]
        )

    def determinant_jakobi(self):
        # found no 'cot' function in numpy, used '1/tan' instead
        determinant = self.beta_Mb * np.cos(self.beta_Mb * self.tau) * math.pow(np.sin(self.omega), 2) * (
            np.cos(self.omega) * (np.tan(self.alpha_m[0]) - np.tan(alpha_m[1])) + np.sin(self.omega) * np.tan(self.alpha_m[0]) - 
            np.tan(self.beta_Mb * self.tau) * np.tan(self.beta_Mb) * (1/(np.tan(self.omega)) + np.sin(self.omega) * np.tan(self.alpha_m[1])) + 1/(np.tan(self.omega))
        )

    # need a better name for that function
    def calc_function_equation(self):
        func_equation = self.h_z * np.array(
            [
                np.cos(self.omega) * np.tan(self.alpha_m[1]) - np.sin(self.omega) * np.sin(self.phi),
                np.cos(self.omega) * np.tan(self.alpha_m[0]) - np.sin(self.omega) * np.sin(self.phi),
                np.sin(self.omega) * np.cos(self.phi)
            ]
        ) - self.h_z * np.array(
            [
                -np.sin(self.beta_Mb * self.tau),
                np.cos(self.beta_Mb * self.tau) * np.tan(self.beta_Mb),
                np.sin(self.beta_Mb * self.tau)
            ]
        )
        return func_equation
    
    def calc_tilting_angle(self):
        pass
    