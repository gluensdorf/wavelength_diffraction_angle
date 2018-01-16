!/usr/bin/env python3
# -*- coding: utf-8 -*-

import math as math
import numpy as np

class foobar:
    def __init__(self, g, m, hz, lambda_min, beta_Mb):
        """
        g is the raster constant
        m is the order
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

        self.theta
        self.phi
        self.tau
        self.alpha_m
        self.A
        self.B
        self.C
        self.determinant
        self.x_1

    # takes lambda_vert and lambda_diag checks if they are within the range of lambda_min and lambda_max (correct?)
    # if so, calculates alpha_m_vert and alpha_m_diag
    # and returns alpha_m[alpha_m_vert, alpha_m_diag]
    def calc_diffraction_angle(self, lambda_vert, lambda_diag):
        # TODO: insert check for lambda_min and lambda_max
        self.alpha_m = [
            np.arcsin((self.m * lambda_vert / self.g) + np.sin(self.alpha_e)),
            np.arcsin((self.m * lambda_diag / self.g) + np.sin(self.alpha_e))
        ]

    def calc_x_0(self):
        alpha_m_vert = self.alpha_m[0]
        alpha_m_diag = self.alpha_m[1]

        if alpha_m_vert > self.beta_Mb:
            a = 1
        else:
            a = -1

        if alpha_m_vert == self.beta_Mb: 
            self.theta = -alpha_m_diag
            self.phi = 0
            self.tau = (-alpha_m_diag / self.beta_Mb) 
        else:
            self.theta = a * np.arctan(np.sqrt(math.pow(np.tan(
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
        self.A = self.beta_Mb * np.sin(self.theta) * np.cos(self.beta_Mb * self.tau) * np.array(
            [
                np.sin(self.phi) * np.tan(self.beta_Mb * self.tau) * np.tan(self.beta_Mb) + np.cos(self.phi),
                -(np.sin(self.phi) + np.cos(self.phi)),
                np.cos(self.phi)*(1 - np.tan(self.beta_Mb * self.tau) * np.tan(self.beta_Mb))
            ]
        )
        self.B = self.beta_Mb * np.sin(self.theta) * np.cos(self.beta_Mb * self.tau) * np.array(
            [
                np.cos(self.phi) * np.tan(self.beta_Mb * self.tau) * np.tan(self.beta_Mb) - np.sin(self.phi) - np.tan(self.theta) * np.tan(self.alpha_m[0]),
                np.sin(self.phi) - np.cos(self.phi) + self.tan(self.theta) * self.tan(self.alpha_m[1]),
                np.tan(self.beta_Mb * self.tau) * np.tan(self.beta_Mb) * (np.tan(self.theta) * np.tan(self.alpha_m[1] + np.sin(self.phi)) - 
                    np.tan(self.theta) * np.tan(self.alpha_m[0] - np.sin(self.phi))
            ]
        self.C = math.pow(np.sin(self.theta), 2) * np.array(
            [
                np.cos(self.theta) + np.sin(self.phi) * self.tan(self.alpha_m[0]),
                -(np.cos(self.theta) + np.sin(self.phi) * np.tan(self.alpha_m[1])),
                np.cos(self.phi) * (np.tan(self.alpha_m[1]) - np.tan(self.alpha_m[0]))
            ]
        )

    def determinant_jakobi(self):
        # found no 'cot' function in numpy, used '1/tan' instead
        self.determinant = self.beta_Mb * np.cos(self.beta_Mb * self.tau) * math.pow(np.sin(self.theta), 2) * (
            np.cos(self.theta) * (np.tan(self.alpha_m[0]) - np.tan(alpha_m[1])) + np.sin(self.theta) * np.tan(self.alpha_m[0]) - 
            np.tan(self.beta_Mb * self.tau) * np.tan(self.beta_Mb) * (1/(np.tan(self.theta)) + np.sin(self.theta) * np.tan(self.alpha_m[1])) + 1/(np.tan(self.theta))
        )

    # need a better name for that function
    def calc_function_equation(self):
        func_equation = self.h_z * np.array(
            [
                np.cos(self.theta) * np.tan(self.alpha_m[1]) - np.sin(self.theta) * np.sin(self.phi),
                np.cos(self.theta) * np.tan(self.alpha_m[0]) - np.sin(self.theta) * np.sin(self.phi),
                np.sin(self.theta) * np.cos(self.phi)
            ]
        ) - self.h_z * np.array(
            [
                -np.sin(self.beta_Mb * self.tau),
                np.cos(self.beta_Mb * self.tau) * np.tan(self.beta_Mb),
                np.sin(self.beta_Mb * self.tau)
            ]
        )
        return func_equation
    
    def noname(self):
        self.x_1 = np.array([self.theta, self.phi, self.tau]) - (
            (np.array([self.A, self.B, self.C]) / self.determinant) * self.calc_function_equation()
        )

    def calc_tilting_angle(self):
        new_theta = self.x_1[0]
        new_phi = self.x_1[1]
        new_tau = self.x_1[2]

        if :
            beta_xy = [0, 0]
        else if new_phi >= 0 and new_phi < np.pi/2:
            beta_xy = np.sign(new_theta) * (
                np.array(
                    [
                        -np.arccos(
                            1/np.sqrt(
                                1 + math.pow(np.tan(new_theta), 2) * math.pow(np.sin(new_phi), 2)
                            )
                        ),
                        np.arccos(
                            1/np.sqrt(
                                1 + math.pow(np.tan(new_theta), 2) * math.pow(np.cos(new_phi), 2)
                            )
                        )
                    ]
                )
            )
        else if new_phi >= np.pi/2 and new_phi < np.pi:
            beta_xy = np.sign(new_theta) * ( # missing '-' infront of np.sign()?
                np.array(
                    [
                        np.arccos(
                            1/np.sqrt(
                                1 + math.pow(np.tan(new_theta), 2) * math.pow(np.sin(new_phi), 2)
                            )
                        ),
                        np.arccos(
                            1/np.sqrt(
                                1 + math.pow(np.tan(new_theta), 2) * math.pow(np.cos(new_phi), 2)
                            )
                        )
                    ]
                )
            )
    