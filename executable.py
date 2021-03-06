# !/usr/bin/env python3
# -*- coding: utf-8 -*-

import math as math
import numpy as np
import sys

"""
TODO:
die main soll alle __init__ parameter uebergeben bekommen wenn die executable ausgefuehrt wird


nachfolgendes nicht in diesem script abfangen sondern in labview(?)

winkelpaare + fehlercode soll immer als ergebnis geliefert werden
    fehlercodes siehe labview
    das script kann an sich vordefinierte fehlercodes schmeissen und den bool (sagt obs ein fehler oder eine warnung ist)
        warnung: bool=FALSE, code != 0
        fehler: bool=TRUE, code != 0 (evtl. reicht es auch dass bool=TRUE damit es ein fehler ist)
        ok:     bool=FALSE, code == 0
    das labview programm muss nochmal erweitert werden um moegliche fehler richtig zu handhaben

winkelpaare sollen innerhalb von lambda_min und dem berechneten lambda_max liegen
    wenn nicht soll ein entsprechender fehlercode geworfen werden
    winkel sind beide auf 99grad

"""

class WavelengthDiffractionAngle:
    def __init__(self, g, m, h_z, lambda_min, beta_Mb):
        """
        g is the raster constant
        m is the order
        h_z 'Abstand opt. Gitter - geom. Muster bzw. Pendel'
        lambda_min is wavelength minimum
        beta_Mb is the 'Messbereich'
        """
        self.g = g
        self.m = m
        self.h_z = h_z
        self.lambda_min = lambda_min
        self.beta_Mb = beta_Mb
        self.alpha_g = 2 * self.beta_Mb # +/- 10°
        self.delta_lambda = (2 * self.g)/self.m * np.sin(self.alpha_g)
        self.lambda_max = self.delta_lambda + self.lambda_min
        self.alpha_e = np.arcsin(
            -1 * (self.m * (self.lambda_min + self.lambda_max)/(2 * self.g))
        )
        self.new_theta = 0.0
        self.new_phi = 0.0
        self.new_tau = 0.0
        self.beta_xy = []
        

    def get_transformation_parameters(self):
        """
        returns new_theta, new_phi and new_tau in a list
        """
        return [self.new_theta, self.new_phi, self.new_tau]

    def set_transformation_parameters(self, theta, phi, tau):
        """
        sets value self.theta, self.phi and self.tau to corresponding 
        given parameter value
        """
        self.theta = theta
        self.phi = phi
        self.tau = tau

    # takes lambda_vert and lambda_diag checks if they are within the range of lambda_min and lambda_max (correct?)
    # if so, calculates alpha_m_vert and alpha_m_diag
    # and returns alpha_m[alpha_m_vert, alpha_m_diag]
    def calc_diffraction_angle(self, lambda_diag, lambda_vert):
        self.alpha_m = np.array([
            np.arcsin((self.m * lambda_vert) / self.g + np.sin(self.alpha_e)),
            np.arcsin((self.m * lambda_diag) / self.g + np.sin(self.alpha_e))
        ])

    def calc_x_0(self):
        alpha_m_vert = self.alpha_m[0]
        alpha_m_diag = self.alpha_m[1]

        if alpha_m_vert > self.beta_Mb:
            a = 1
        else:
            a = -1

        # TODO: with the test values the if-case is never executed
        # if alpha_m_vert == self.beta_Mb: 
        if abs(alpha_m_vert - self.beta_Mb) <= 10**(-16): 
            self.theta = -alpha_m_diag
            self.phi = 0
            self.tau = -alpha_m_diag / self.beta_Mb 
            # print('if-case')
            # print('type(self.theta): ', type(self.theta))
            # print('type(self.phi): ', type(self.phi))
        else:
            self.theta = a * np.arctan(
                np.sqrt(
                    math.pow(
                        np.tan(
                            self.beta_Mb - alpha_m_vert), 2) + math.pow(np.tan(
                                alpha_m_vert - alpha_m_diag - self.beta_Mb), 2) 
                    )
                ),
            self.phi = np.arccos((a * np.tan(
                    alpha_m_vert - alpha_m_diag - self.beta_Mb)) / 
                    (np.sqrt(math.pow(np.tan(self.beta_Mb - alpha_m_vert), 2) + 
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
        self.B = self.beta_Mb * np.cos(self.theta) * np.cos(self.beta_Mb * self.tau) * np.array(
            [
                np.cos(self.phi) * np.tan(self.beta_Mb * self.tau) * np.tan(self.beta_Mb) - np.sin(self.phi) - np.tan(self.theta) * np.tan(self.alpha_m[0]),
                np.sin(self.phi) - np.cos(self.phi) + np.tan(self.theta) * np.tan(self.alpha_m[1]),
                np.tan(self.beta_Mb * self.tau) * np.tan(self.beta_Mb) * (np.tan(self.theta) * np.tan(self.alpha_m[1]) + np.sin(self.phi)) - 
                    np.tan(self.theta) * np.tan(self.alpha_m[0]) - np.sin(self.phi)
            ]
        )
        self.C = math.pow(np.sin(self.theta), 2) * np.array(
            [
                np.cos(self.theta) + np.sin(self.phi) * np.tan(self.alpha_m[0]),
                -(np.cos(self.theta) + np.sin(self.phi) * np.tan(self.alpha_m[1])),
                np.cos(self.phi) * (np.tan(self.alpha_m[1]) - np.tan(self.alpha_m[0]))
            ]
        )

    def determinant_jakobi(self):
        # found no 'cot' function in numpy, used '1/tan' instead
        self.determinant = self.beta_Mb * np.cos(self.beta_Mb * self.tau) * math.pow(np.sin(self.theta), 2) * ( #looks correct
            np.cos(self.phi) * (np.tan(self.alpha_m[0]) - np.tan(self.alpha_m[1])) + np.sin(self.phi) * np.tan(self.alpha_m[0]) - 
            np.tan(self.beta_Mb * self.tau) * np.tan(self.beta_Mb) * (1/(np.tan(self.theta)) + np.sin(self.phi) * np.tan(self.alpha_m[1])) + 1/(np.tan(self.theta))
        )

    # need a better name for that function
    def calc_function_equation(self):
        self.func_equation = np.array(
            [
                [np.cos(self.theta[0]) * np.tan(self.alpha_m[1]) - np.sin(self.theta[0]) * np.sin(self.phi[0])],
                [np.cos(self.theta[0]) * np.tan(self.alpha_m[0]) - np.sin(self.theta[0]) * np.sin(self.phi[0])],
                [np.sin(self.theta[0]) * np.cos(self.phi[0])]
            ]
        ) - np.array(
            [
                [-np.sin(self.beta_Mb * self.tau)],
                [np.cos(self.beta_Mb * self.tau) * np.tan(self.beta_Mb)],
                [np.sin(self.beta_Mb * self.tau)]
            ]
        )
    
    def noname(self):
        foo = np.array([[self.theta[0]], [self.phi[0]], [self.tau]]) # is 3x1
        ABC = np.hstack((self.A, self.B, self.C))
        foobar = np.transpose(ABC) / self.determinant
        bar = foobar.dot(self.func_equation) # is 3x3
        self.x_1 = foo - bar

    def calc_tilting_angle(self):
        self.new_theta = self.x_1[0][0]
        self.new_phi = self.x_1[1][0]
        self.new_tau = self.x_1[2][0]

        if self.new_theta == 0:
            self.beta_xy = [0, 0]
            return
        else:
            sign_x = -1 * np.sign(self.new_theta)
            if self.new_phi >= 0 and self.new_phi < np.pi/2:
                sign_y = -1 * sign_x
            else:
                sign_y = sign_x

            beta_x = sign_x * np.arccos(
                1/np.sqrt( 
                    1 + math.pow(np.tan(self.new_theta), 2) * math.pow(np.sin(self.new_phi), 2)
                    )
                )
            
            beta_y = sign_y * np.arccos( 
                1/np.sqrt( 
                    1 + math.pow(np.tan(self.new_theta), 2) * math.pow(np.cos(self.new_phi), 2)
                    )
                )
            self.beta_xy = np.array([beta_x, beta_y])

data_x = sys.argv[1]
data_y = sys.argv[2]

def main():
    diag = float(data_x)
    vert = float(data_y)
    _wlda = WavelengthDiffractionAngle(
        g=10/13, m=1, h_z=0.769231, 
        lambda_min=0.4300, beta_Mb=np.radians(5)
    )

    _wlda.calc_diffraction_angle(lambda_vert=vert/1000, lambda_diag=diag/1000)
    _wlda.calc_x_0()
    _wlda.calc_jokabi_matrix()
    _wlda.determinant_jakobi()
    _wlda.calc_function_equation()
    _wlda.noname()
    _wlda.calc_tilting_angle()

    trans_params = _wlda.get_transformation_parameters()

    # print(np.degrees(_wlda.beta_xy[0]))
    # print(np.degrees(_wlda.beta_xy[1]))
    return _wlda.beta_xy

main()