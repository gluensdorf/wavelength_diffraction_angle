#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import math as math
import numpy as np

def class __init__(self, g, m, hz, lambda_min, beta_Mb):
    """
    g is the raster constant
    m is the order)
    hz 'Abstand opt. Gitter - geom. Muster bzw. Pendel'
    lambda_min is wavelength minimum
    beta_Mb is the 'Messbereich'
    """
    self.g = g
    self.m = m
    self.hz = hz
    self.lambda_min = lambda_min
    self.beta_Mb = beta_Mb
    self.alpha_g = 2*beta_Mb # +/- 10°
    self.delta_lambda = (2 * g)/m * np.sin(alpha_g)
    self.lambda_max = delta_lambda + lambda_min
    self.alpha_e = np.arcsin(-(m/(2 * g))*(lambda_min + lambda_max))
    pass

def calc_diffraction_angle(self):
    # we need to specify 'interval'
    # we propably need a loop for the interval - need more information on how it will be calculated
    interval = [self.lambda_min, self.lambda_max]
    self.alpha_m = np.arcsin((self.m * self.interval / self.g) + np.sin(self.alpha_e))
    # alpha_m includes diag and vert version of it
    return alpha_m

def calc_x_0(self):
    a = 'what am I? how big am I?' # what value does 'a' have? 
    alpha_m_vert = alpha_m[0]
    alpha_m_diag = alpha_m[1]

    if alpha_m_vert == self.beta_Mb: 
        x_0 = [-self.alpha_m_diag, 0, (-self.alpha_m_diag / self.beta_Mb)]
    else:
        x_0 = [
            a * np.arctan(np.sqrt(math.pow(np.tan(self.beta_Mb - alpha_m_vert), 2) + math.pow(np.tan(alpha_m_vert - alpha_m_diag - self.beta_Mb), 2) )),
            np.arccos((a * np.tan(alpha_m_vert - alpha_m_diag - self.beta_Mb)) / (np.sqrt(math.pow(np.tan(self.beta_Mb - alpa_m_vert)) + 
                math.pow(np.tan(alpha_m_vert - alpha_m_diag - self.beta_Mb), 2)))),
            ((alpha_m_vert - alpha_m_diag) / self.beta_Mb) - 1
        ]
    return x_0