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
    self.alpha_m = np.arcsin(self.m * self.)
    pass