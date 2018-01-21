# !/usr/bin/env python3
# -*- coding: utf-8 -*-

import math as math
import numpy as np

from wavelength_diffraction_angle import WavelengthDiffractionAngle as wlda

# -5.000 | -5.000 |   7.053 |  135.000 | -0.996 |    696.15911  |    696.40479  
# -2.500 | -2.500 |   3.533 |  135.000 | -0.500 |    630.49157  |    663.66575  
# g is in micro-meter
# lambda_min is in nano-meters
foo = wlda(
    g=10/13, m=1, h_z=0.769231, 
    lambda_min=0.4300, beta_Mb=np.radians(5)
    )

foo.calc_diffraction_angle(lambda_diag=0.69640479, lambda_vert=0.69615911) # first sample
# foo.calc_diffraction_angle(lambda_diag=0.63049157, lambda_vert=0.66366575)
foo.calc_x_0()
foo.calc_jokabi_matrix()
foo.determinant_jakobi()
foo.calc_function_equation()
foo.noname()

bar = foo.calc_tilting_angle()
print(bar)
a = np.degrees(bar[0])
b = np.degrees(bar[1])
print(a)
print(b)