!/usr/bin/env python3
# -*- coding: utf-8 -*-

import math as math
import numpy as np

from wavelength_diffraction_angle import WavelengthDiffractionAngle as wlda

foo = wlda(
    g=0.769231, m=1, hz=80, 
    lambda_min=430.0000, beta_Mb=5
    )

foo.calc_diffraction_angle(lambda_vert=696.15911, lambda_diag=696.40479)
foo.calc_x_0()
foo.calc_jokabi_matrix()
foo.determinant_jakobi()
foo.calc_function_equation()
foo.noname()

bar = foo.calc_tilting_angle()
print(bar)