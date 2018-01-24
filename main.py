# !/usr/bin/env python3
# -*- coding: utf-8 -*-

import math as math
import numpy as np

from wavelength_diffraction_angle import WavelengthDiffractionAngle as wlda

#  Sensorneigung  |   Transformationsparameter  |Wellenl�nge am |Wellenl�nge am
#   bx   |   by   |  Theta  |    Phi        t   |diagonalen Stab|vertikalen Stab
# [grad] | [grad] |  [grad] |   [grad]     []   |     [nm]      |    [nm]
# ------------------------------------------------------------------------------
# -5.000 | -5.000 |   7.053 |  135.000 | -0.996 |    696.15911  |    696.40479  
# -2.500 | -2.500 |   3.533 |  135.000 | -0.500 |    630.49157  |    663.66575  
# g is in micro-meter
# lambda_min is in nano-meters

f = open('data/Testwerte 21-12-2017.txt', 'r')
data = list(f)[13:len(list(f))-1]
f.close()
value = []
for x in data:
    x.split('|', 6)[-2:]
    value.append([float(bla) for bla in x.split('|', 6)])

for pair in value:
    diag = pair[-2] / 1000
    vert = pair[-1] / 1000
    foo = wlda(
        g=10/13, m=1, h_z=0.769231, 
        lambda_min=0.4300, beta_Mb=np.radians(5)
        )

    # foo.calc_diffraction_angle(lambda_vert=0.69640479, lambda_diag=0.69615911) # first sample
    foo.calc_diffraction_angle(lambda_vert=vert, lambda_diag=diag)
    foo.calc_x_0()
    foo.calc_jokabi_matrix()
    foo.determinant_jakobi()
    foo.calc_function_equation()
    foo.noname()

    bar = foo.calc_tilting_angle()
    print(f'lambda_vert: {vert}\nlambda_diag: {diag}')
    print(bar)
    a = np.degrees(bar[0])
    b = np.degrees(bar[1])
    print(a)
    print(b)