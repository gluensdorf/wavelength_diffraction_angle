# !/usr/bin/env python3
# -*- coding: utf-8 -*-

import math as math
import numpy as np
import helper as hp

from wavelength_diffraction_angle import WavelengthDiffractionAngle as wlda

# data = hp.load_file('data/Testwerte 21-12-2017.txt')
data = hp.load_file('data/new_test')
value = hp.get_values_from_data(data)
result = [hp.get_header()]
num_iteration = int()

for pair in value:
    num_iteration = 0
    diag = pair[-2] / 1000
    vert = pair[-1] / 1000
    _wlda = wlda(
        g=10/13, m=1, h_z=0.769231, 
        lambda_min=0.4300, beta_Mb=np.radians(5)
    )

    _wlda.calc_diffraction_angle(lambda_vert=vert, lambda_diag=diag)
    _wlda.calc_x_0()
    _wlda.calc_jokabi_matrix()
    _wlda.determinant_jakobi()
    _wlda.calc_function_equation()
    _wlda.noname()
    _wlda.calc_tilting_angle()

    result.append(
        hp.format_result(
            _wlda=_wlda, 
            theta=pair[2], phi=pair[3], tau=pair[4], bx=pair[0], by=pair[1], 
            num_iteration=num_iteration, vert=vert, diag=diag, diff=False
        )
    )

hp.write_into_file("diff_output.txt", result)