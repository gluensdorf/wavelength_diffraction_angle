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

# single_wlda = hp.single_calculation(563.57552/1000, 630.61840/1000)
# print(np.degrees(single_wlda.phi))

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
            _wlda=_wlda, #vert, diag, 
            theta=pair[2], phi=pair[3], tau=pair[4], bx=pair[0], by=pair[1], 
            num_iteration=num_iteration, vert=vert, diag=diag, diff=False
        )
    )
    """
    # BEGIN NEXT ITERATION 
    num_iteration = 1
    trans_params = _wlda.get_transformation_parameters()

    # _wlda.calc_diffraction_angle(lambda_vert=vert, lambda_diag=diag)
    _wlda.set_transformation_parameters(
        theta=(trans_params[0],),
        phi=(trans_params[1],),
        tau=trans_params[2]
        )
    _wlda.calc_jokabi_matrix()
    _wlda.determinant_jakobi()
    _wlda.calc_function_equation()
    _wlda.noname()
    _wlda.calc_tilting_angle()

    result.append(
        hp.format_result(
            _wlda=_wlda, #vert, diag, 
            theta=pair[2], phi=pair[3], tau=pair[4], 
            bx=pair[0], by=pair[1], 
            num_iteration=num_iteration, vert=vert, diag=diag, diff=False
        )
    )

    num_iteration = 2
    trans_params = _wlda.get_transformation_parameters()

    # _wlda.calc_diffraction_angle(lambda_vert=vert, lambda_diag=diag)
    _wlda.set_transformation_parameters(
        theta=(trans_params[0],),
        phi=(trans_params[1],),
        tau=trans_params[2]
        )
    _wlda.calc_jokabi_matrix()
    _wlda.determinant_jakobi()
    _wlda.calc_function_equation()
    _wlda.noname()
    _wlda.calc_tilting_angle()

    result.append(
        hp.format_result(
            _wlda=_wlda, #vert, diag, 
            theta=pair[2], phi=pair[3], tau=pair[4], 
            bx=pair[0], by=pair[1], 
            num_iteration=num_iteration, vert=vert, diag=diag, diff=False
        )
    )
    """

hp.write_into_file("diff_output.txt", result)