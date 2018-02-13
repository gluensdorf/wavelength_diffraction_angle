# !/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from wavelength_diffraction_angle import WavelengthDiffractionAngle as wlda

def load_file(path_to_file):
    f = open(path_to_file, 'r')
    data = list(f)[13:len(list(f))-1]
    f.close()
    return data

def get_values_from_data(data):
    value = []
    for x in data:
        x.split('|', 6)[-2:]
        value.append([float(bla) for bla in x.split('|', 6)])
    return value

# change header to:
# bx by given | theta phi t | bx by calculated
# labview python integration, install labview
# fix script, last step still an issue?
# use new testvalues - still need to get it from the email
# adopted get_header() for new formating style
#
# THETA PHI and TAU are calculated with too much error

def get_header():
    header = ' Sensorneigung  |   Transformationsparameter  |\n'\
    ' g bx  | g by   |  Theta  |    Phi        t   | calc bx       | calc by       |Iteration  \n'\
    '[grad] | [grad] |  [grad] |   [grad]     []   |    [grad]     |    [grad]     |   []      \n'\
    '------------------------------------------------------------------------------------------\n'
    return header

def write_into_file(path_to_file, result):
    with open(path_to_file, "a") as myfile:
        for line in result:
            myfile.write(line)
    myfile.close()

def single_calculation(diag, vert):
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
    return _wlda

def single_calculation_two_iterations(diag, vert):
    num_iteration = int()
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

    trans_params = _wlda.get_transformation_parameters()

    result = []
    result.append(
        format_result(
            _wlda=_wlda, #vert, diag, 
            theta=trans_params[0], phi=trans_params[1], tau=trans_params[2], bx=_wlda.beta_xy[0], by=_wlda.beta_xy[1], 
            num_iteration=num_iteration, vert=vert, diag=diag, diff=False
        )
    )
    # BEGIN NEXT ITERATION 
    print('######### NEXT ITERATION ###########')
    num_iteration = 1
    _wlda = wlda(
        g=10/13, m=1, h_z=0.769231, 
        lambda_min=0.4300, beta_Mb=np.radians(5)
    )

    _wlda.calc_diffraction_angle(lambda_vert=vert, lambda_diag=diag)
    _wlda.set_transformation_parameters(
        theta=trans_params[0],
        phi=trans_params[1],
        tau=trans_params[2]
        )
    _wlda.calc_jokabi_matrix()
    _wlda.determinant_jakobi()
    _wlda.calc_function_equation()
    _wlda.noname()
    _wlda.calc_tilting_angle()

    trans_params = _wlda.get_transformation_parameters()
    result.append(
        format_result(
            _wlda=_wlda, #vert, diag, 
            theta=trans_params[0], phi=trans_params[1], tau=trans_params[2], bx=_wlda.beta_xy[0], by=_wlda.beta_xy[1], 
            num_iteration=num_iteration, vert=vert, diag=diag, diff=False
        )
    )

    write_into_file("diff_output.txt", result)
    return _wlda

def format_result(_wlda, theta, phi, tau, bx, by, num_iteration,vert, diag, diff=False): #vert, diag, 
    width = 20
    precision = 10
    if diff:
        result = f'{bx:^+{5}.{3}f}|'\
            f'{by:^+{5}.{3}f}|'\
            f'{abs(theta - np.degrees(_wlda.new_theta[0])):^+{width}.{precision}f}|'\
            f'{abs(phi - np.degrees(_wlda.new_phi[0])):^+{width}.{precision}f}|'\
            f'{abs(tau - _wlda.new_tau):^+{width}.{precision}f}|'\
            f'{abs(bx - np.degrees(_wlda.beta_xy[0])):^+{width}.{precision}f}|'\
            f'{abs(by - np.degrees(_wlda.beta_xy[1])):^+{width}.{precision}f}|'\
            f'{num_iteration:^10}\n'
            # f'{diag*1000:15.5f}|'\
            # f'{vert*1000:15.5f}|'\
    else:
        result = f'{bx:^+{5}.{3}f}|'\
            f'{by:^+{5}.{3}f}|'\
            f'{np.degrees(_wlda.theta[0]):^+{width}.{precision}f}|'\
            f'{np.degrees(_wlda.phi[0]):^+{width}.{precision}f}|'\
            f'{_wlda.tau:^+{width}.{precision}f}|'\
            f'{np.degrees(_wlda.beta_xy[0]):^+{width}.{precision}f}|'\
            f'{np.degrees(_wlda.beta_xy[1]):^+{width}.{precision}f}|'\
            f'{num_iteration:^10}\n'
            # f'{diag*1000:15.5f}|'\
            # f'{vert*1000:15.5f}|'\
    return result
