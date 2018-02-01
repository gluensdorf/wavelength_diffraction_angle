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

def write_output(_wlda, vert, diag, theta, phi, tau, bx, by, diff=False):
    if diff:
        result = f'{abs(bx - np.degrees(_wlda.beta_xy[0])):^+7.3f}|'\
            f'{abs(by - np.degrees(_wlda.beta_xy[1])):^+8.3f}|'\
            f'{abs(theta - np.degrees(_wlda.new_theta)):^+9.3f}|'\
            f'{abs(phi - np.degrees(_wlda.new_phi)):^+10.3f}|'\
            f'{abs(tau - _wlda.new_tau):^+8.3f}|'\
            f'{vert*1000:15.5f}|'\
            f'{diag*1000:15.5f}|'\
            f'{num_iteration:^10}\n'
    else:
        result = f'{np.degrees(_wlda.beta_xy[0]):^+7.3f}|'\
            f'{np.degrees(_wlda.beta_xy[1]):^+8.3f}|'\
            f'{np.degrees(_wlda.new_theta):^+9.3f}|'\
            f'{np.degrees(_wlda.new_phi):^+10.3f}|'\
            f'{_wlda.new_tau:^+8.3f}|'\
            f'{vert*1000:15.5f}|'\
            f'{diag*1000:15.5f}|'\
            f'{num_iteration:^10}\n'
    return result

'''
def write_diff_output(_wlda, vert, diag, theta, phi, tau, bx, by):
    result = f'{np.degrees(_wlda.beta_xy[0]):^+7.3f}|'\
        f'{np.degrees(_wlda.beta_xy[1]):^+8.3f}|'\
        f'{np.degrees(_wlda.new_theta):^+9.3f}|'\
        f'{np.degrees(_wlda.new_phi):^+10.3f}|'\
        f'{_wlda.new_tau:^+8.3f}|'\
        f'{vert*1000:15.5f}|'\
        f'{diag*1000:15.5f}|'\
        f'{num_iteration:^10}\n'
    return result
'''

f = open('data/Testwerte 21-12-2017.txt', 'r')
data = list(f)[13:len(list(f))-1]
f.close()
value = []
header = ' Sensorneigung  |   Transformationsparameter  |Wellenl�nge am |Wellenl�nge am |           \n'\
    '  bx   |   by   |  Theta  |    Phi        t   |diagonalen Stab|vertikalen Stab|Iteration  \n'\
    '[grad] | [grad] |  [grad] |   [grad]     []   |     [nm]      |    [nm]       |   []      \n'\
    '------------------------------------------------------------------------------------------\n'
result = [header]
num_iteration = int()

for x in data:
    x.split('|', 6)[-2:]
    value.append([float(bla) for bla in x.split('|', 6)])

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

    bar = _wlda.calc_tilting_angle()

    result.append(
        write_output(
            _wlda, vert, diag, pair[2], pair[3], 
            pair[4], pair[0], pair[1], diff=True
        )
    )
    """ BEGIN NEXT ITERATION """
    num_iteration = 1
    trans_params = _wlda.get_transformation_parameters()

    _wlda.calc_diffraction_angle(lambda_vert=vert, lambda_diag=diag)
    _wlda.set_transformation_parameters(
        theta=[trans_params[0]],
        phi=[trans_params[1]],
        tau=trans_params[2]
        )
    _wlda.calc_jokabi_matrix()
    _wlda.determinant_jakobi()
    _wlda.calc_function_equation()
    _wlda.noname()

    bar = _wlda.calc_tilting_angle()

    result.append(
        write_output(
            _wlda, vert, diag, pair[2], pair[3], 
            pair[4], pair[0], pair[1], diff=True
        )
    )

with open("foo.txt", "a") as myfile:
    for line in result:
        myfile.write(line)
myfile.close()
"""
with open("output.txt", "a") as myfile:
    for line in result:
        myfile.write(line)
myfile.close()
"""