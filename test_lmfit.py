# !/usr/bin/env python3
# -*- coding: utf-8 -*-

import math as math
import numpy as np
import helper as hp
import matplotlib.pyplot as plt
# from lmfit.models import PseudoVoigtModel as PVM
from lmfit import models

from wavelength_diffraction_angle import WavelengthDiffractionAngle as wlda

# data = hp.load_file('data/make_params')
data = np.loadtxt('data/07_Serie_Prototyp/13-N-0Grad2.SSM', skiprows=2, unpack=True)
# value = hp.get_values_from_data(data)
x = data[0]
y = data[1]

# for record in value:
    # x.append(record[-2])
    # y.append(record[-1])prefixes

# x = np.array(x)
# y = np.array(y)

# print(x)
"""
x = np.array([
    696.15911, 689.68461, 683.19140, 676.67985, 670.15033,
    663.60325, 657.03900, 650.45797, 643.86058, 637.24725,
    630.61840, 623.97447, 617.31589, 610.64312, 603.95661,
    597.25682, 590.54423, 583.81931, 577.08256, 570.33446,
    563.57552 #, 689.68461, 683.18245, 676.66290, 670.12634,
    # 663.57317, 657.00380, 650.41865, 643.81813, 637.20267,
    # 630.57270, 623.92867, 617.27102, 610.60021, 603.91669,
    # 597.22094, 590.51342, 583.79463, 577.06505, 570.32518,
    # 563.57552, 556.81658 
])
y = np.array([
    696.40479, 696.40479, 696.40479, 696.40479, 696.40479,
    696.40479, 696.40479, 696.40479, 696.40479, 696.40479,
    696.40479, 696.40479, 696.40479, 696.40479, 696.40479,
    696.40479, 696.40479, 696.40479, 696.40479, 696.40479,
    696.40479 #, 689.88437, 689.88437, 689.88437, 689.88437,
    # 689.88437, 689.88437, 689.88437, 689.88437, 689.88437,
    # 689.88437, 689.88437, 689.88437, 689.88437, 689.88437,
    # 689.88437, 689.88437, 689.88437, 689.88437, 689.88437,
    # 689.88437, 689.88437
])
"""
# want to find 3 peaks and assign them a prefix each
peak1 = models.PseudoVoigtModel(prefix='p1_') # (nan_policy='propagate')
peak2 = models.PseudoVoigtModel(prefix='p2_')
peak3 = models.PseudoVoigtModel(prefix='p3_')
model = peak1 + peak2 + peak3
params = model.make_params(p1_sigma=0.701, p2_sigma=0.602, p3_sigma=0.803)
params = peak1.guess(y, x=x)
params.update(peak2.guess(y, x=x))
params.update(peak3.guess(y, x=x))

params['p1_center'].set(min=350, max=800)
params['p2_center'].set(min=350, max=800)
params['p3_center'].set(min=350, max=800)
result = model.fit(y, params, x=x)
# print(result.fit_report())

print(result.fit_report(min_correl=0.25))
result.plot()

plt.plot(x, y, 'bo')
plt.plot(x, result.init_fit, 'k--')
plt.plot(x, result.best_fit, 'r-')
plt.show()
