#!/usr/bin/env python3
#-*- coding: utf-8 -*-
"""
The purpose of this script is to read relevant variables from the output of eiktest and do further calculations using those variables.
"""

import numpy as np
from matplotlib import pyplot as plt
import pdb
from scipy.integrate import cumtrapz as ctrap


arr1 = np.loadtxt('./balls/eik5.out', skiprows = 1)
arr2 = np.loadtxt('./balls/eik6.out', skiprows = 3)



theta_col = arr1[:, 0]
aprime = arr2[:, 6]
#aprime = aprime[theta_col >= 0]

#theta_col = theta_col[theta_col >= 0]


bmag  = arr1[:, 2]

B_p = arr2[:, 7]


R = arr2[:, 1]
Z = arr2[:, 2]

dl = np.sqrt(derm(R, ch='l', par ='e')**2 + derm(Z, ch='l', par='o')**2)
L = np.cumsum(dl)

dtdl = B_p/bmag

intgrl = ctrap(dtdl, L, initial=0)
intgrl = intgrl - (intgrl[-1]/2)
theta_eqarc = (intgrl/intgrl[-1])*np.pi

##plt.plot(theta_col[theta_col >=0 ], arr1[:, 6][theta_col >= 0]); plt.show() 

#plt.plot(theta_col, local_shear)
local_shear = B_p*derm(aprime, ch ='l', par = 'e')/dl
pdb.set_trace()
