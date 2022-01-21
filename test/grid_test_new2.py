#!/usr/bin/env python3
#python script to compare various eikcoefs with the GS2 output for the given set of theta_grid_knob variables
#ntheta = 128
#nperiod = 2
#rhoc = 0.454
#qinp = 0.961
#shift = -0.060
#s_hat_input = 0.164
#R_geo = 0.7297266  # GS2 doesn't take Bunit so R_geo is taken from the local eikcoefs gen script
#Rmaj =  3.24564035 #1.850/a_N
#akappa = 1.398
#akappri = 0.0102
#tri_gs2 = -0.057
#tripri_gs2 = 0.026
#beta_prime_input = -0.1845 

import pdb
import os
import sys
import numpy as np
from netCDF4 import Dataset
from matplotlib import pyplot as plt


parnt_dir_name = os.path.dirname(os.getcwd())

fname_in_txt2 = '{0}/{1}/{2}'.format(parnt_dir_name, 'test', 'gs2_test_input.out.nc') 

rtg = Dataset(fname_in_txt2, 'r')
# read from eiktest output
tgrid     = rtg.variables['theta'][:].data
gbdrift   = rtg.variables['gbdrift'][:].data
gradpar   = rtg.variables['gradpar'][:].data
grho      = rtg.variables['grho'][:].data
cvdrift   = rtg.variables['cvdrift'][:].data
gds2      = rtg.variables['gds2'][:].data
bmag      = rtg.variables['bmag'][:].data
gds21     = rtg.variables['gds21'][:].data
gds22     = rtg.variables['gds22'][:].data
cvdrift0  =  rtg.variables['cvdrift0'][:].data
gbdrift0  = rtg.variables['gbdrift0'][:].data
aprime    = rtg.variables['aprime'][:].data
jacob     = rtg.variables['jacob'][:].data
Rplot     = rtg.variables['Rplot'][:].data
Zplot     = rtg.variables['Zplot'][:].data

pdb.set_trace()

plt.plot(Rplot, Zplot, '-r', linewidth=2)
plt.show()

