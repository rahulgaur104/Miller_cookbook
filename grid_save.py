#!/usr/bin/env python3
"""
This script saves the eikcoefs data coming from local_eikcoefs_gen.py
"""

import os
import sys
import pdb
import numpy as np
import pickle
from utils import *

eikcoefs_dict = sys.argv[1]
parnt_dir_nam = os.getcwd()

dict_file = open(eikcoefs_dict, 'rb')
eikcoefs_dict = pickle.load(dict_file)


theta    = eikcoefs_dict['theta_ex']
nperiod  = eikcoefs_dict['nperiod']
qfac     = eikcoefs_dict['qfac']
shat     = eikcoefs_dict['shat']
dpsidrho = eikcoefs_dict['dpsidrho']
gradpar  = eikcoefs_dict['gradpar_ex']
R        = eikcoefs_dict['R_ex']
B        = eikcoefs_dict['B_ex']
gds21    = eikcoefs_dict['gds21_ex']
gds22    = eikcoefs_dict['gds22_ex']
gds2     = eikcoefs_dict['gds2_ex']
grho     = eikcoefs_dict['grho_ex']
gbdrift  = eikcoefs_dict['gbdrift_ex']
cvdrift  = eikcoefs_dict['cvdrift_ex']
gbdrift0 = eikcoefs_dict['gbdrift0_ex']
cvdrift0 = gbdrift0


gradpar_ball = reflect_n_append(gradpar, 'e')
theta_ball = reflect_n_append(theta, 'o')
cvdrift_ball = reflect_n_append(cvdrift, 'e')
gbdrift_ball = reflect_n_append(gbdrift, 'e')
gbdrift0_ball = reflect_n_append(gbdrift0, 'o')
B_ball = reflect_n_append(B, 'e')
gds2_ball = reflect_n_append(gds2, 'e')
gds21_ball = reflect_n_append(gds21, 'o')
gds22_ball = reflect_n_append(gds22 , 'e')
grho_ball = reflect_n_append(grho , 'e')

ntheta   = len(theta_ball)


pdb.set_trace()

A1 = []
A2 = []
A3 = [] 
A4 = [] 
A5 = [] 

for i in range(ntheta):
	A2.append('    %.9f    %.9f    %.9f    %.9f\n'%(gbdrift_ball[i], gradpar_ball[i], grho_ball[i], theta_ball[i]))
	A3.append('    %.9f    %.9f    %.12f    %.9f\n'%(cvdrift_ball[i], gds2_ball[i], B_ball[i], theta_ball[i]))
	A4.append('    %.9f    %.9f    %.9f\n'%(gds21_ball[i], gds22_ball[i], theta_ball[i]))
	A5.append('    %.9f    %.9f    %.9f\n'%(gbdrift0_ball[i], gbdrift0_ball[i], theta_ball[i]))


A1.append([A2, A3, A4, A5])
A1 = A1[0]

char = "grid.out_Miller_nperiod_%d_nt%d"%(nperiod, ntheta)


fname_in_txt_rescaled = '{0}/{1}/{2}'.format(parnt_dir_nam, 'output_grid_files', char)
g = open(fname_in_txt_rescaled, 'w')


headings = ['ntgrid nperiod ntheta drhodpsi rmaj shat kxfac q\n', 'gbdrift gradpar grho tgrid\n', 'cvdrift gds2 bmag tgrid\n', 'gds21 gds22 tgrid\n', 'cvdrift0 gbdrift0 tgrid\n']

g.writelines(headings[0])
g.writelines('  %d    %d    %d   %0.3f   %0.1f    %.7f   %.1f   %.4f\n'%((ntheta-1)/2, 1, (ntheta-1), np.abs(1/dpsidrho), 1.0, shat, 1.0, qfac))

for i in np.arange(1, len(headings)):
    g.writelines(headings[i])
    for j in range(ntheta):
            g.write(A1[i-1][j])
g.close()

print('file successfully save! returning back to the main script')


























