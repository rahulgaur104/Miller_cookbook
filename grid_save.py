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
from matplotlib import pyplot as plt

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
Z        = eikcoefs_dict['Z_ex']
B        = eikcoefs_dict['B_ex']
gds21    = eikcoefs_dict['gds21_ex']
gds22    = eikcoefs_dict['gds22_ex']
gds2     = eikcoefs_dict['gds2_ex']
grho     = eikcoefs_dict['grho_ex']
gbdrift  = eikcoefs_dict['gbdrift_ex']
cvdrift  = eikcoefs_dict['cvdrift_ex']
gbdrift0 = eikcoefs_dict['gbdrift0_ex']
cvdrift0 = gbdrift0
aplot    = eikcoefs_dict['aplot']
aprime   = eikcoefs_dict['aprime']
fac      = eikcoefs_dict['fac']
file_idx = eikcoefs_dict['file_idx']


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


##################################################################################################################
###########################---------------------GRYFX SAVE----------------------------------######################
##################################################################################################################

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

char = "grid.gryfx_out_Miller_%d_nperiod_%d_nt%d"%(int(file_idx), nperiod, ntheta)


fname_in_txt_rescaled = '{0}/{1}/{2}'.format(parnt_dir_nam, 'output_grid_files', char)
g = open(fname_in_txt_rescaled, 'w')


headings = ['ntgrid nperiod ntheta drhodpsi rmaj shat kxfac q\n', 'gbdrift gradpar grho tgrid\n', 'cvdrift gds2 bmag tgrid\n',\
            'gds21 gds22 tgrid\n', 'cvdrift0 gbdrift0 tgrid\n']

g.writelines(headings[0])
g.writelines('  %d    %d    %d   %0.3f   %0.1f    %.7f   %.1f   %.4f\n'%((ntheta-1)/2, 1, (ntheta-1), np.abs(1/dpsidrho), 1.0,\
            shat, 1.0, qfac))

for i in np.arange(1, len(headings)):
    g.writelines(headings[i])
    for j in range(ntheta):
            g.write(A1[i-1][j])
g.close()


print('gryfx file saved successfully!\n')

##################################################################################################################
###########################---------------------GS2 SAVE----------------------------------########################
##################################################################################################################

if np.max(theta) > np.pi:
    B0 = B_ball[theta_ball <= np.pi]
    theta_ball0 = theta_ball[theta_ball <= np.pi]
    B0 = B0[theta_ball0 >= 0]
    theta_ball0 = theta_ball0[theta_ball >= 0]
else:
    B0 = B
    theta_ball0 = theta

temp1 =  find_peaks(-reflect_n_append(B0, 'e'))[0]

assert len(temp1) <= 1, "more than 1 mag_wells(grid_save) Can't write GS2 grid file"

#thresh = np.mean(np.diff(theta_ball0))*(int(np.log(ntheta))*0.3)
thresh = np.mean(np.diff(B0))*(int(np.log(ntheta))*0.3+1)
thresh_sum = 0
theta_uni = [theta_ball0[0]]
B1 = [B0[0]]

for i in np.arange(1, len(theta_ball0)-1):
    #if (theta_ball0[i+1] - theta_ball0[i]) < thresh:
    if np.abs(B0[i+1] - B0[i]) < thresh and thresh_sum <= int(0.8*np.log(ntheta)+1)*thresh:
        thresh_sum = thresh_sum + thresh
    else:
        theta_uni.append(theta_ball0[i])
        B1.append(B0[i+1])
        thresh_sum = 0.

if np.max(np.array(B1)) != np.max(B0):
    theta_uni.append(theta_ball0[-1])
    B1.append(np.max(B0[-1]))


#if len(temp1) == 1 and temp1.item() == int((ntheta-3)/2):
#    lambda_arr = lambda_create(B_ball, fac)
#else:
lambda_arr = lambda_create(np.array(B1), fac)
#lambda_arr = lambda_create(B0, fac)

nlambda = len(lambda_arr)
print("The number of lambda points = %d\n"%int(nlambda))

response = input("Do you want to save a GS2 grid file with %d lambda points? (y/n)\n"%(nlambda))

#pdb.set_trace()
plt.plot(np.array(theta_uni), np.array(B1) , '-sg', ms=3)
#plt.plot(theta_ball0, B0 , '-sg', ms=3)
plt.hlines(1/(lambda_arr) ,  xmin=-10, xmax=10)
plt.show()

if response == 'n':
    sys.exit()
    print('returning back to the main script \n')


char = "grid.GS2_out_Miller_idx_%d_nperiod_%d_nl%d_nt%d"%(int(file_idx), nperiod, nlambda, ntheta)
fname_in_txt_rescaled = '{0}/{1}/{2}'.format(parnt_dir_nam, 'output_grid_files', char)

g = open(fname_in_txt_rescaled, 'w')
headings = ['nlambda\n', 'lambda\n', 'ntgrid nperiod ntheta drhodpsi rmaj shat kxfac q\n', 'gbdrift gradpar grho tgrid\n',\
            'cvdrift gds2 bmag tgrid\n', 'gds21 gds22 tgrid\n', 'cvdrift0 gbdrift0 tgrid\n', 'Rplot Rprime tgrid\n',\
            'Zplot Zprime tgrid\n', 'aplot aprime tgrid\n']
g.write(headings[0])
g.write('%d\n'%(nlambda))
g.writelines(headings[1])

for i in range(nlambda):
    g.writelines('%.19f\n'%(lambda_arr[i]))

g.writelines(headings[2])
g.writelines('  %d    %d    %d   %0.4f   %0.4f    %.6f   %.1f   %.5f\n'%((ntheta-1)/2, 1, (ntheta-1), np.abs(1/dpsidrho), 1.0,\
            shat, 1.0, qfac))

for i in np.arange(3, len(headings)):
    g.writelines(headings[i])
    for j in range(ntheta):
        if i <= 6:
            g.write(A1[i-3][j])
        else:
            g.write(A1[-1][j])
g.close()

print('GS2 file saved successfully! returning back to the main script')



























