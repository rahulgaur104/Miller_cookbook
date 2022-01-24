#!/usr/bin/env python3
"""
This script saves the eikcoefs data coming from local_eikcoefs_gen_norm.py
"""

import os
import sys
import pdb
import numpy as np
import pickle
import netCDF4 as nc
from utils import *


eikcoefs_dict = sys.argv[1]
parnt_dir_nam = os.getcwd()

dict_file = open(eikcoefs_dict, 'rb')
eikcoefs_dict = pickle.load(dict_file)

theta       = eikcoefs_dict['theta_ex']
nperiod     = eikcoefs_dict['nperiod']
qfac        = eikcoefs_dict['qfac']
shat0       = eikcoefs_dict['shat']
dpsidrho    = eikcoefs_dict['dpsidrho']
gradpar     = eikcoefs_dict['gradpar_ex']
R           = eikcoefs_dict['R_ex']
Z           = eikcoefs_dict['Z_ex']
B           = eikcoefs_dict['B_ex']
gds21       = eikcoefs_dict['gds21_ex']
gds22       = eikcoefs_dict['gds22_ex']
gds2        = eikcoefs_dict['gds2_ex']
grho        = eikcoefs_dict['grho_ex']
gbdrift     = eikcoefs_dict['gbdrift_ex']
cvdrift     = eikcoefs_dict['cvdrift_ex']
gbdrift0    = eikcoefs_dict['gbdrift0_ex']
cvdrift0    = gbdrift0
jacob       = eikcoefs_dict['jacob']
aplot       = eikcoefs_dict['aplot']
aprime      = eikcoefs_dict['aprime']
fac         = eikcoefs_dict['fac']
file_idx    = eikcoefs_dict['file_idx']
lambda_knob = eikcoefs_dict['lambda_knob']
u_ML        = eikcoefs_dict['u_ML']

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

Rplot_ball = reflect_n_append(R, 'e')
Rprime_ball = reflect_n_append(nperiod_data_extend(np.sin(u_ML[theta <= np.pi]), nperiod, istheta=0, par='e'), 'e')

Zplot_ball = reflect_n_append(Z, 'o')
Zprime_ball = -reflect_n_append(nperiod_data_extend(np.cos(u_ML[theta <= np.pi]), nperiod, istheta=0, par='o'), 'o')

aplot_ball = reflect_n_append(aplot, 'o')
aprime_ball = reflect_n_append(aprime, 'o')

jacob_ball = reflect_n_append(jacob, 'e')


ntheta   = len(theta_ball)


##################################################################################################################
###########################---------------------GX NETCDF SAVE-------------------------------######################
##################################################################################################################

# creating equispaced equal-arc data
ntheta2       = ntheta - 1
theta_ball2   = np.delete(theta_ball, int(ntheta)-1)
gradpar_sav   = np.interp(theta_ball2, theta_ball,  gradpar_ball)
bmag_sav      = np.interp(theta_ball2, theta_ball,  B_ball)
jacob_sav     = np.interp(theta_ball2, theta_ball,  jacob_ball)
grho_sav      = np.interp(theta_ball2, theta_ball,  grho_ball)
gbdrift_sav   = np.interp(theta_ball2, theta_ball,  gbdrift_ball)
gbdrift0_sav  = np.interp(theta_ball2, theta_ball,  gbdrift0_ball)
cvdrift_sav   = np.interp(theta_ball2, theta_ball,  cvdrift_ball)
gds21_sav     = np.interp(theta_ball2, theta_ball,  gds21_ball)
gds2_sav      = np.interp(theta_ball2, theta_ball,  gds2_ball)
gds22_sav     = np.interp(theta_ball2, theta_ball,  gds22_ball)


Rplot_sav     = np.interp(theta_ball2, theta_ball,  Rplot_ball)
Zplot_sav     = np.interp(theta_ball2, theta_ball,  Zplot_ball)
aplot_sav     = np.interp(theta_ball2, theta_ball,  aplot_ball)
Rprime_sav    = np.interp(theta_ball2, theta_ball,  Rprime_ball)
Zprime_sav    = np.interp(theta_ball2, theta_ball,  Zprime_ball)
aprime_sav    = np.interp(theta_ball2, theta_ball,  aprime_ball)


fn = "./output_grid_files/gx_out_%d_nperiod_%d_nt%d.nc"%(int(file_idx), nperiod, ntheta2)

ds = nc.Dataset(fn, 'w')

z_nc = ds.createDimension('z', ntheta2)
#scalar = ds.createDimension('scalar', 1)


theta_nc    = ds.createVariable('theta', 'f8', ('z',))
bmag_nc     = ds.createVariable('bmag', 'f8', ('z',))
gradpar_nc  = ds.createVariable('gradpar', 'f8', ('z',))
grho_nc     = ds.createVariable('grho', 'f8', ('z',))
gds2_nc     = ds.createVariable('gds2', 'f8', ('z',))
gds21_nc    = ds.createVariable('gds21', 'f8', ('z',))
gds22_nc    = ds.createVariable('gds22', 'f8', ('z',))
gbdrift_nc  = ds.createVariable('gbdrift', 'f8', ('z',))
gbdrift0_nc = ds.createVariable('gbdrift0', 'f8', ('z',))
cvdrift_nc  = ds.createVariable('cvdrift', 'f8', ('z',))
cvdrift0_nc = ds.createVariable('cvdrift0', 'f8', ('z',))
jacob_nc    = ds.createVariable('jacob', 'f8', ('z',))


Rplot_nc    = ds.createVariable('Rplot', 'f8', ('z',))
Zplot_nc    = ds.createVariable('Zplot', 'f8', ('z',))
aplot_nc    = ds.createVariable('aplot', 'f8', ('z',))
Rprime_nc   = ds.createVariable('Rprime', 'f8', ('z',))
Zprime_nc   = ds.createVariable('Zprime', 'f8', ('z',))
aprime_nc   = ds.createVariable('aprime', 'f8', ('z',))


drhodpsi_nc = ds.createVariable('drhodpsi', 'f8', )
kxfac_nc    = ds.createVariable('kxfac', 'f8', )
Rmaj_nc     = ds.createVariable('Rmaj', 'f8', )
q           = ds.createVariable('q', 'f8', )
shat        = ds.createVariable('shat', 'f8', )


theta_nc[:]    = theta_ball2
bmag_nc[:]     = bmag_sav
gradpar_nc[:]  = gradpar_sav
grho_nc[:]     = grho_sav
gds2_nc[:]     = gds2_sav
gds21_nc[:]    = gds21_sav
gds22_nc[:]    = gds22_sav
gbdrift_nc[:]  = gbdrift_sav
gbdrift0_nc[:] = gbdrift0_sav
cvdrift_nc[:]  = cvdrift_sav
cvdrift0_nc[:] = gbdrift0_sav
jacob_nc[:]    = jacob_sav

Rplot_nc[:]    = Rplot_sav
Zplot_nc[:]    = Zplot_sav
aplot_nc[:]    = aplot_sav

Rprime_nc[:]   = Rprime_sav
Zprime_nc[:]   = Zprime_sav
aprime_nc[:]   = aprime_sav

drhodpsi_nc[0] = 1/dpsidrho
kxfac_nc[0]    = 1.
Rmaj_nc[0]     = (np.max(Rplot_nc) + np.min(Rplot_nc))/2
q[0]           = qfac
shat[0]        = shat0

ds.close()

print('GX file saved succesfully in the dir output_grid_files\n')

##################################################################################################################
###########################---------------------GS2 SAVE----------------------------------########################
##################################################################################################################

if np.max(theta) > np.pi:
    B0 = B_ball[theta_ball <= np.pi]
    theta_ball0 = theta_ball[theta_ball <= np.pi]
    B0 = B0[theta_ball0 >= 0]
    theta_ball0 = theta_ball0[theta_ball0 >= 0]
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
if lambda_knob == 1:
    from matplotlib import pyplot as plt
    plt.plot(np.array(theta_uni), np.array(B1) , '-sg', ms=3)
    #plt.plot(theta_ball0, B0 , '-sg', ms=3)
    plt.hlines(1/(lambda_arr) ,  xmin=-1, xmax=4)
    plt.show()

if response == 'n':
    sys.exit()
    print('returning back to the main script \n')


A1 = []
A2 = []
A3 = []
A4 = []
A5 = []
A6 = []
A7 = []
A8 = [] 

for i in range(ntheta):
	A2.append('    %.9f    %.9f    %.9f    %.9f\n'%(gbdrift_ball[i], gradpar_ball[i], grho_ball[i], theta_ball[i]))
	A3.append('    %.9f    %.9f    %.12f    %.9f\n'%(cvdrift_ball[i], gds2_ball[i], B_ball[i], theta_ball[i]))
	A4.append('    %.9f    %.9f    %.9f\n'%(gds21_ball[i], gds22_ball[i], theta_ball[i]))
	A5.append('    %.9f    %.9f    %.9f\n'%(gbdrift0_ball[i], gbdrift0_ball[i], theta_ball[i]))
	A6.append('    %.9f    %.9f    %.9f\n'%(Rplot_ball[i], Rprime_ball[i], theta_ball[i]))
	A7.append('    %.9f    %.9f    %.9f\n'%(Zplot_ball[i], Zprime_ball[i], theta_ball[i]))
	A8.append('    %.9f    %.9f    %.9f\n'%(aplot_ball[i], aprime_ball[i], theta_ball[i]))


A1.append([A2, A3, A4, A5, A6, A7, A8])
A1 = A1[0]


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
        g.write(A1[i-3][j])
g.close()

print('GS2 file saved successfully! returning back to the main script')



























