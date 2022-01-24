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
from matplotlib import pyplot as plt


parnt_dir_name = os.path.dirname(os.getcwd())

fname_in_txt = '{0}/{1}/{2}'.format(parnt_dir_name, 'output_grid_files', 'grid.GS2_out_Miller_idx_42_nperiod_2_nl5_nt763')

arr1 = np.loadtxt('./eik5.out', skiprows = 1)
arr2 = np.loadtxt('./eik6.out', skiprows = 3)


f = open(fname_in_txt, 'r')
lines = f.readlines()

nlambda = eval(lines[1].split()[0])
lambda_arr = np.zeros((nlambda,))

for i in range(nlambda):
    lambda_arr[i] = [eval(val) for val in lines[3+i].split()][0] 

ntgrid, nperiod, ntheta, drhodpsi, rmaj, shat,  kxfac, q =  [eval(val) for val in lines[4+nlambda].split()]

# read from the python script
#5 + nlambda is useless becaue it contains variable names
gbdrift = np.zeros((ntheta+1, ))
gradpar = np.zeros((ntheta+1, ))
grho = np.zeros((ntheta+1, ))
tgrid = np.zeros((ntheta+1, ))

for i in range(ntheta+1):
    gbdrift[i], gradpar[i], grho[i], tgrid[i] = [eval(val) for val in lines[6+nlambda+i].split()]   

# read from eiktest output
gbdrift2 = arr1[:, 4]
gradpar2 = arr1[:, 1]
grho2 = arr1[:, 3]
tgrid2 = arr1[:, 0]

# read from the python script
cvdrift = np.zeros((ntheta+1, ))
gds2 = np.zeros((ntheta+1, ))
bmag = np.zeros((ntheta+1, )) 
tgrid = np.zeros((ntheta+1, ))

for i in range(ntheta+1):
    cvdrift[i], gds2[i], bmag[i], tgrid[i] = [eval(val) for val in lines[6+nlambda+ntheta+2+i].split()]   

# read from the eiktest output
cvdrift2 = arr1[:, 6]
gds2_2 = arr1[:, 8]
bmag2 = arr1[:, 2]
tgrid2 = arr1[:, 0]

# read from the python script
gds21 = np.zeros((ntheta+1, ))
gds22 = np.zeros((ntheta+1, ))
tgrid = np.zeros((ntheta+1, ))

for i in range(ntheta+1):
    gds21[i], gds22[i], tgrid[i] = [eval(val) for val in lines[6+nlambda+2*(ntheta+2)+i].split()]   

# read from the eiktest output
gds21_2 = arr1[:, 9]
gds22_2 = arr1[:, 10]
tgrid2 = arr1[:, 0]


# read from the python script
cvdrift0 = np.zeros((ntheta+1, ))
gbdrift0 = np.zeros((ntheta+1, ))
tgrid = np.zeros((ntheta+1, ))

for i in range(ntheta+1):
    cvdrift0[i], gbdrift0[i], tgrid[i] = [eval(val) for val in lines[6+nlambda+4*(ntheta+2)+i].split()]   


# read from the eiktest output
cvdrift0_2 = arr1[:, 7]
gbdrift0_2 = arr1[:, 5]
tgrid2 = tgrid2

#Rplot and Zplot are wrong in GS2 for all versions <= 8.0.6. THe new version 8.1-RC fixes this problem.
#Rplot = np.zeros((ntheta, ))
#Rprime = np.zeros((ntheta, ))
#tgrid = np.zeros((ntheta, ))
#
#for i in range(ntheta):
#    Rplot[i], Rprime[i], tgrid2[i] = [eval(val) for val in lines[6+nlambda+4*(ntheta+2)+i].split()]   
#
#
#Zplot = np.zeros((ntheta, ))
#Zprime = np.zeros((ntheta, ))
#tgrid = np.zeros((ntheta, ))
#
#for i in range(ntheta):
#    Zplot[i], Zprime[i], tgrid[i] = [eval(val) for val in lines[6+nlambda+5*(ntheta+2)+i].split()]   
#
#
aplot = np.zeros((ntheta+1, ))
aprime = np.zeros((ntheta+1, ))
tgrid = np.zeros((ntheta+1, ))

for i in range(ntheta+1):
    aplot[i], aprime[i], tgrid[i] = [eval(val) for val in lines[6+nlambda+6*(ntheta+2)+i].split()]   


aplot2 = arr2[:, 3]
aprime2 = arr2[:, 6]
tgrid2 = arr2[:, 0]

fig, ((plt1, plt2, plt3, plt4), (plt5, plt6, plt7, plt8)) = plt.subplots(2,4)

plt.suptitle('GS2-local_eikcoefs_gen_norm-comparison', fontsize=16)

pdb.set_trace()
plt1.plot(tgrid, cvdrift,'-r', tgrid2, cvdrift2, '-g')
plt1.set_xlabel('theta_equal_arc', fontsize=16)
plt1.set_ylabel('cvdrift', fontsize=16)
#plt1.text(0.1, 0.2, 'scale=%.3f'%(np.max(gbdrift)/np.max(gbdrift2)))
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt1.legend(['GS2', 'loc_eik'], fontsize=14)


plt2.plot(tgrid, bmag,'-r', tgrid2, bmag2, '-g')
plt2.set_xlabel('theta_equal_arc', fontsize=16)
plt2.set_ylabel('bmag', fontsize=16)
#plt2.text(0.9, 0.2, 'scale=%.3f'%(np.max(bmag)/np.max(bmag2)))
#plt2.set_title('rho=0.25', fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt2.legend(['GS2', 'local_eik'], fontsize=14)

plt5.plot(tgrid, gds22, '-r', tgrid2, gds22_2, '-g')
plt5.set_xlabel('theta_equal_arc', fontsize=16)
plt5.set_ylabel('gds22', fontsize=16)
#plt2.text(0.9, 0.2, 'scale=%.3f'%(np.max(bmag)/np.max(bmag2)))
#plt2.set_title('rho=0.25', fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt5.legend(['GS2', 'local_eik'], fontsize=14)


plt7.plot(tgrid, gradpar, '-r', tgrid2, gradpar2, '-g')
plt7.set_xlabel('theta_equal_arc', fontsize=16)
plt7.set_ylabel('gradpar2', fontsize=16)
#plt2.text(0.9, 0.2, 'scale=%.3f'%(np.max(bmag)/np.max(bmag2)))
#plt2.set_title('rho=0.25', fontsize=16)
#plt7.set_ylim([gradpar[0]-0.1, gradpar[0]+0.1])
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt7.legend(['GS2', 'local_eik'], fontsize=14)


plt3.plot(tgrid, gds2,'-r', tgrid2, gds2_2, '-g')
plt3.set_xlabel('theta_equal_arc', fontsize=16)
plt3.set_ylabel('gds2', fontsize=16)
#plt3.text(0.1, 0.1, 'scale=%.3f'%(np.min(gds2)/np.min(gds2_2)))
#plt3.set_title('rho=0.25', fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt3.legend(['GS2', 'local_eikc'], fontsize=14)


plt4.plot(tgrid, gds21,'-r', tgrid2, gds21_2, '-g')
plt4.set_xlabel('theta_equal_arc', fontsize=16)
plt4.set_ylabel('gds21', fontsize=16)
#plt4.text(0.9, 0.1, 'scale=%.3f'%(np.max(gds21)/np.max(gds21_2)))
#plt4.set_title('rho=0.25', fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt4.legend(['GS2', 'local_eik'], fontsize=14)


plt6.plot(tgrid, grho, '-r', tgrid2, grho2, '-g')
plt6.set_xlabel('theta_equal_arc', fontsize=16)
plt6.set_ylabel('grho', fontsize=16)
#plt2.text(0.9, 0.2, 'scale=%.3f'%(np.max(bmag)/np.max(bmag2)))
#plt2.set_title('rho=0.25', fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt6.legend(['GS2', 'local_eik'], fontsize=14)


plt8.plot(tgrid, aprime, '-r', tgrid2, aprime2, '-g')
plt8.set_xlabel('theta_equal_arc', fontsize=16)
plt8.set_ylabel('aprime', fontsize=16)
#plt2.text(0.9, 0.2, 'scale=%.3f'%(np.max(bmag)/np.max(bmag2)))
#plt2.set_title('rho=0.25', fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt8.legend(['GS2', 'local_eik'], fontsize=14)

plt.show()

pdb.set_trace()


