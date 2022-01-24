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

fname_in_txt = '{0}/{1}/{2}'.format(parnt_dir_name, 'output_grid_files', 'gx_out_42_nperiod_2_nt763.nc')

fname_in_txt2 = '{0}/{1}/{2}'.format(parnt_dir_name, 'test', 'gs2_test_input.out.nc') 

rtg = Dataset(fname_in_txt, 'r')

#reading output of the local_eikcoefs_gen script
tgrid     = rtg.variables['theta'][:].data
gbdrift   = rtg.variables['gbdrift'][:].data
gradpar   = rtg.variables['gradpar'][:].data
grho      = rtg.variables['grho'][:].data
cvdrift   = rtg.variables['cvdrift'][:].data
gds2      = rtg.variables['gds2'][:].data
bmag      = rtg.variables['bmag'][:].data
gds21     = rtg.variables['gds21'][:].data
gds22     = rtg.variables['gds22'][:].data
cvdrift0  = rtg.variables['cvdrift0'][:].data
gbdrift0  = rtg.variables['gbdrift0'][:].data
aplot     = rtg.variables['aplot'][:].data
aprime    = rtg.variables['aprime'][:].data
Rplot     = rtg.variables['Rplot'][:].data
Rprime    = rtg.variables['Rprime'][:].data
Zplot     = rtg.variables['Zplot'][:].data
Zprime    = rtg.variables['Zprime'][:].data
jacob     = rtg.variables['jacob'][:].data

rtg2 = Dataset(fname_in_txt2, 'r')

# read from GS2 output
tgrid2     = rtg2.variables['theta'][:].data
gbdrift2   = rtg2.variables['gbdrift'][:].data
gradpar2   = rtg2.variables['gradpar'][:].data
grho2      = rtg2.variables['grho'][:].data
cvdrift2   = rtg2.variables['cvdrift'][:].data
gds2_2     = rtg2.variables['gds2'][:].data
bmag2      = rtg2.variables['bmag'][:].data
gds21_2    = rtg2.variables['gds21'][:].data
gds22_2    = rtg2.variables['gds22'][:].data
cvdrift0_2 = rtg2.variables['cvdrift0'][:].data
gbdrift0_2 = rtg2.variables['gbdrift0'][:].data
aplot2     = rtg2.variables['aplot'][:].data
aprime2    = rtg2.variables['aprime'][:].data
Rplot2     = rtg2.variables['Rplot'][:].data
Rprime2    = rtg2.variables['Rprime'][:].data
Zplot2     = rtg2.variables['Zplot'][:].data
Zprime2    = rtg2.variables['Zprime'][:].data
jacob2     = rtg2.variables['jacob'][:].data



fig, ((plt1, plt2, plt3, plt4), (plt5, plt6, plt7, plt8)) = plt.subplots(2,4)

plt.suptitle('GS2-local_eikcoefs_gen_norm-comparison', fontsize=16)

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

plt5.plot(tgrid, gds2, '-r', tgrid2, gds2_2, '-g')
plt5.set_xlabel('theta_equal_arc', fontsize=16)
plt5.set_ylabel('gds2', fontsize=16)
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


#plt3.plot(tgrid, jacob,'-r', tgrid2, jacob2, '-g')
plt3.plot(tgrid, gds22,'-r', tgrid2, gds22_2, '-g')
plt3.set_xlabel('theta_equal_arc', fontsize=16)
#plt3.set_ylabel('jacob', fontsize=16)
plt3.set_ylabel('gds22', fontsize=16)
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


plt6.plot(tgrid, gbdrift, '-r', tgrid2, gbdrift2, '-g')
plt6.set_xlabel('theta_equal_arc', fontsize=16)
#plt6.set_ylabel('grho', fontsize=16)
plt6.set_ylabel('gbdrift', fontsize=16)
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

fig.set_size_inches(16, 8, forward=True)
fig.tight_layout()
plt.savefig('GS2-local-eikcoefs-comparison.png', dpi=300)
print("png saved!")

#plt.show()
#pdb.set_trace()


