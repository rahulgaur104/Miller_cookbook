"""
The purpose of this script is to generate a local Miller equilibrium and compare various parameters of interest with 
eiktest(old routine on GS2) for the same equilibrium. In some ways, this script is the pythonized version of eiktest.

Additionally, it also performs a ballooning stability analysis based on Newcomb's theorem both in the collocation and 
the flux(straight-field-line) theta grid.

The derivatives are calculated usign a central-difference method. The integrals are performed using a trapezoidal sum.

Some of the figures of merit(FOMs) we calculate here are:
    total and polidal magnetic field
    geometric components of the grad-B and curvature drifts
    curvature of the field line
    local magnetic shear
    angle of the flux tube cross-section(the cross-section is a parallelogram so one angle defines the rest)

These FOMs can be plotted w.r.t different types of thetas:
    geometric theta
    collocation theta
    straight-field-line theta
    equal-arc theta

Figures of merit to be added in the future
    bounce-averaged curvature drift(precession drift)
"""
import os   
import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import cumtrapz as ctrap
from scipy.signal import find_peaks
from scipy.interpolate import InterpolatedUnivariateSpline as linspl
from scipy.interpolate import CubicSpline as cubspl
from scipy.integrate import cumtrapz as ctrap
import pdb 
import time
import multiprocessing as mp

parnt_dir_nam = os.path.dirname(os.getcwd())
# definitions of the following variables can be obtained from 
# https://gyrokinetics.gitlab.io/gs2/page/namelists/#theta_grid_eik_knobs

# Further information can be found from Miller's paper. Refernce provided
# in the repo
ntheta = 512

rhoc = 0.97
qinp = 8.3
shift = -0.66
s_hat_input = 8.3
R_geo = 1.47
Rmaj = 1.47
akappa = 2.2
akappri = 0.4


#rhoc = 0.98
#qinp = 5
#shift = -0.2
#s_hat_input = 5
#R_geo = 1.9
#Rmaj = 1.5
#akappa = 1.7
#akappri = 1
tri_gs2 = 0.5
tripri_gs2 = 0

delrho = 0.005


# Note that tri in gs2 is actually the sin(delta).
tri = np.sin(tri_gs2) # tri is tri Miller
tripri = (np.sin(tri_gs2+tripri_gs2*delrho) - np.sin(tri_gs2-tripri_gs2*delrho))/(2*delrho)


#beta_prime_input =  -0.24
beta_prime_input =  -0.2

no_of_surfs = 3
# note that this theta is neither geometric nor flux. It's just used to generate the surfaces. We will call it colocation theta
# GS2 uses this theta for the grad par calculation
theta = np.linspace(0, np.pi, ntheta)


# primary lowest level calculations
R_0 = np.array([Rmaj+np.abs(shift)*delrho, Rmaj, Rmaj-np.abs(shift)*delrho])
#drhodpsi = -1.38205 # always check with the value from eik6.out Probably related to B_r*a**2(a = rhoc/2), changes with shift. Why?

rho = np.array([rhoc - delrho, rhoc, rhoc + delrho])

qfac = np.array([qinp-s_hat_input*(qinp/rhoc)*delrho, qinp, qinp+s_hat_input*(qinp/rhoc)*delrho]) 
kappa = np.array([akappa-akappri*delrho, akappa, akappa+akappri*delrho])
delta = np.array([tri-tripri*delrho, tri, tri+tripri*delrho])

#R_mag_ax can be anything as long as it's inside the annulus. 
R_mag_ax = Rmaj

dpdrho = beta_prime_input/2  #This definiton with a factor of 2 has been taken directly from geometry.f90.


R= np.array([R_0[i] + (rho[i])*np.cos(theta +np.arcsin(delta[i])*np.sin(theta)) for i in range(no_of_surfs)])
Z = np.array([kappa[i]*(rho[i])*np.sin(theta) for i in range(no_of_surfs)])

# theta array with a common magnetic axis
theta_comn_mag_ax = np.array([np.arctan2(Z[i], R[i]-R_mag_ax) for i in range(no_of_surfs)])

for i in range(no_of_surfs):
	R[i] = np.interp(theta_comn_mag_ax[1], theta_comn_mag_ax[i], R[i])
	Z[i] = np.interp(theta_comn_mag_ax[1], theta_comn_mag_ax[i], Z[i])
	theta_comn_mag_ax[i] = theta_comn_mag_ax[1]



def derm(arr, ch, par='e'):
    # Finite difference subroutine
    # ch = 'l' means difference along the flux surface
    # ch = 'r' mean difference across the flux surfaces
    # par corresponds to parity of the equilibrium quantities, i.e., up-down symmetry or anti-symmetry
    # par = 'e' means even parity of the arr. PARITY OF THE INPUT ARRAY 
    # par = 'o' means odd parity
    # par is only useful for ch = 'l'
    # arr must be in the range [-pi, pi]
    # This routine is only valid for up-down symmetric Miller equilibria 

    temp = np.shape(arr)
    if len(temp) == 1 and ch == 'l': #finite diff along the flux surface for a single array
    	#pdb.set_trace()
    	if par == 'e':
    		d1, d2 = np.shape(arr)[0], 1
    		arr = np.reshape(arr, (d2,d1))
    		diff_arr = np.zeros((d2,d1))
    		diff_arr[0, 0] =  0. #(arr_theta_0- - arr_theta_0+)  = 0
    		diff_arr[0, -1] = 0. #(arr_theta_pi- - arr_theta_pi+)  = 0
    		diff_arr[0, 1:-1] = np.diff(arr[0,:-1], axis=0) + np.diff(arr[0,1:], axis=0)  
    	else:
    		d1, d2 = np.shape(arr)[0], 1
    		arr = np.reshape(arr, (d2,d1))
    		diff_arr = np.zeros((d2,d1))
    		diff_arr[0, 0] = 2*(arr[0, 1] - arr[0, 0]) 
    		diff_arr[0, -1] = 2*(arr[0, -1] - arr[0, -2]) 
    		diff_arr[0, 1:-1] = np.diff(arr[0,:-1], axis=0) + np.diff(arr[0,1:], axis=0)  

    elif len(temp) == 1 and ch == 'r': # across surfaces for a single array
    	#pdb.set_trace()
    	d1, d2 = np.shape(arr)[0], 1
    	diff_arr = np.zeros((d1,d2))
    	arr = np.reshape(arr, (d1,d2))
    	diff_arr[0, 0] = 2*(arr[1, 0] - arr[0, 0]) # single dimension arrays like psi, F and q don't have parity
    	diff_arr[-1, 0] = 2*(arr[-1, 0] - arr[-2, 0]) 
    	diff_arr[1:-1, 0] = np.diff(arr[:-1,0], axis=0) + np.diff(arr[1:,0], axis=0)  


    else:
    	d1, d2 = np.shape(arr)[0], np.shape(arr)[1]

    	diff_arr = np.zeros((d1,d2))
    	if ch == 'r': # across surfaces for multi-dim array
    		#pdb.set_trace()
    		diff_arr[0, :] = 2*(arr[1,:] - arr[0,:]) 
    		diff_arr[-1, :] = 2*(arr[-1,:] - arr[-2,:]) 
    		diff_arr[1:-1, :] = (np.diff(arr[:-1,:], axis=0) + np.diff(arr[1:,:], axis=0))  

    	else: #along a surface for a multi-dim array
    		#pdb.set_trace()
    		if par == 'e':
    			#pdb.set_trace()
    			diff_arr[:, 0] = np.zeros((d1,))
    			diff_arr[:, -1] = np.zeros((d1,))
    			diff_arr[:, 1:-1] = (np.diff(arr[:,:-1], axis=1) + np.diff(arr[:,1:], axis=1))  
    		else:
    			diff_arr[:, 0] = 2*(arr[:, 1] - arr[:, 0]) 
    			diff_arr[:, -1] = 2*(arr[:, -1] - arr[:, -2]) 
    			diff_arr[:, 1:-1] = (np.diff(arr[:,:-1], axis=1) + np.diff(arr[:,1:], axis=1))  

    arr = np.reshape(diff_arr, temp)
    return diff_arr


def intersection_chk(R, Z, R_mag_ax):
    # Inputs: R and Z arrays on a common theta grid
    # Output: number informing whether the surfaces intersect
    r = np.sqrt((R-R_mag_ax)**2 + Z**2)
    diffr = derm(r, ch = 'r')
    if np.min(diffr) <= 0:
        return 1
    else:
        return 0


# to check is the equilibrium sufaces intersect with each other
if intersection_chk(R, Z, R_mag_ax) != 0:
	print("WARNING! Surfaces intersect...")
else:
	print("curve intersection check passed... surfaces do not intersect")


def nperiod_data_extend(arr, nperiod, istheta=0, par='e'):
    # the purpose of this routine is to extend the values of various parameters from their usual range of [0, pi]
    # to [0, (2*nperiod-1)*pi]
    # par = 'e' is the parity of arr. It can either be up-down symmetric(par = 'e') or up-down asymmetric(par='o')
    if nperiod > 1:
        if istheta: #for istheta par='o'
            arr_dum = arr
            for i in range(nperiod-1):
                arr_app = np.concatenate((2*np.pi*(i+1)-arr_dum[::-1][1:], 2*np.pi*(i+1)+arr_dum[1:]))
                arr = np.concatenate((arr, arr_app))
        else:
            if par == 'e':
                arr_app = np.concatenate((arr[::-1][1:], arr[1:]))
                for i in range(nperiod-1):
                    arr = np.concatenate((arr, arr_app))
            else:
                arr_app = np.concatenate((-arr[::-1][1:], arr[1:]))
                for i in range(nperiod-1):
                    arr = np.concatenate((arr, arr_app))
    return arr       

def reflect_n_append(arr, ch):
	"""
	The purpose of this function is to increase the span of an array from [0, (2*nperiod-1)*np.pi] to [-(2*nperiod-1)*pi, (2*nperiod-1)*pi].
        ch can either be 'e'(even) or 'o'(odd) depending upon the parity of the input array.
	"""
	rows = 1
	brr = np.zeros((2*len(arr)-1, ))
	if ch == 'e':
		for i in range(rows):
                    brr = np.concatenate((arr[::-1][:-1], arr[0:]))
	else :
		for i in range(rows):
                    brr = np.concatenate((-arr[::-1][:-1],np.array([0.]), arr[1:]))
	return brr

def find_optim_theta_arr(arr, theta_arr, res_par):
        # The purpose of this routine is to optimize the size of the theta array so that one can keep all the important features 
        # with the minimum number of theta points. This routine is only used when one wants to save a grid.out file for a GS2 run
        rows, colms = np.shape(arr)
        idx = []
        idx2 = []
        idx3 = []
        idx4 = []
        for i in range(rows):
            peaks, _ = find_peaks(arr[i], height=-1E10)
            peaks = peaks.astype(np.int)
            idx.append(np.ndarray.tolist(peaks))
            peaks2, _ = find_peaks(-arr[i], height=-1E10)
            idx.append(np.ndarray.tolist(peaks2))

        idx.append([0, len(theta_arr)-1])    
        idx = np.sum(idx)
        idx = list(set(idx))
        idx.sort()
        comb_peaks = np.array(idx)
        diff_peaks = np.sort(np.unique(np.diff(np.sort(comb_peaks))))
        diff_peaks = diff_peaks[diff_peaks>8]
        #pdb.set_trace()
        diff = int(diff_peaks[0]/2)
        comb_peaks = np.sort(np.abs(np.concatenate((peaks-diff, peaks, peaks+diff, peaks2-diff, peaks2, peaks2+diff, np.array([0, len(theta_arr)-1-diff])))))
        diff2 = int(np.mean(np.diff(comb_peaks)))-2

        comb_peaks_diff = np.diff(comb_peaks)
        idx_gt_diff2 = np.where(comb_peaks_diff>diff2)[0][:]
        for i in idx_gt_diff2:
            j = comb_peaks[i]
            #pdb.set_trace()
            while j < comb_peaks[i+1]:
                idx2.append(j+diff2)
                j = j + diff2
        comb_peaks = np.concatenate((comb_peaks, np.array(idx2)))
        comb_peaks = np.concatenate((comb_peaks, np.array([len(theta_arr)-1])))
        comb_peaks = comb_peaks[comb_peaks < len(theta_arr)]
        comb_peaks = np.sort(np.unique(comb_peaks))
        #pdb.set_trace()

        return theta_arr[comb_peaks]

def lambda_create(arr):
    arr1 = np.sort(np.unique(1/arr))
    diff_arr1 = np.diff(arr1)
    req_diff = np.mean(diff_arr1)/10
    idx = [arr1[0], arr1[-1]]
    
    diff_arr_sum = 0
    i = 1
    while i < len(arr1)-1:
        if diff_arr1[i] <= req_diff:
            diff_arr_sum = diff_arr_sum + req_diff
        else:
            idx.append(arr1[i])
            diff_arr_sum = 0
        i = i + 1

    return np.unique(np.array(idx))
        


dRj = np.zeros((no_of_surfs, ntheta))
dZj = np.zeros((no_of_surfs, ntheta))
dRi = np.zeros((no_of_surfs, ntheta))
dZi = np.zeros((no_of_surfs, ntheta))
dt = np.zeros((no_of_surfs, ntheta))
dti = np.zeros((no_of_surfs, ntheta))
dt_st_i = np.zeros((no_of_surfs, ntheta))
dBr_ML = np.zeros((no_of_surfs, ntheta))
dr = np.zeros((no_of_surfs, ntheta))
theta_u = np.zeros((no_of_surfs, ntheta))
theta_d = np.zeros((no_of_surfs, ntheta))
theta_st = np.zeros((no_of_surfs, ntheta))
phi_n = np.zeros((no_of_surfs, ntheta))
u_ML = np.zeros((no_of_surfs, ntheta))



###################-----------------GRADIENTS ON GEOMETRIC THETA GRID------------------#################

dl = np.sqrt(derm(R,'l','e')**2 + derm(Z,'l','o')**2)
dt = derm(theta_comn_mag_ax, 'l', 'o')


rho_diff = derm(rho, 'r')
# partial derivatives of R and Z on the exact rho and theta_geometric grid
dR_drho = derm(R, 'r')/rho_diff
dR_dt = derm(R, 'l', 'e')/dt
dZ_drho = derm(Z, 'r')/rho_diff
dZ_dt = derm(Z, 'l', 'o')/dt

jac = dR_drho*dZ_dt - dZ_drho*dR_dt

# partial derivatives of psi and theta_geometric on the cartesian grid
drhodR = dZ_dt/jac
drhodZ = -dR_dt/jac
dt_dR = -dZ_drho/jac
dt_dZ =  dR_drho/jac


dpsidrho_arr = -R_geo/np.abs(2*np.pi*qfac/(2*ctrap(jac/R, theta_comn_mag_ax)[:, -1])) # determining dpsidrho from the safety factor relation
dpsidrho = dpsidrho_arr[1]
#dpsidrho = -0.515198
F = np.ones((3,))*R_geo
drhodpsi = 1/dpsidrho

dpdpsi = dpdrho*drhodpsi

#pdb.set_trace()
psi = np.array([1-delrho/drhodpsi, 1, 1+delrho/drhodpsi])
psi_diff = derm(psi, 'r')

# partial theta/partial rho (radial component of grad theta)
dtdr_geo = np.sign(psi_diff)*(dt_dR*drhodR + dt_dZ*drhodZ)/np.sqrt(drhodR**2 + drhodZ**2)

B_p = np.abs(dpsidrho)*np.array([np.sqrt(drhodR[i]**2 + drhodZ[i]**2)/R[i] for i in range(no_of_surfs)])
B_t = np.array([np.reshape(F, (-1,1))[i]/R[i] for i in range(no_of_surfs)])
B2 = np.array([B_p[i]**2 + B_t[i]**2 for i in range(no_of_surfs)])
B = np.sqrt(B2)

# grad psi from the cartesian grid 
grad_psi_cart = dpsidrho*np.sqrt(drhodR**2 + drhodZ**2)

# gradpar_0 is b.grad(theta) where theta = collocation theta
gradpar_0 = 1/(R*B)*np.array([np.abs(dpsidrho_arr[i])*np.sqrt(drhodR[i]**2 + drhodZ[i]**2) for i in range(no_of_surfs)])*(derm(theta, 'l', 'o')/dl) 
# To reiterate, this theta is neither the geometric nor flux theta
# This calculation of gradpar_0 is only meaningful on the central surface as theta = collocation theta is only known as a function of 
# geometric theta on the central surface. On the adjacent surfaces, we don't know the relation b/w geometric and collocation theta.


################------------------GRADIENTS ON FLUX THETA GRID------------------------##################

# Calculating theta_f or theta_st from the cartesian derivatives
# Note that this theta_st is only meaningful for the central surface. This happens because we only know the exact
# value of F on the central surface. We cannot get the values of F on the adjacent surfaces without using the 
# Bishop recipe which comes later in this script. Please refer to the second reference in the repo README.md

for i in range(no_of_surfs):
	 theta_st[i, 1:] = ctrap(np.abs(np.reshape(F,(-1,1))[i]*(1/dpsidrho_arr[i])*jac[i]/R[i]), theta_comn_mag_ax[i])
	 #theta_st[i, 1:] = ctrap(np.abs(np.reshape(F,(-1,1))[i]*(1/dpsidrho)*jac[i]/R[i]), theta_comn_mag_ax[i])
	 theta_st[i, 1:] = theta_st[i, 1:]/theta_st[i, -1]
	 theta_st[i, 1:] = np.pi*theta_st[i, 1:]

#spline object b/w flux theta and collocation theta
spl1 = linspl(theta_st[1], theta)
#spline object b/w geometric theta and flux theta
th_geo_st_spl = linspl(theta_comn_mag_ax[1], theta_st[1], k = 1)

#Before we take gradients on the theta_st grid we need to interpolate R, Z and B on to a uniform theta_st grid
theta_st_new = np.linspace(0, np.pi, ntheta)*np.reshape(np.ones((no_of_surfs,)),(-1,1))
#theta_st_new = theta_st
theta_comn_mag_ax_new = np.zeros((no_of_surfs, ntheta))
B1 = np.zeros((1, ntheta))
B1 = B[1].copy()
# gradpar1 is b.grad(theta) where theta = flux theta
gradpar1 = 1/(B1)*(B_p[1])*(derm(theta_st[1], 'l', 'o')/dl[1])

for i in range(no_of_surfs):
        R[i] = np.interp(theta_st_new[i], theta_st[i], R[i])
        Z[i] = np.interp(theta_st_new[i], theta_st[i], Z[i])
        B[i] = np.interp(theta_st_new[i], theta_st[i], B[i])
        B_p[i] = np.interp(theta_st_new[i], theta_st[i], B_p[i])
        gradpar_0[i] = np.interp(theta_st_new[i], theta_st[i], gradpar_0[i])
        dtdr_geo[i] = np.interp(theta_st_new[i], theta_st[i], dtdr_geo[i])
        #dtdr_st[i] = np.interp(theta_st_new[i], theta_st[i], dtdr_st[i])
        theta_comn_mag_ax_new[i] = np.arctan2(Z[i], R[i]-R_mag_ax)


# partial derivatives of R and Z on the exact psi and theta_f grid
dt_st_l = derm(theta_st_new, 'l', 'o')
dR_dpsi = derm(R, 'r')/psi_diff
dR_dt = derm(R, 'l', 'e')/dt_st_l
dZ_dpsi = derm(Z, 'r')/psi_diff
dZ_dt = derm(Z, 'l', 'o')/dt_st_l

jac = dR_dpsi*dZ_dt - dZ_dpsi*dR_dt

# partial derivatives of psi and theta_f on the cartesian grid
dpsidR = dZ_dt/jac
dpsidZ = -dR_dt/jac
dt_dR = -dZ_dpsi/jac
dt_dZ =  dR_dpsi/jac


dtdr_st0 = (dt_dR*dpsidR + dt_dZ*dpsidZ)/np.sqrt(dpsidR**2 + dpsidZ**2)

# Recalculate dl on the new grid
dl = np.sqrt(derm(R,'l', 'e')**2 + derm(Z,'l', 'o')**2)

dt = derm(theta_comn_mag_ax_new, 'l', 'o')
for i in range(no_of_surfs):
	 dRj[i, :] = derm(R[i,:], 'l', 'e')
	 dZj[i, :] = derm(Z[i,:], 'l', 'o') 
	 phi = np.arctan2(dZj[i,:], dRj[i,:])
	 phi = np.concatenate((phi[phi>=0]-np.pi/2, phi[phi<0]+3*np.pi/2)) 
	 phi_n[i,:] = phi

u_ML = np.arctan2(derm(Z, 'l', 'o'), derm(R, 'l', 'e'))

# du_ML/dl is negative and dphi = -du_ML so R_c = -du_ML/dl > 0
#R_c = dl/(2*np.concatenate((np.diff(phi_n, axis=1), np.reshape(np.diff(phi_n)[:, -1],(-1,1))), axis=1))
R_c = dl/derm(phi_n, 'l', 'o')

#pdb.set_trace()
gradpar2 = 1/(B[1])*(B_p[1])*(derm(theta_st_new[1], 'l', 'o')/dl[1]) # gradpar is b.grad(theta)


nperiod = 1
B_p_ex = nperiod_data_extend(np.abs(B_p[1]), nperiod)
B_ex = nperiod_data_extend(B[1], nperiod)
R_ex = nperiod_data_extend(R[1], nperiod)
theta_col = spl1(theta_st_new[1])
theta_col_ex = nperiod_data_extend(theta_col, nperiod, istheta=1)
theta_st_new_ex = nperiod_data_extend(theta_st_new[1], nperiod, istheta=1)
u_ML_ex = nperiod_data_extend(u_ML[1], nperiod)
R_c_ex = nperiod_data_extend(R_c[1], nperiod)
dl_ex = nperiod_data_extend(dl[1], nperiod)


diffrho = derm(rho, 'r')


# Since we are calculating these coefficients in straight field line theta, we can use the fact that F[1]*jac[1]/R[1] = qfac[1]

a_s = -(2*qfac[1]/F[1]*theta_st_new_ex + 2*F[1]*qfac[1]*ctrap(1/(R_ex**2*B_p_ex**2), theta_st_new_ex, initial=0))  
b_s = -(2*qfac[1]*ctrap(1/(B_p_ex**2), theta_st_new_ex, initial=0))
c_s =  (2*qfac[1]*ctrap((2*np.sin(u_ML_ex)/R_ex - 2/R_c_ex)*1/(R_ex*B_p_ex), theta_st_new_ex, initial=0))


# calculating the exact dFdpsi on the surface from relation 21 in Miller's paper.
dFdpsi = (-s_hat_input/(rho[1]*(psi_diff[1]/diffrho[1])*(1/(2*np.pi*qfac[1]*(2*nperiod-1))))-(b_s[-1]*dpdpsi - c_s[-1]))/a_s[-1]

#pdb.set_trace()
# psi_diff[1]/2 is essential
F[0], F[1], F[2]= F[1]-dFdpsi*(psi_diff[1]/2), F[1], F[1]+dFdpsi*(psi_diff[1]/2)

# Calculating the current from the relation (21) in Miller's paper(involving shat) and comparing it with F = q*R^2/J, where J = R*jac is the flux theta jacobian 
F_chk = np.array([np.abs(np.mean(qfac[i]*R[i]/jac[i])) for i in range(no_of_surfs)])



if np.max(np.abs((dt_dR[1]*dpsidZ[1] - dt_dZ[1]*dpsidR[1])/np.sqrt(dpsidR[1]**2 + dpsidZ[1]**2) - dt_st_l[1]/dl[1])) > 1E-11:
	print("grad theta_st along l doesn't match...")
else:
	print("grad theta_st along the surface test passed...")


if  np.abs(np.max((-dt_dR[1]*dpsidZ[1] + dpsidR[1]*dt_dZ[1])*jac[1]) - 1.0) > 1E-11:
	print("theta dot grad theta = 1 test failed...")
else:
	print("theta dot grad theta = 1 test passed...")


dpsi_dr = np.zeros((no_of_surfs, ntheta))
dpsi_dr = np.sign(psi_diff)*np.sqrt(dpsidR**2 + dpsidZ**2)
B_p1 = np.array([np.sqrt(dpsidR[i]**2 + dpsidZ[i]**2)/R[i] for i in range(no_of_surfs)])
B_p1_ex = nperiod_data_extend(B_p1[1], nperiod)

B_p = np.abs(dpsi_dr)/R
B_t = np.array([np.reshape(F, (-1,1))[i]/R[i] for i in range(no_of_surfs)])

B2 = B_p**2 + B_t**2
B = np.sqrt(B2)

B_p_ex = nperiod_data_extend(B_p[1], nperiod)
B_ex = nperiod_data_extend(B[1], nperiod)
B2_ex = nperiod_data_extend(B2[1], nperiod)
#pdb.set_trace()

#B2 = B_ex**2

#diffP = derm(pres, 'r')
dB2l = derm(B2, 'l')
dBl = derm(B, 'l')
diffq = derm(qfac, 'r')

#dB2l_ex = nperiod_data_extend(dB2l[1], nperiod, istheta=0, par='o')
#dBl_ex = nperiod_data_extend(dBl[1], nperiod, istheta=0, par='o')

#B_ex_1 = B_ex*np.ones((3,len(theta_st_new_ex)))
dB2l_ex = derm(B_ex**2, 'l')[0] # not dB[1]2l zero because the higher dimensional array
dBl_ex =  derm(B_ex, 'l')[0]



dpsi_dr_ex = nperiod_data_extend(dpsi_dr[1], nperiod)
gds22 = (diffq/diffrho)**2*np.abs(dpsi_dr_ex)**2  
alpha = -np.reshape(qfac,(-1,1))*theta_st_new_ex
grho = np.reshape(rho,(-1,1))*(diffrho)/(rhoc*psi_diff/dpsi_dr_ex)

dFdpsi =  derm(F, 'r')/psi_diff


dqdr = diffq*dpsi_dr_ex/psi_diff
dpdr = dpdpsi*dpsi_dr_ex

dpsidR_ex = nperiod_data_extend(dpsidR[1], nperiod)
dt_dR_ex = nperiod_data_extend(dt_dR[1], nperiod, par = 'o')
dt_dZ_ex = nperiod_data_extend(dt_dZ[1], nperiod)
dpsidZ_ex = nperiod_data_extend(dpsidZ[1], nperiod, istheta=0, par = 'o')

#dtdr_st_ex = (dt_dR_ex*dpsidR_ex + dt_dZ_ex*dpsidZ_ex)/dpsi_dr_ex 
#dtdr_st_ex = nperiod_data_extend(dtdr_st[1], nperiod)
dt_st_l_ex = nperiod_data_extend(dt_st_l[1], nperiod, par='e')
#pdb.set_trace()
#aprime = (dqdr*theta_st_new + np.reshape(qfac,(-1,1))*dtdr_st)/drhodpsi

#plt.plot(theta, np.interp(theta_comn_mag_ax[1], theta_comn_mag_ax_new[1],aprime[1])); plt.show()

# gradpar = b.grad(theta) with theta = flux theta
gradpar_ex = -1/(R_ex*B_ex)*(dpsi_dr_ex)*(dt_st_l_ex/dl_ex) 

#gradpar with theta = colocation theta
gradpar_col_ex = -1/(R_ex*B_ex)*(dpsi_dr_ex)*(nperiod_data_extend(derm(theta_col, 'l', 'o')[0], nperiod)/dl_ex) 



aprime_bish = -R_ex*B_p_ex*(a_s*dFdpsi[1] +b_s*dpdpsi - c_s)/(2*np.abs(drhodpsi))
#pdb.set_trace()
#plt.plot(theta, np.interp(theta_comn_mag_ax[1], theta_comn_mag_ax_new[1],aprime_bish)); plt.show()
#dtdr_st = diffrho/psi_diff*(aprime_bish - dqdr*theta_st_new)/np.reshape(qfac, (-1,1))

gds21 = diffq/diffrho*(-dpsi_dr_ex)*aprime_bish

dtdr_st_ex = (aprime_bish*drhodpsi - dqdr*theta_st_new_ex)/np.reshape(qfac, (-1,1))

#plt.plot(theta, np.interp(theta_comn_mag_ax[1], theta_comn_mag_ax_new[1],dtdr_st[1]))

gds2 =  (psi_diff/diffrho)**2*(1/R_ex**2 + (dqdr*theta_st_new_ex)**2 + \
        (np.reshape(qfac,(-1,1)))**2*(dtdr_st_ex**2 + (dt_st_l_ex/dl_ex)**2)+ 2*np.reshape(qfac,(-1,1))*dqdr*theta_st_new_ex*dtdr_st_ex)
#gds2_col =  (psi_diff/diffrho)**2*(1/R_ex**2 + (dqdr*theta_col_ex)**2 + \
#        (np.reshape(qfac,(-1,1)))**2*(dtdr_st_ex**2 + (dt_st_l_ex/dl_ex)**2)+ 2*np.reshape(qfac,(-1,1))*dqdr*theta_st_new_ex*dtdr_st_ex)

#plt.plot(theta, np.interp(theta_comn_mag_ax[1], theta_comn_mag_ax_new[1], gds2[1]))

#plt.figure()
gbdrift0 =  1/(B2_ex**2)*psi_diff[1]/diffrho[1]*F[1]/R_ex*(dqdr[1]*dB2l_ex/dl_ex)

#plt.plot(theta, np.interp(theta_comn_mag_ax[1], theta_comn_mag_ax_new[1], gbdrift0[1]))
#theta_eq_arc = eqarc_creator(gradpar, theta_st_new_ex)


####################-------------------dBr calculation-----------------------------####################

#We use Miller's equations to find dBdr using the information given on the middle surface. Miller and Bishop subscripts have been used interchangeably
#dBdr_bish = (B_p**2/B*(1/R_c + dpdpsi*R/(B_p) + F*dFdpsi/dpsi_dr) + B_t**2/(R*B)*(np.sin(u_ML) - dFdpsi/F*R*dpsi_dr))
dBdr_bish = B_p_ex/B_ex*(-B_p_ex/R_c_ex + dpdpsi*R_ex - F[1]**2*np.sin(u_ML_ex)/(R_ex**3*B_p_ex))
#dBdr_bish_2 = B_p_ex/B_ex*(B_p_ex/R_c_ex + dpdpsi*R_ex - F[1]**2*np.sin(u_ML_ex)/(R_ex**3*B_p_ex))
dBdr = dBdr_bish

jac0 =  1

gbdrift = 1/np.abs(drhodpsi*B_ex**3)*(2*B2_ex*dBdr/dpsi_dr_ex + 1/jac0*aprime_bish*drhodpsi*F[1]/R_ex*dB2l_ex/dl_ex*1/B_ex)
#gbdrift = dpsidrho*(-2/B_ex*dBdr_bish/dpsi_dr_ex + 2*aprime*F/R_ex*1/B_ex**3*dBl_ex/dl_ex)

#cvdrift =  1/np.abs(drhodpsi*B**3)*(2*B2*(dpdr/(2*B) + dBdr)/dpsi_dr + 1/jac*aprime_bish*drhodpsi*F/R*dB2l/dl*1/B)
cvdrift = 1/np.abs(drhodpsi*B_ex**3)*(2*B_ex*dpdpsi) + gbdrift 


#pdb.set_trace()
#eqaual-arc theta calculation
gradpar_lim = gradpar_ex[theta_st_new_ex <= np.pi]
B_lim = B_ex[theta_st_new_ex <= np.pi]
B_p_lim = B_p_ex[theta_st_new_ex <= np.pi]
theta_lim = theta_st_new_ex[theta_st_new_ex <= np.pi]
L_eqarc = ctrap(B_p_lim/(B_lim*gradpar_lim), theta_lim, initial=0)
gradpar_eqarc = np.pi/ctrap(1/(gradpar_lim), theta_lim, initial=0)[-1]
#gradpar_eqarc = np.pi/L_eqarc[-1]
#maxval = ctrapz(1/gradpar, theta_st[1])[-1]
#fin_gradpar = np.pi/maxval

theta_eqarc = ctrap(B_lim/B_p_lim*gradpar_eqarc, L_eqarc, initial=0)
theta_eqarc_ex = nperiod_data_extend(theta_eqarc, nperiod, istheta=1)


#pdb.set_trace()
#gradpar_eqarc_ex = gradpar_eqarc*np.ones((len(theta_eqarc_ex),))
#gds21_eqarc = np.interp(theta_eqarc_ex, theta_st_new_ex, gds21[1])
#gds22_eqarc = np.interp(theta_eqarc_ex, theta_st_new_ex, gds22[1])
#gds2_eqarc = np.interp(theta_eqarc_ex, theta_st_new_ex, gds2[1])
#gbdrift0_eqarc = np.interp(theta_eqarc_ex, theta_st_new_ex, gbdrift0)
#B_eqarc = np.interp(theta_eqarc_ex, theta_st_new_ex, B_ex)
#gbdrift_eqarc = np.interp(theta_eqarc_ex, theta_st_new_ex, gbdrift)
#gbdrift0_eqarc = np.interp(theta_eqarc_ex, theta_st_new_ex, gbdrift0)
#cvdrift_eqarc = np.interp(theta_eqarc_ex, theta_st_new_ex, cvdrift)

#gradpar_eqarc_ex = gradpar_col_ex
gradpar_eqarc_ex = gradpar_ex
R_eqarc = R_ex
gds21_eqarc = gds21[1]
gds22_eqarc = gds22[1]
gds2_eqarc = gds2[1]
gbdrift0_eqarc = gbdrift0
B_eqarc = B_ex
gbdrift_eqarc = gbdrift
gbdrift0_eqarc = gbdrift0
cvdrift_eqarc = cvdrift


##theta_st_new_ex_uniq = find_optim_theta_arr(np.vstack((gradpar_eqarc_ex, cvdrift_eqarc, gbdrift_eqarc, gbdrift0_eqarc, B_eqarc, gds2_eqarc, gds21_eqarc, gds22_eqarc)), theta_st_new_ex, 5)
##theta_st_new_ex_uniq = find_optim_theta_arr(np.vstack((gradpar_eqarc_ex, cvdrift_eqarc, gbdrift_eqarc, gbdrift0_eqarc, B_eqarc, gds2_eqarc, gds21_eqarc, gds22_eqarc)), theta_col_ex, 5)
#theta_st_new_ex_uniq = find_optim_theta_arr(np.vstack((gradpar_eqarc_ex, cvdrift, gbdrift, gbdrift0, B_ex, gds2, gds21, gds22)), theta_col_ex, 5)
theta_st_new_ex_uniq = theta_st_new_ex

#plt.plot(theta_st_com_ex, B_ex, '-og', theta_st_com_ex_uniq, B_ex_uniq, '-or', ms=2.2); plt.show() 
#plt.plot(theta_st_com_ex, gbdrift[2], '-og', theta_st_com_ex_uniq, gbdrift_uniq, '-or', ms=2.2); plt.show() 

gradpar_uniq = np.interp(theta_st_new_ex_uniq, theta_st_new_ex, gradpar_eqarc_ex) 
R_uniq = np.interp(theta_st_new_ex_uniq, theta_st_new_ex, R_ex)
cvdrift_uniq = np.interp(theta_st_new_ex_uniq, theta_st_new_ex, cvdrift_eqarc) 
gbdrift_uniq = np.interp(theta_st_new_ex_uniq, theta_st_new_ex, gbdrift_eqarc) 
gbdrift0_uniq = np.interp(theta_st_new_ex_uniq, theta_st_new_ex, gbdrift0_eqarc) 
B_ex_uniq =  np.interp(theta_st_new_ex_uniq, theta_st_new_ex, B_eqarc)
gds2_uniq =  np.interp(theta_st_new_ex_uniq, theta_st_new_ex, gds2_eqarc)
gds21_uniq =  np.interp(theta_st_new_ex_uniq, theta_st_new_ex, gds21_eqarc)
gds22_uniq =  np.interp(theta_st_new_ex_uniq, theta_st_new_ex, gds22_eqarc)

#print("no. of points reduced %.2f%%"%(100*(len(theta_st_new_ex)-len(theta_st_new_ex_uniq))/len(theta_st_new_ex)))
#print("||B_ex_uniq-B_ex|| = %.4f"%(np.linalg.norm(np.interp(theta_st_new_ex, theta_st_new_ex_uniq, B_ex_uniq)-B_ex)))
#print("||gbdrift_uniq-gbdrift|| = %.4f"%(np.linalg.norm(np.interp(theta_st_new_ex, theta_st_new_ex_uniq, gbdrift_uniq)-gbdrift)))
#print("||gds2_uniq-gds2|| = %.4f"%(np.linalg.norm(np.interp(theta_st_new_ex, theta_st_new_ex_uniq, gds2_uniq)-gds2)))


eps = 0.8E-15
#pdb.set_trace()
theta1 = theta_st_new_ex_uniq[theta_st_new_ex_uniq >=0-eps]
theta1 = theta1[theta1 <= np.pi+eps]

theta2 = theta_st_new_ex_uniq[theta_st_new_ex_uniq >= np.pi-eps]
theta2 = theta2[theta2 <= 2*np.pi+eps]

if  len(theta1) > len(theta2):
    #np.sort(theta1)[-1]-np.pi == 0:
    theta_st_new_uniq_sym = theta1
else:
    theta2 = np.abs(theta2 - 2*np.pi)[::-1]
    theta_st_new_uniq_sym = theta2

#theta2 = np.abs(theta2 - 2*np.pi)[::-1]
#theta_st_new_uniq_sym = np.sort(np.unique(np.concatenate((theta1, theta2))))
#pdb.set_trace()

theta_st_new_ex_uniq_sym = nperiod_data_extend(theta_st_new_uniq_sym, nperiod, istheta=1)

eps = 0.2E-15
theta1 = theta_st_new_ex_uniq_sym[theta_st_new_ex_uniq_sym >=0-eps]
theta1 = theta1[theta1 <= np.pi+eps]

theta2 = theta_st_new_ex_uniq_sym[theta_st_new_ex_uniq_sym >= np.pi-eps]
theta2 = theta2[theta2 <= 2*np.pi+eps]

if np.linalg.norm((theta1 + theta2[::-1])/(2*np.pi)-1) > 0:
    print("warning...something is wrong with the theta symmetrization operation")


####################################################################################
###############----------LOCAL SHEAR PLOT----------------###########################
####################################################################################


#plt.plot(theta_eqarc, -B_p_ex*derm(aprime_bish/(R_ex*B_p_ex), 'l', 'o')[0]/dl[1], theta_eqarc, aprime_bish); plt.show()

#plt.plot(theta_eqarc, -B_p[1]*derm(aprime_bish/(R[1]*B_p[1]), 'l', par='o')[0]/dl[1]); plt.show()
#plt.plot(theta_comn_mag_ax_new[1], B_p[1], theta_comn_mag_ax_new[1], B[1])
#plt.show()

#plt.plot(theta_comn_mag_ax_new[1], -B_p[1]*derm(aprime_bish/(R[1]*B_p[1]), 'l', par='o')[0]/dl[1]); plt.show()
#plt.plot(theta_eqarc, -B_p_ex*derm(aprime_bish/(R_ex*B_p_ex), 'l', 'o')[0]/dl[1], theta_eqarc, aprime_bish); plt.show()
#plt.plot(theta_comn_mag_ax_new[1], -B_p_ex*derm(aprime_bish/(R_ex*B_p_ex), 'l', 'o')[0]/dl[1], theta_comn_mag_ax_new[1], aprime_bish); plt.show()
#print('integrated local shear = %.4f'%(-(aprime_bish*2/(R_ex*B_p_ex)*rho[1]/qinp*1/(2*np.pi))[-1]))

# Check is the integrated local shear matches the global shear
print((-(aprime_bish[-1]/(R[1][-1]*B_p[1][-1]))+(aprime_bish[0]/(R[1][0]*B_p[1][0])))*rho[1]/qinp*1/np.pi)

#pdb.set_trace()

a_N = 1
B_N = 1

######################------BALLOONING CALCULATION IN FLUX COORDINATES----------------##############


diff = 0.
one_m_diff = 1. - diff

#gradpar_ex = gradpar_ex[theta_st_new_ex <= np.pi]
#gradpar_ex = gradpar_ex[theta_st_new_ex >=0]
gradpar_ball = reflect_n_append(gradpar_ex, 'e')

theta_ball = reflect_n_append(theta_st_new_ex_uniq, 'o')
#cvdrift = cvdrift[theta_st_new_ex <= np.pi]
#cvdrift = cvdrift[theta_st_new_ex >= 0]
cvdrift_ball = reflect_n_append(cvdrift_uniq, 'e')
R_ball = reflect_n_append(R_uniq, 'e')
B_ball = reflect_n_append(B_ex_uniq, 'e')
#gds2_1 = gds2[1][theta_st_new_ex <=np.pi]
gds2_ball = reflect_n_append(gds2_uniq, 'e')
#gds2 = gds2[theta_st_new_ex >=0]
#gds2_ball = reflect_n_append(gds2_1, 'e')

theta_ball2 = np.linspace(-(2*nperiod-1)*np.pi, (2*nperiod-1)*np.pi, 1001)
gradpar_ball = np.interp(theta_ball2, theta_ball, gradpar_ball)
cvdrift_ball = np.interp(theta_ball2, theta_ball, cvdrift_ball)
R_ball = np.interp(theta_ball2, theta_ball, R_ball)
B_ball = np.interp(theta_ball2, theta_ball, B_ball)
gds2_ball = np.interp(theta_ball2, theta_ball, gds2_ball)



delthet = np.diff(theta_ball2)
ntheta = len(theta_ball2)

#average of dbetadrho
#dbetadrho = -np.sum((cvdrift-gbdrift)*bmag**2)/ntheta 

#dbetadrho = beta_prime_input


#g = gds2_ball/((R_ball*B_ball)**2)
#c = dpdrho*cvdrift_ball*R_ball**2*qfac[1]**2/(F[1]**2)
g = gds2_ball*gradpar_ball/B_ball
c = 1*dpdrho*cvdrift_ball/(B_ball*gradpar_ball)


ch = np.zeros((ntheta,))
gh = np.zeros((ntheta,))

for i in np.arange(1, ntheta):
    ch[i] = 0.5*(c[i] + c[i-1]) 
    gh[i] = 0.5*(g[i] + g[i-1])

cflmax = np.max(np.abs(delthet**2*ch[1:]/gh[1:]))


c1 = np.zeros((ntheta,))
for ig in np.arange(1, ntheta-1):
    c1[ig] = -delthet[ig]*(one_m_diff*c[ig]+0.5*diff*ch[ig+1])\
             -delthet[ig-1]*(one_m_diff*c[ig]+0.5*diff*ch[ig])\
             -delthet[ig-1]*0.5*diff*ch[ig]
    c1[ig]=0.5*c1[ig]


c2 = np.zeros((ntheta,))
g1 = np.zeros((ntheta,))
g2 = np.zeros((ntheta,))

for ig in np.arange(1, ntheta):
    c2[ig] = -0.25*diff*ch[ig]*delthet[ig-1]
    g1[ig] = gh[ig]/delthet[ig-1]
    g2[ig] = 1.0/(0.25*diff*ch[ig]*delthet[ig-1]+gh[ig]/delthet[ig-1])


psi = np.zeros((ntheta,))
psi[0] = 0.
#psi[1]=delthet[0]+0.0001
psi[1]=delthet[0]+0.1
#psi_prime=psi[1]/g2[1]
psi_prime = 0.01

for ig in np.arange(1,ntheta-1):
    psi_prime=psi_prime+c1[ig]*psi[ig]+c2[ig]*psi[ig-1]
    psi[ig+1]=(g1[ig+1]*psi[ig]+psi_prime)*g2[ig+1]

for ig in np.arange(1,ntheta-1):
    if(psi[ig]*psi[ig+1] <= 0 ):
           print("instability detected... please choose a different equilibrium")


plt.plot(theta_ball2, psi); plt.show()


######################-------------BALLOONING CALCULATION---------------------##############

def check_ball(shat_n, dpdpsi_n):

        dpsidrho = 1/drhodpsi
        #dFdpsi_n = (-shat_n*2*np.pi*qfac[1]/(rho[1]*dpsidrho) - b_s[-1]*dpdpsi + c_s[-1])/a_s[-1]
        dFdpsi_n = (-shat_n/(rho[1]*(psi_diff[1]/diffrho[1])*(1/(2*np.pi*qfac[1]*(2*nperiod-1))))-(b_s[-1]*dpdpsi - c_s[-1]))/a_s[-1]
        #dFdpsi_n = -(shat_n/(rho[1]*dpsidrho*(1/(2*np.pi*qfac[1])))-(b_s[-1]*dpdpsi_n - c_s[-1]))/a_s[-1]
        dqdpsi_n = shat_n*qfac[1]/rho[1]*1/dpsidrho
        dqdr_n = dqdpsi_n*dpsi_dr
        
        aprime_n = -R_ex*B_p_ex*(a_s*dFdpsi_n +b_s*dpdpsi_n - c_s)*0.5/np.abs(drhodpsi)
        dtdr_st_n = (aprime_n*drhodpsi - dqdr_n*theta_st_new_ex)/qfac[1]
        #dtdr_st_n = (aprime_n - dqdr_n*theta_st_com)/q_vmec_half[rel_surf_idx]
        gds2_n =  (dpsidrho)**2*(1/R**2 + (dqdr_n*theta_st_new_ex)**2 + (np.reshape(qfac,(-1,1)))**2*(dtdr_st_n**2 + (dt_st_l/dl)**2)+ 2*np.reshape(qfac,(-1,1))*dqdr_n*theta_st_new_ex*dtdr_st_n)
        # BAD PRACTICE F[1]
        dBdr_bish_n = B_p_ex/B_ex*(-B_p_ex/R_c_ex + dpdpsi_n*R_ex - F[1]**2*np.sin(u_ML_ex)/(R_ex**3*B_p_ex))
        gbdrift_n = 1/np.abs(drhodpsi*B**3)*(2*B2*dBdr_bish_n/dpsi_dr + aprime_n*drhodpsi*F[1]/R*dB2l/dl*1/B)
        dPdr_n = dpdpsi_n*dpsi_dr
        #cvdrift_n =  dpsidrho/np.abs(B)*(-2*(2*dPdr_n/(2*B))/dpsi_dr) + gbdrift
        cvdrift_n =  1/np.abs(drhodpsi*B**3)*(2*B2*(2*dPdr_n/(2*B))/dpsi_dr) + gbdrift_n
        #pdb.set_trace()
        #plt.plot(theta_st_new_ex, cvdrift_n[1], theta_st_new[0], cvdrift)
        
        
        diff = 0.
        one_m_diff = 1. - diff
        
        
        theta_ball = reflect_n_append(theta_st_new[1], 'o')
        #theta_ball = reflect_n_append(theta_eqarc, 'o')
        #theta_ball = reflect_n_append(theta_col_ex, 'o')
        cvdrift_ball = reflect_n_append(cvdrift_n[1], 'e')
        #gradpar_ball = reflect_n_append(gradpar_eqarc_ex, 'e')
        #gradpar_ball = reflect_n_append(gradpar_0[1], 'e')
        gradpar_ball = reflect_n_append(gradpar_ex, 'e')
        R_ball = reflect_n_append(R_ex, 'e')
        B_ball = reflect_n_append(B_ex, 'e')
        gds2_ball = reflect_n_append(gds2_n[1], 'e')
        #ntheta = len(theta_ball)
        
        #delthet = np.diff(theta_ball)
        theta_ball2 = np.linspace(-(2*nperiod-1)*np.pi, (2*nperiod-1)*np.pi, 901)
        gradpar_ball = np.interp(theta_ball2, theta_ball, gradpar_ball)
        cvdrift_ball = np.interp(theta_ball2, theta_ball, cvdrift_ball)
        R_ball = np.interp(theta_ball2, theta_ball, R_ball)
        B_ball = np.interp(theta_ball2, theta_ball, B_ball)
        gds2_ball = np.interp(theta_ball2, theta_ball, gds2_ball)
        
        delthet = np.diff(theta_ball2)
        ntheta = len(theta_ball2)
        #average of dbetadrho
        #dbetadrho = -np.sum((cvdrift-gbdrift)*bmag**2)/ntheta 
        
        #g = gds2_ball/((R_ball*B_ball)**2)
        #c = -dpdpsi_n*dpsidrho*cvdrift_ball*qfac[1]**2/(F[1]/R_ball)**2
        #c = -dpdpsi_n*cvdrift_ball*qfac[1]**2/(F[1]/R_ball)**2
        
        g = gds2_ball*gradpar_ball/B_ball
        c = dpdrho*cvdrift_ball/(B_ball*gradpar_ball)

        ch = np.zeros((ntheta,))
        gh = np.zeros((ntheta,))
        
        for i in np.arange(1, ntheta):
            ch[i] = 0.5*(c[i] + c[i-1]) 
            gh[i] = 0.5*(g[i] + g[i-1])
        
        cflmax = np.max(np.abs(delthet**2*ch[1:]/gh[1:]))
        
        
        c1 = np.zeros((ntheta,))
        for ig in np.arange(1, ntheta-1):
            c1[ig] = -delthet[ig]*(one_m_diff*c[ig]+0.5*diff*ch[ig+1])\
                     -delthet[ig-1]*(one_m_diff*c[ig]+0.5*diff*ch[ig])\
                     -delthet[ig-1]*0.5*diff*ch[ig]
            c1[ig]=0.5*c1[ig]
        
        
        c2 = np.zeros((ntheta,))
        g1 = np.zeros((ntheta,))
        g2 = np.zeros((ntheta,))
        
        for ig in np.arange(1, ntheta):
            c2[ig] = -0.25*diff*ch[ig]*delthet[ig-1]
            g1[ig] = gh[ig]/delthet[ig-1]
            g2[ig] = 1.0/(0.25*diff*ch[ig]*delthet[ig-1]+gh[ig]/delthet[ig-1])
        
        
        psi_dum = np.zeros((ntheta,))
        psi_dum[0] = 0.
        psi_dum[1]=delthet[0]
        psi_prime=psi_dum[1]/g2[1]
        
        for ig in np.arange(1,ntheta-1):
            psi_prime=psi_prime+c1[ig]*psi_dum[ig]+c2[ig]*psi_dum[ig-1]
            psi_dum[ig+1]=(g1[ig+1]*psi_dum[ig]+psi_prime)*g2[ig+1]
        
        isunstable = 0
        for ig in np.arange(1,ntheta-1):
            if(psi_dum[ig]*psi_dum[ig+1] <= 0 ):
                isunstable = 1
                #print("instability detected... please choose a different equilibrium")
        return isunstable

shat = s_hat_input
#print("isunstable for the original eq =",check_ball(shat, dpdpsi))



want_a_scan =0
if want_a_scan == 1:
	len1 = 100
	len2 = 100
	#shat_grid = np.linspace(-shat*0.5, 2*shat, len1)
	shat_grid = np.linspace(2, 20, len1)
	#shat_grid = np.linspace(0.10, 0.25, len1)
	#dp_dpsi_grid = np.linspace(-dpdpsi*0.0001, 0.005*dpdpsi, len2)
	#dp_dpsi_grid = np.linspace(-dpdpsi*0.2, 5*dpdpsi, len2)
	dp_dpsi_grid = np.linspace(0, 2, len2)
	#dp_dpsi_grid = np.linspace(0.08, 0.12, len2)
	ball_scan_arr = np.zeros((len1, len2))
	#for i in range(len1):
	#	for j in range(len2):
	#		ball_scan_arr[i,j] = check_ball(shat_grid[i], dp_dpsi_grid[j]) 


	tic5 = time.perf_counter()
	nop = int(48)
	pool = mp.Pool(processes=nop)
	tic5 = time.perf_counter()
	results = np.array([[pool.apply_async(check_ball, args=(shat_grid[i], dp_dpsi_grid[j])) for i in range(len1)] for j in range(len2)])
	for i in range(len2):
		ball_scan_arr[:, i] = np.array([results[i, j].get() for j in range(len1)])
	
	toc5 = time.perf_counter()
	print("Time taken %f"%(toc5-tic5))

	from matplotlib import pyplot as plt
	
	X, Y = np.meshgrid(shat_grid, dp_dpsi_grid)
	Z = ball_scan_arr
	#ax = plt.axes(projection='3d')
	cs = plt.contour(X.T, Y.T*np.abs(dpsidrho)*2, ball_scan_arr, levels=[0.])  
	data = cs.collections[0].get_paths()[0].vertices
	plt.plot(data[:, 0], data[:, 1], '-og', ms=1.5)
	plt.plot(shat, dpdpsi, 'or', ms=3.0)  
	plt.xlabel('shat', fontsize=14)
	plt.ylabel('dbetadrho', fontsize=14)
	plt.title('betaprimeinp=%.3f, s_hat_input=%.2f, theta0=0'%(dpdpsi*np.abs(dpsidrho)*2, s_hat_input))
	del cs
	#ax.plot3d(shat, dpdpsi, 1, 'or')
	#ax.contour3D(X.T, Y.T, ball_scan_arr, 1, cmap='binary')
	plt.show()

pdb.set_trace()



