#!/usr/bin/env python3
"""
This script contains all the functions that are called by the main script local_eikcoefs_gen.py
"""


import numpy as np

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

def dermv(arr, brr, ch, par='e'):
    # Finite difference subroutine
    # brr is the independent variable arr. Needed for weighted finite-difference
    # ch = 'l' means difference along the flux surface
    # ch = 'r' mean difference across the flux surfaces
    # par = 'e' means even parity of the arr. PARITY OF THE INPUT ARRAY 
    # par = 'o' means odd parity
    #pdb.set_trace()
    temp = np.shape(arr)
    if len(temp) == 1 and ch == 'l': #finite diff along the flux surface for a single array
        if par == 'e':
                d1, d2 = np.shape(arr)[0], 1
                arr = np.reshape(arr, (d2,d1))
                brr = np.reshape(brr, (d2,d1))
                diff_arr = np.zeros((d2,d1))
                diff_arr[0, 0] =  0. #(arr_theta_-0 - arr_theta_+0)  = 0
                diff_arr[0, -1] = 0. 
                #diff_arr[0, 1:-1] = np.diff(arr[0,:-1], axis=0) + np.diff(arr[0,1:], axis=0)  
                for i in range(1, d1-1):
                    h1 = (brr[0, i+1] - brr[0, i])
                    h0 = (brr[0, i] - brr[0, i-1])
                    diff_arr[0, i] =  (arr[0, i+1]/h1**2 + arr[0, i]*(1/h0**2 - 1/h1**2) - arr[0, i-1]/h0**2)/(1/h1 + 1/h0) 
        else:
                #pdb.set_trace()
                d1, d2 = np.shape(arr)[0], 1
                arr = np.reshape(arr, (d2,d1))
                brr = np.reshape(brr, (d2,d1))
                diff_arr = np.zeros((d2,d1))
                
                h1 = (np.abs(brr[0, 1]) - np.abs(brr[0, 0]))
                h0 = (np.abs(brr[0, -1]) - np.abs(brr[0, -2]))
                diff_arr[0, 0] =  (4*arr[0, 1]-3*arr[0, 0]-arr[0,2])/(2*(brr[0, 1]-brr[0,0]))

                #diff_arr[0, -1] = (-4*arr[0,-1]+3*arr[0, -2]+arr[0, -3])/(2*(brr[0, -1]-brr[0, -2]))
                diff_arr[0, -1] = (-4*arr[0,-2]+3*arr[0, -1]+arr[0, -3])/(2*(brr[0, -1]-brr[0, -2]))
                #diff_arr[0, -1] = 2*(arr[0, -1] - arr[0, -2])/(2*(brr[0, -1] - brr[0, -2])) 
                #diff_arr[0, 1:-1] = np.diff(arr[0,:-1], axis=0) + np.diff(arr[0,1:], axis=0)  
                for i in range(1, d1-1):
                    h1 = (brr[0, i+1] - brr[0, i])
                    h0 = (brr[0, i] - brr[0, i-1])
                    diff_arr[0, i] =  (arr[0, i+1]/h1**2 + arr[0, i]*(1/h0**2 - 1/h1**2) - arr[0, i-1]/h0**2)/(1/h1 + 1/h0) 

    elif len(temp) == 1 and ch == 'r': # across surfaces for a single array
        pdb.set_trace()
        d1, d2 = np.shape(arr)[0], 1
        diff_arr = np.zeros((d1,d2))
        arr = np.reshape(arr, (d1,d2))
        diff_arr[0, 0] = 2*(arr[1, 0] - arr[0, 0])/(2*(brr[1, 0] - brr[0, 0])) # single dimension arrays like psi, F and q don't have parity
        diff_arr[-1, 0] = 2*(arr[-1, 0] - arr[-2, 0])/(2*(brr[-1, 0] - brr[-2, 0])) 
        #diff_arr[1:-1, 0] = np.diff(arr[:-1,0], axis=0) + np.diff(arr[1:,0], axis=0)  
        for i in range(1, d1-1):
            h1 = (brr[i+1, 0] - brr[i, 0])
            h0 = (brr[i, 0] - brr[i-1, 0])
            diff_arr[i, 0] =  (arr[i+1, 0]/h1**2 - arr[i, 0]*(1/h0**2 - 1/h1**2) - arr[i-1, 0]/h0**2)/(1/h1 + 1/h0) 



    else:
        d1, d2 = np.shape(arr)[0], np.shape(arr)[1]
        
        diff_arr = np.zeros((d1,d2))
        if ch == 'r': # across surfaces for multi-dim array
                #pdb.set_trace()
                diff_arr[0, :] = 2*(arr[1,:] - arr[0,:])/(2*(brr[1, :] - brr[0, :])) 
                diff_arr[-1, :] = 2*(arr[-1,:] - arr[-2,:])/(2*(brr[-1, :] - brr[-2, :])) 
                #diff_arr[1:-1, :] = (np.diff(arr[:-1,:], axis=0) + np.diff(arr[1:,:], axis=0))  
                for i in range(1, d1-1):
                    h1 = (brr[i+1, :] - brr[i, :])
                    h0 = (brr[i, :] - brr[i-1, :])
                    diff_arr[i, :] =  (arr[i+1, :]/h1**2 + arr[i, :]*(1/h0**2 - 1/h1**2) - arr[i-1, :]/h0**2)/(1/h1 + 1/h0) 
            
        else: #along a surface for a multi-dim array
            #pdb.set_trace()
            if par == 'e':
                #pdb.set_trace()
                diff_arr[:, 0] = np.zeros((d1,))
                diff_arr[:, -1] = np.zeros((d1,))
                #diff_arr[:, 1:-1] = (np.diff(arr[:,:-1], axis=1) + np.diff(arr[:,1:], axis=1))  
                for i in range(1, d2-1):
                    h1 = (brr[:, i+1] - brr[:, i])
                    h0 = (brr[:, i] - brr[:, i-1])
                    diff_arr[:, i] =  (arr[:, i+1]/h1**2 + arr[:, i]*(1/h0**2 - 1/h1**2) - arr[:, i-1]/h0**2)/(1/h1 + 1/h0) 
                    #pdb.set_trace()
            else:
                #pdb.set_trace()
                diff_arr[:, 0] = 2*(arr[:, 1] - arr[:, 0])/(2*(brr[:, 1] - brr[:, 0])) 
                diff_arr[:, -1] = 2*(arr[:, -1] - arr[:, -2])/(2*(brr[:, -1] - brr[:, -2]))
                #diff_arr[:, 1:-1] = (np.diff(arr[:,:-1], axis=1) + np.diff(arr[:,1:], axis=1))  
                for i in range(1, d2-1):
                    h1 = (brr[:, i+1] - brr[:, i])
                    h0 = (brr[:, i] - brr[:, i-1])
                    diff_arr[:, i] =  (arr[:, i+1]/h1**2 + arr[:, i]*(1/h0**2 - 1/h1**2) - arr[:, i-1]/h0**2)/(1/h1 + 1/h0) 
        
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
