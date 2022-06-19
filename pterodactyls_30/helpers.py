import numpy as np
from lightkurve.utils import running_mean

"""
Function (adapted from lightkurve) to calculate rms over a given window (in hours)
"""
def calculate_rms(detrended_lc, transit_duration = 7):
	if not isinstance(transit_duration, int):
		raise ValueError("transit_duration must be an integer in units "
					"number of cadences, got {}.".format(transit_duration))
	mean = running_mean(data = detrended_lc, window_size = transit_duration)
	rms_ppm = np.std(mean) * 1e6
	return rms_ppm	


"""
Finds nearest value in array
"""
def find_nearest(array, value):
	array = np.asarray(array)
	idx = (np.abs(array - value)).argmin()
	return idx		


"""
Creates a Period-Radius grid for the injections, with n injections per bin
"""
def PeriodRadius_grid(Pbins, Rpbins, n_per_bin):
	nP = len(Pbins)-1
	nR = len(Rpbins)-1
	#n_inj= nP*nR*n_per_bin
	
	logP_left = np.log(Pbins[:-1])
	logP_right = np.log(Pbins[1:])
	logRp_left = np.log(Rpbins[:-1])
	logRp_right = np.log(Rpbins[1:])
	
	# newaxis magic
	logP_3D = np.random.uniform(logP_left[np.newaxis, :, np.newaxis], 
					   logP_right[np.newaxis, :, np.newaxis], 
					   (nR,nP, n_per_bin) )
	logRp_3D = np.random.uniform(logRp_left[:, np.newaxis, np.newaxis], 
					   logRp_right[:, np.newaxis, np.newaxis], 
					   (nR,nP, n_per_bin) )
	
	P_injs = np.exp(logP_3D.flatten())
	Rp_injs = np.exp(logRp_3D.flatten())

	P0_injs = np.random.uniform(0, P_injs)
	
	return P_injs, Rp_injs, P0_injs