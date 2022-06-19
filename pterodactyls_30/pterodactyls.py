import numpy as np
import os
import eleanor
from eleanor.utils import *
from wotan import t14, flatten
from wotan.gaps import get_gaps_indexes
from transitleastsquares import (
	transitleastsquares,
	cleaned_array,
	catalog_info,
	transit_mask
	)

from .helpers import *
from .lightcurve import LightCurve
from astropy.timeseries import LombScargle
import matplotlib.pyplot as plt
import pandas as pd
import math
from .plotter import *

xls = 'Pterodactyls_30mins.xlsx'
psf_fixes = pd.read_excel(xls, sheet_name = 'eleanor_psf_masking', header = 0)
high_amp = pd.read_excel(xls, sheet_name = 'High_amp_sig3_cut', header = 0)
sem_mask = pd.read_excel(xls, sheet_name = 'SEM_mask', header = 0)
offset_fix = pd.read_excel(xls, sheet_name = 'offset_fix', header = 0)
disco_fix = pd.read_excel(xls, sheet_name = 'discontinuity_fix', header = 0)
mask_fix = pd.read_excel(xls, sheet_name = 'mask_fix', header = 0)


class dactyl:
# Uses TIC ID and TLS to get stellar radii, stellar mass, and limb darkening parameters, when available for a given target. 
# Has the option to save it all to disk for optimization.
# Args:
#     tic_id (int): TIC ID of the target.
#     cluster (str): Folder name to save Pterodactyl Object of Interest (POI) plots; default: NoCluster
#     Verbose (bool): Prints out warnings and any other relevant statements; default: False
#     save_on_disk (bool): Saves stellar parameters on disk for optimization; default: False 
#     Overwrite (bool): Overwrites the lightcurve, stellar radii, stellar mass, and limb darkening parameters that were saved to disk; default: False
	def __init__(self, tic_id, cluster = 'NoCluster', Verbose = False, 
		save_on_disk = False, Overwrite = False):
		self.tic_id = tic_id
		self.cluster = cluster        
		self.Verbose = Verbose
		self.path = 'targets/{}/'.format(tic_id)
		self.fratio_path = 'flux_ratios/{}/'.format(tic_id)
		self.Injection = False
		self.candidate_path = 'POIs/{}/{}/'.format(cluster,tic_id)

		if save_on_disk:
			if not os.path.exists(self.path): os.makedirs(self.path)
			fstar = self.path + 'star.npz'

		# read previously saved data
		if save_on_disk and os.path.isfile(fstar) and not Overwrite:
			if Verbose: print('(load {})'.format(fstar))
			npz = np.load(fstar)
			self.mass = npz['mass'][()] # array -> float
			self.radius = npz['radius'][()]
			self.ab = npz['ab']
			npz.close()
		else:

			# get stellar properties using TLS function
			self.ab, self.mass, mass_min, mass_max, self.radius, radius_min, radius_max = \
				catalog_info(TIC_ID = tic_id)

			''' do some checks on mass and radius '''	
			if np.isfinite(self.radius) == False:
				print('WARNING! TIC {}'.format(self.tic_id) + ' does not have a stellar radius value!')
				self.radius = 1.0
           			
			if np.isfinite(self.mass) == False:
				print('WARNING! TIC {}'.format(self.tic_id) + ' does not have a stellar mass value!')
				self.mass = 1.0
				
			if np.isfinite(self.ab[0]) == False & np.isfinite(self.ab[1]) == False:
				print('WARNING! TIC {}'.format(self.tic_id) + ' does not have a stellar limb darkening parameters!')
				self.ab[0] = 0.4804 #a G2V star in the Kepler bandpass
				self.ab[1] = 0.1867

			if save_on_disk:
				np.savez_compressed(fstar, mass = self.mass, radius = self.radius, ab = self.ab)
				print ('(saved {})'.format(fstar))

		if self.Verbose:
			print ('Init TIC {}'.format(self.tic_id))
			print ('  R_star: {:.3g} Rsun'.format(self.radius))
			print ('  M_star: {:.3g} Msun'.format(self.mass))
			print ('  Limb darkening parameters (a,b): {}'.format(self.ab))

	def lightcurve_from_disk(self, Verbose = True):

        # Returns light curves for each sector stored on disk.
        # Args:
        #     Verbose (bool): Prints out warnings and any other relevant statements; default: True

		# read in which sectors have data stored
		if Verbose: print('(load {})'.format(self.path + 'sectors.npz'))
		npz = np.load(self.path + 'sectors.npz')
		
		self.sectors = npz['sectors']
		self.nsectors = len(self.sectors)
		if self.nsectors==1: self.sector= self.sectors[0]

		if 'Tmag' in npz:
			self.Tmag = npz['Tmag'][()]

		npz.close()

		# read in the light curve of each sector
		lcs = {} # dictionary to hold the different curves (raw, cor, pca, psf)
		lcs['apertures'] = []
		lcs['chippos'] = []

		for sector in self.sectors:
			fsector = self.path + 'sector{}.npz'.format(sector)
			if Verbose: print ('load sector {} from tmp/'.format(sector))
			npz = np.load(fsector)
			for key in ['raw','cor','psf','pca']:
				if key+'_flux' in npz:
					if not key in lcs: lcs[key] = []
					lcs[key].append( LightCurve(npz['time'], npz[key+'_flux'], 
										lctype = key, sector = sector ))
			lcs['apertures'].append(npz['aperture'])
			lcs['chippos'].append(npz['chippos'])
			npz.close()

		return lcs

	def lightcurve_from_eleanor(self, do_psf = True, do_pca = True, save_on_disk = False, Verbose = True):

        # """
        # Returns light curves for each sector from eleanor
        # Args:
        #     do_psf (bool): Downloads eleanor's psf light curves for the target; default: True
        #     do_pca (bool): Downloads eleanor's pca light curves for the target; default: True
        #     save_on_disk (bool): Saves light curve on disk for optimization; default: False 
        #     Verbose (bool): Prints out warnings and any other relevant statements; default: True
        # """

		print('Extracting LightCurve of TIC {}'.format(self.tic_id))

		# lookup if a lightcurve exists and download target pixel file
		stars = eleanor.multi_sectors(tic = self.tic_id, sectors = 'all', tc = False)
		
		if type(stars) is list:
			self.nsectors = len(stars)
			self.sectors = [s.sector for s in stars]
			self.Tmag= stars[0].tess_mag
		else:
			self.sector = stars.sector
			self.sectors = [self.sector]
			self.nsectors = 1
			self.Tmag= stars.tess_mag

		if save_on_disk:
			if not os.path.exists(self.path): os.makedirs(self.path)
			np.savez_compressed(self.path + 'sectors.npz', 
				sectors = self.sectors, Tmag = self.Tmag)
			print ('(saved {})'.format(self.path + 'sectors.npz'))

		# select TPF
		apt = dict(height = 13, width = 13)

		# extract lightcurve for each sector
		data_all = []
		for star in stars:
			print ('eleanor.TargetData sector {}'.format(star.sector))
			data = eleanor.TargetData(star, do_psf = do_psf, do_pca = do_pca, regressors = 'corner')
			for sector in self.sectors: 
			    if sector == 13:
			        #replace problematic values with Nans
			        data.tpf[540:550]=np.nan
			        data.tpf_err[540:550]=np.nan
			        data.flux_bkg[540:550]=np.nan
                    #create new array that masks Nans
			        mask_tpf = np.ma.masked_where(data.tpf==np.nan, data.tpf)
			        mask_tpf_err = np.ma.masked_where(data.tpf_err==np.nan, data.tpf_err)
			        mask_flux_bkg = np.ma.masked_where(data.flux_bkg==np.nan, data.flux_bkg)
			        #converts from masked array to np.array
			        new_tpf = mask_tpf.filled() 
			        new_tpf_err = mask_tpf_err.filled()
			        new_flux_bkg = mask_flux_bkg.filled()
			        #Pass new stuff to psf function
			        data.psf_lightcurve(new_tpf, new_tpf_err, new_flux_bkg, model='gaussian')

			    s5_mask = psf_fixes['s5_tics'].tolist()
			    if (sector == 5) & (self.tic_id in s5_mask):

 			        value3 = 1463.93945
 			        idx_3 = find_nearest(data.time, value3) # Camera 1 Guiding Disabled

 			        value4 = 1464.40056
 			        idx_4 = find_nearest(data.time, value4) #end of segment 2
                    
 			        #replace problematic values with Nans
                    
 			        data.tpf[idx_3:idx_4]=np.nan
 			        data.tpf_err[idx_3:idx_4]=np.nan
 			        data.flux_bkg[idx_3:idx_4]=np.nan
                    
                    #create new array that masks Nans
 			        mask_tpf = np.ma.masked_where(data.tpf==np.nan, data.tpf)
 			        mask_tpf_err = np.ma.masked_where(data.tpf_err==np.nan, data.tpf_err)
 			        mask_flux_bkg = np.ma.masked_where(data.flux_bkg==np.nan, data.flux_bkg)
 			        #converts from masked array to np.array
 			        new_tpf = mask_tpf.filled() 
 			        new_tpf_err = mask_tpf_err.filled()
 			        new_flux_bkg = mask_flux_bkg.filled()
 			        #Pass new stuff to psf function
 			        data.psf_lightcurve(new_tpf, new_tpf_err, new_flux_bkg, model='gaussian')  

			    s8_mask = psf_fixes['s8_tics'].tolist()
			    if (sector == 8) & (self.tic_id in s8_mask):
 			        value1 = 1517.34150
 			        idx_1 = find_nearest(data.time, value1) # start of segment 1

 			        value2 = 1517.39566
 			        idx_2 = find_nearest(data.time, value2) #before guiding camera was enabled for segment 1

 			        value3 = 1530.25816
 			        idx_3 = find_nearest(data.time, value3) # start of segment 2

 			        value4 = 1535.00264
 			        idx_4 = find_nearest(data.time, value4) #before guiding camera was enabled for segment 2
                    
 			        #replace problematic values with Nans
 			        data.tpf[idx_1+1:idx_2]=np.nan
 			        data.tpf_err[idx_1+1:idx_2]=np.nan
 			        data.flux_bkg[idx_1+1:idx_2]=np.nan
                    
 			        data.tpf[idx_3:idx_4]=np.nan
 			        data.tpf_err[idx_3:idx_4]=np.nan
 			        data.flux_bkg[idx_3:idx_4]=np.nan
                    
                    #create new array that masks Nans
 			        mask_tpf = np.ma.masked_where(data.tpf==np.nan, data.tpf)
 			        mask_tpf_err = np.ma.masked_where(data.tpf_err==np.nan, data.tpf_err)
 			        mask_flux_bkg = np.ma.masked_where(data.flux_bkg==np.nan, data.flux_bkg)
 			        #converts from masked array to np.array
 			        new_tpf = mask_tpf.filled() 
 			        new_tpf_err = mask_tpf_err.filled()
 			        new_flux_bkg = mask_flux_bkg.filled()
 			        #Pass new stuff to psf function
 			        data.psf_lightcurve(new_tpf, new_tpf_err, new_flux_bkg, model='gaussian')                     

			    s9_mask = psf_fixes['s9_tics'].tolist()                   
			    if (sector == 9) & (self.tic_id in s9_mask):
 			        value1 = 1543.21648
 			        idx_1 = find_nearest(data.time, value1) # start of segment 1

 			        value2 = 1543.75080
 			        idx_2 = find_nearest(data.time, value2) #before guiding camera was enabled for segment 1

 			        value3 = 1556.72344
 			        idx_3 = find_nearest(data.time, value3) # start of segment 2

 			        value4 = 1557.00080
 			        idx_4 = find_nearest(data.time, value4) #before guiding camera was enabled for segment 2
                    
 			        #replace problematic values with Nans
 			        data.tpf[idx_1+1:idx_2]=np.nan
 			        data.tpf_err[idx_1+1:idx_2]=np.nan
 			        data.flux_bkg[idx_1+1:idx_2]=np.nan
                    
 			        data.tpf[idx_3:idx_4]=np.nan
 			        data.tpf_err[idx_3:idx_4]=np.nan
 			        data.flux_bkg[idx_3:idx_4]=np.nan
                    
                    #create new array that masks Nans
 			        mask_tpf = np.ma.masked_where(data.tpf==np.nan, data.tpf)
 			        mask_tpf_err = np.ma.masked_where(data.tpf_err==np.nan, data.tpf_err)
 			        mask_flux_bkg = np.ma.masked_where(data.flux_bkg==np.nan, data.flux_bkg)
 			        #converts from masked array to np.array
 			        new_tpf = mask_tpf.filled() 
 			        new_tpf_err = mask_tpf_err.filled()
 			        new_flux_bkg = mask_flux_bkg.filled()
 			        #Pass new stuff to psf function
 			        data.psf_lightcurve(new_tpf, new_tpf_err, new_flux_bkg, model='gaussian')                    

			    s10_mask = psf_fixes['s10_tics'].tolist()                    
			    if (sector == 10) & (self.tic_id in s10_mask):
 			        value1 = 1569.43176
 			        idx_1 = find_nearest(data.time, value1) # start of segment 1

 			        value2 = 1570.87620
 			        idx_2 = find_nearest(data.time, value2) #before guiding camera was enabled for segment 1

 			        value3 = 1582.76231
 			        idx_3 = find_nearest(data.time, value3) # start of segment 2

 			        value4 = 1584.72342
 			        idx_4 = find_nearest(data.time, value4) #before guiding camera was enabled for segment 2
                    
 			        #replace problematic values with Nans
 			        data.tpf[idx_1+1:idx_2]=np.nan
 			        data.tpf_err[idx_1+1:idx_2]=np.nan
 			        data.flux_bkg[idx_1+1:idx_2]=np.nan
                    
 			        data.tpf[idx_3:idx_4]=np.nan
 			        data.tpf_err[idx_3:idx_4]=np.nan
 			        data.flux_bkg[idx_3:idx_4]=np.nan
                    
                    #create new array that masks Nans
 			        mask_tpf = np.ma.masked_where(data.tpf==np.nan, data.tpf)
 			        mask_tpf_err = np.ma.masked_where(data.tpf_err==np.nan, data.tpf_err)
 			        mask_flux_bkg = np.ma.masked_where(data.flux_bkg==np.nan, data.flux_bkg)
 			        #converts from masked array to np.array
 			        new_tpf = mask_tpf.filled() 
 			        new_tpf_err = mask_tpf_err.filled()
 			        new_flux_bkg = mask_flux_bkg.filled()
 			        #Pass new stuff to psf function
 			        data.psf_lightcurve(new_tpf, new_tpf_err, new_flux_bkg, model='gaussian')
                     
			    s11_mask = psf_fixes['s11_tics'].tolist()
			    if (sector == 11) & (self.tic_id in s11_mask):
 			        value1 = 1596.77203
 			        idx_1 = find_nearest(data.time, value1) # start of segment 1

 			        value2 = 1599.94148
 			        idx_2 = find_nearest(data.time, value2) #before guiding camera was enabled for segment 1

 			        value3 = 1610.77620
 			        idx_3 = find_nearest(data.time, value3) # start of segment 2

 			        value4 = 1614.19842
 			        idx_4 = find_nearest(data.time, value4) #before guiding camera was enabled for segment 2
                    
 			        #replace problematic values with Nans
 			        data.tpf[idx_1+1:idx_2]=np.nan
 			        data.tpf_err[idx_1+1:idx_2]=np.nan
 			        data.flux_bkg[idx_1+1:idx_2]=np.nan
                    
 			        data.tpf[idx_3:idx_4]=np.nan
 			        data.tpf_err[idx_3:idx_4]=np.nan
 			        data.flux_bkg[idx_3:idx_4]=np.nan
                    
                    #create new array that masks Nans
 			        mask_tpf = np.ma.masked_where(data.tpf==np.nan, data.tpf)
 			        mask_tpf_err = np.ma.masked_where(data.tpf_err==np.nan, data.tpf_err)
 			        mask_flux_bkg = np.ma.masked_where(data.flux_bkg==np.nan, data.flux_bkg)
 			        #converts from masked array to np.array
 			        new_tpf = mask_tpf.filled() 
 			        new_tpf_err = mask_tpf_err.filled()
 			        new_flux_bkg = mask_flux_bkg.filled()
 			        #Pass new stuff to psf function
 			        data.psf_lightcurve(new_tpf, new_tpf_err, new_flux_bkg, model='gaussian')                     
                
			data_all.append(data)
		self.data_all =data_all
		
		# create a dictionary to hold the different curves (raw, cor, pca, psf)
		lcs = {}
		lcs['raw'] = []
		lcs['cor'] = []
		if do_psf: lcs['psf'] = []
		if do_pca: lcs['pca'] = []

		lcs['apertures'] = []
		lcs['chippos'] = []

		# store each sector separately
		for data, sector in zip(data_all, self.sectors):
			q = (data.quality == 0)# | (data.quality == 4096) | (data.quality == 128)
			lcs['raw'].append( LightCurve(data.time[q], data.raw_flux[q], 
											lctype = 'raw', sector = sector) )
			lcs['cor'].append( LightCurve(data.time[q], data.corr_flux[q], 
											lctype = 'cor', sector = sector) )
			if do_psf: lcs['psf'].append( LightCurve(data.time[q], data.psf_flux[q], 
											lctype = 'psf', sector = sector) )
			if do_pca: lcs['pca'].append( LightCurve(data.time[q], data.pca_flux[q], 
											lctype = 'pca', sector = sector) )

			lcs['apertures'].append(np.asarray(data.aperture)) 
			lcs['chippos'].append([data.header['CHIPPOS1'], data.header['CHIPPOS2']])

			if save_on_disk:
				    fsector = self.path + 'sector{}.npz'.format(sector)
				    print ('(saved {})'.format(fsector))
				    fluxes =dict(raw_flux = data.raw_flux[q], cor_flux = data.corr_flux[q])					
				    if do_psf: fluxes['psf_flux'] = data.psf_flux[q]
				    if do_pca: fluxes['pca_flux'] = data.pca_flux[q]
				    np.savez_compressed(fsector, time = data.time[q], 
					aperture = np.asarray(data.aperture), 
					chippos = [data.header['CHIPPOS1'],data.header['CHIPPOS2']],
					**fluxes)

		return lcs

	def extract_lightcurve(self, lctype = 'psf',
		do_psf = True, do_pca = True, store_all_lightcurves = True,
		save_on_disk = True, Overwrite = False, 
		NormalizeSegments = True, Verbose = True):

   #      """
   #      Extracts all lightcurves for each sector and stores the combined curve in .lc
   #      Args:
   #          lctype (str): What type of eleanor light curve will be used for analysis i.e. 'psf', 'pca', 'raw' or 'corr'; default: 'psf' (recommended)
   #          do_psf (bool): Downloads eleanor's psf light curves for the target; default: True
   #          do_pca (bool): Downloads eleanor's pca light curves for the target; default: True
			# store_all_lightcurves (bool): Stores each lightcurve type in each sector; default: True
   #          save_on_disk (bool): Saves stellar parameters on disk for optimization; default: False 
   #          Overwrite (bool): Overwrites the lightcurves that were saved to disk; default: False 
   #          NormalizeSegments (bool): Normalize all segments to median flux. Segments are parts of light curves sepearated by breaks in the data; default: True
   #          Verbose (bool): Prints out warnings and any other relevant statements; default: True         
   #      """

		# read lightcurves from disk for faster processing
		if save_on_disk and os.path.isfile(self.path + 'sectors.npz') and not Overwrite:
			lcs = self.lightcurve_from_disk(Verbose = Verbose)
		else:
			lcs = self.lightcurve_from_eleanor(save_on_disk = save_on_disk,
				do_psf = do_psf, do_pca = do_pca, Verbose = Verbose)

		# Normalize all segments to median flux        
		combined_lc = np.concatenate([lc.flux for lc in lcs[lctype]])
		combined_time = np.concatenate([lc.time for lc in lcs[lctype]])
		if NormalizeSegments:
			segs = get_gaps_indexes(combined_time, 0.5)
			segs[-1] -= 2  # remove endpoint padding
			stops = []
			starts = []
			new_flux= []
			new_time = []
			for seg in range(len(segs)-1):
				start = combined_time[segs[seg]]
				stop = combined_time[segs[seg+1]-1]
				starts.append(start)
				stops.append(stop)
			for i,j in zip(starts,stops):
				mask = (combined_time>i) & (combined_time<j)
				masked_flux = combined_lc[mask]
				masked_time = combined_time[mask]
				fluxes = masked_flux/np.median(masked_flux)
				new_flux.append(fluxes)
				new_time.append(masked_time)
			flux = np.concatenate(new_flux)
			time = np.concatenate(new_time)
		else:
			flux = combined_lc
			time = combined_time



		# each light curve type in each sector
		if store_all_lightcurves:
			self.lcs = lcs

		# store one combined lightcurve for all sectors
		self.lctype = lctype
		self.lc = LightCurve(
            time = time, 
            flux = flux,         
				lctype = self.lctype, sector = [lc.sector for lc in lcs[lctype]],
				)
		self.apertures = lcs['apertures']
		self.chippos = lcs['chippos']

	def rotation_rate(self, nterms = 2, method = 'auto', minimum_frequency = 0.1, maximum_frequency = 20, samples_per_peak = 5):

        # """
        # Calculates rotation rate from the light curve using astropy's Lomb-Scargle Periodogram; see: https://docs.astropy.org/en/stable/api/astropy.timeseries.LombScargle.html#astropy.timeseries.LombScargle
        # Args:
        # 	nterms (int): Number of terms to use in the Fourier fit; default: 2
        # 	method (str): Specify the lomb scargle implementation to use. Options are: 'auto', 'slow', 'chi2', 'cython', 'fast', 'fastchi2', 'scipy'; default: 'auto'
        # 	minimum_frequency (float): default: 0.1
        # 	maximum_frequency (float): default: 20
        # 	samples_per_peak (int): The approximate number of desired samples across the typical peak; default: 5       
        # """

		LS = LombScargle(self.lc.time, self.lc.flux, nterms = nterms)

		freq, power = LS.autopower(method = method, minimum_frequency = minimum_frequency, maximum_frequency = maximum_frequency, samples_per_peak = samples_per_peak)
		rot_rate = freq[np.argmax(power)] #in 1/days
		
		# store the rotation rate
		self.freq = freq
		self.power = power
		self.rot_rate = rot_rate
		self.rotation = 1/self.rot_rate
		# find the rotation trend
		fit = LS.model(self.lc.time, self.rot_rate)

		# store the rotation trend in the LightCurve object
		self.lc.rot_fit = fit
		self.lc.rot_flat= self.lc.flux / fit

	def detrend_lightcurve(self, method = 'penalized', break_tolerance = 0.5, max_splines = None, stdev_cut = 2, nterms = 2, edge_cutoff = 0.5, window_length = 0.5, Verbose = True):

        # """
        # Detrends the light curve using routines from Wötan (https://wotan.readthedocs.io/en/latest/)
        # Args:
        # 	method (str): Specify which detrending method to use. Options are: 'penalized', 'Huber', 'robust' ; default: 'penalized' (recommended)
        # 	break_tolerance (float): Split light curve into segments at breaks longer than this value in days; default: 0.5
      		# max_splines (int): The maximum number of knots for peanlized spline to use; default: None (since it uses the rotation rate to calculate the max_splines)
      		# stdev_cut (int): data points more than these many standard deviations away from the fit are removed; default: 2 (3 for high amplitude stars)
      		# nterms (int): Number of terms to use in the Fourier fit to calculate the rotation rate; default: 2
      		# edge_cutoff (float): Removes the value (in days) long edges of the light curve segments; default: 0.5
      		# window_length (float): For robust and Huber splines, the knot distance in units of time; default: 0.5
        #     Verbose (bool): Prints out warnings and any other relevant statements; default: True 
        # """

		print('Detrending LightCurve of TIC {}'.format(self.tic_id))

		''' Optimize detrending using rotation rates'''
		if method == 'penalized':

			if max_splines is None:
				''' Calculate rotation rate to set max splines '''
				self.rotation_rate(nterms=nterms)
			
				if self.rotation <= 2.0: 
					if Verbose: print('Fast Rotator')
					self.max_splines = 200
				elif (self.rotation > 2.0) and (self.rotation < 10.0):
					if Verbose: print('Not a fast rotator')
					self.max_splines = int(200/self.rotation)
                    
				elif self.rotation == 10.0:
 					if Verbose: print('Not a fast rotator')
 					self.max_splines = 25

			else:
				''' Manually set number of splines '''
				self.rotation_rate(nterms=nterms)                
				self.max_splines = max_splines

		''' Use better sigma clipping for high amplitude rotators'''
		if method == 'penalized':
			sig3_tics = high_amp['tics'].tolist()
			if self.tic_id in sig3_tics:
				print ('TIC {} has a high amplitude rotation signal! Increasing sigma clipping to 3!'.format(self.tic_id))                
				self.stdev_cut = 3
			else:
				print ('Using a default sigma clipping of 2 for TIC {}!'.format(self.tic_id))                 
				self.stdev_cut = stdev_cut              
		
		self.break_tolerance = break_tolerance
		self.edge_cutoff = edge_cutoff
		self.return_nsplines = True

		# define the detrending function
		if method == 'huber':
			kwargs=dict(method = 'hspline', 
				window_length = self.window_length, 
				break_tolerance = self.break_tolerance, return_trend = True)
			
		elif method == 'penalized':
			kwargs=dict(method = 'pspline', 
				break_tolerance = self.break_tolerance, 
				edge_cutoff = self.edge_cutoff, stdev_cut=self.stdev_cut, max_splines = self.max_splines, 
				return_nsplines = self.return_nsplines, return_trend = True)   

		elif method == 'robust':
			kwargs=dict(method = 'rspline', 
				window_length = self.window_length, 
				break_tolerance = self.break_tolerance, return_trend = True)  
		
		else:
			print ('{} detrending method from Wötan has not been implemented yet. Please contact us and we will work with you to get it running...}'.format(method))

				
		try:
			if method == 'penalized':
				if self.return_nsplines: #for number of knots in pspline
					flat, trend, nsplines = flatten(self.lc.time, self.lc.flux, **kwargs)
					self.nsplines = nsplines	
					if Verbose:	 
						print('Lightcurve was split into: ', len(self.nsplines), 'segments')
						print('Chosen number of splines for each segment: ', self.nsplines)
					# Get segment information
					segs = get_gaps_indexes(self.lc.time, self.break_tolerance)
					segs[-1] -= 2  # remove endpoint padding
					durations = []
					self.durations = durations
					for seg in range(len(segs)-1):
						start = self.lc.time[segs[seg]]
						stop = self.lc.time[segs[seg+1]-1]
						duration = stop - start
						durations.append(duration)
						if Verbose:
						    print('Segment:', seg+1, '\nStart time: ', '{:04.2f}'.format(start), 
					'days \tStop time: ', '{:04.2f}'.format(stop), 
					'days \tSegment Length: ', '{:04.2f}'.format(duration), 'days')
				else:
					flat, trend = flatten(self.lc.time, self.lc.flux, **kwargs)
			else:
				flat, trend = flatten(self.lc.time, self.lc.flux, **kwargs)
		except Exception as err:
			err_message = 'Detrending failed'
			print (err_message)
			print (kwargs)
			raise ValueError(err_message) from err

		# add detrended flux to lc object
		self.lc.add_detrend(flat, trend)		

	def planet_mask(self, period, duration, T0, use_detrended_lc = False, DoPlot = False):

   #      """
   #      Used to mask known planets for multiplanet search
   #      Args:
			# period (float): orbital period of the planet you want to mask
			# duration (float): transit duration of the planet you want to mask. We recommend you mask 2-3 times the transit duration to make sure the entire transit is masked
			# T0 (float): Mid-transit time of the first transit within the time series
			# use_detrended_lc (bool): Masks detrended light curve instead of psf/pca light curve; default: False (recommended)
			# DoPlot (bool): Plots masked region of the psf/pca light curve for visualization; default: True

   #      """        

		if use_detrended_lc:
			mask = transit_mask(self.lc.time_detrend, period = period, 
				duration = duration, T0 = T0)

			self.lc.flat = self.lc.flat[~mask]
			self.lc.time_detrend = self.lc.time_detrend[~mask]
			self.lc.trend_nz = self.lc.trend_nz[~mask]
		else:
			intransit = transit_mask(self.lc.time, period = period, duration = duration, T0 = T0)
			new_time, new_flux = cleaned_array(self.lc.time[~intransit], self.lc.flux[~intransit])
			if DoPlot:
				plt.figure(figsize = (15, 5))
				plt.scatter(self.lc.time[intransit], self.lc.flux[intransit], color = 'red', zorder = 7)  # in-transit points in red
				plt.scatter(new_time, new_flux, color = 'blue')  # other points in blue

			self.lc= LightCurve(time = new_time, flux = new_flux, 
				lctype = self.lc.lctype, sector = self.lc.sector)	

	def mask(self, time_min, time_max, use_detrended_lc = True, DoPlot = False):
   #      """
   #      Mask parts of the light curve
   #      Args:
   #      	time_min (float): Start time location of mask
   #      	time_min (float): End time location of mask
			# use_detrended_lc (bool): Masks detrended light curve instead of psf/pca light curve; default: False (recommended)
			# DoPlot (bool): Plots masked region of the psf/pca light curve for visualization; default: True
   #      """

		if use_detrended_lc:
			mask = (self.lc.time_detrend<time_max) & (self.lc.time_detrend>time_min)
			if DoPlot:
				plt.plot(self.lc.time_detrend, self.lc.flat, 'g.')
				plt.plot(self.lc.time_detrend[mask], self.lc.flat[mask], 'r.', label = 'masked')
				plt.legend(loc = 1, prop = {'size': 16})
				plt.xlabel('Time (in days)', size = 18)
				plt.ylabel('Flux', size = 18)
				plt.show()
				plt.close()

			self.lc.flat = self.lc.flat[~mask]
			self.lc.time_detrend = self.lc.time_detrend[~mask]
			self.lc.trend_nz = self.lc.trend_nz[~mask]

            
		else:
			mask = (self.lc.time<time_max) & (self.lc.time>time_min)
			if DoPlot:
				plt.plot(self.lc.time, self.lc.flux, 'g.')
				plt.plot(self.lc.time[mask], self.lc.flux[mask], 'r.', label = 'masked')
				plt.legend(loc = 1, prop = {'size': 18})
				plt.title('TIC {}'.format(self.tic_id), size = 18)
				plt.xlabel('Time (in days)', size = 18)
				plt.ylabel('Flux', size = 18)
				plt.show()
				plt.close()
			new_time, new_flux = cleaned_array(self.lc.time[~mask], self.lc.flux[~mask])

			self.lc= LightCurve(time = new_time, flux = new_flux, 
				lctype = self.lc.lctype, sector = self.lc.sector)

	def disco_check(self, min_win, max_win):
        # """
        # Check for breaks in the light curve between two points in time
        # Args:
        # 	min_win (float): Start time
        # 	max_win (float): End time
        # """

		check = (self.lc.time > min_win) & (self.lc.time < max_win)
		time_window = self.lc.time[check]
		flux_window = self.lc.flux[check]
		diffs = [t - s for s, t in zip(flux_window, flux_window[1:])]
		max_diffs = max(diffs, key = abs)
		max_idx = (np.where(diffs == max_diffs))
		disc_idx = (max_idx[0][0])
		offset_flux = flux_window[disc_idx]
		offset_time = time_window[disc_idx]
		return np.round(offset_time, decimals = 2)

	def fix_offset(self, offset_loc = None,  mid_start = None, mid_end = None, 
			mid_sector = False, DoPlot = False, Verbose=True):
  #   """
  #   Fixes broken light curves based on known location or unknown location when used in conjunction with disco_check
  #   Args:
  #   	offset_loc (float): Time location of offset; default: None
		# mid_start (float): Start time location of offset if a part in the middle of the light curve is offset; default: None
		# mid_end (float): End time location of offset if a part in the middle of the light curve is offset; default: None
		# mid_sector (bool): Set to True if a part in the middle of the light curve is offset; default: False
		# DoPlot (bool): Plots fixed region of the light curve for visualization; default: True
		# Verbose (bool): Prints out warnings and any other relevant statements; default: True 		
  #   """

		if (offset_loc != None):
 			lc1 = self.lc.time<offset_loc

 			lc2 = self.lc.time>offset_loc
          
 			if lc1.sum()==0:
 				if Verbose: print ('Skipping: No data before flux drop', str(self.tic_id))
 				return
 			elif lc2.sum()==0:
 				if Verbose: print ('Skipping: No data after flux drop', str(self.tic_id))
 				return
 			else: 
	 			#get flux difference
	 			diff = self.lc.flux[lc1][-1] - self.lc.flux[lc2][1]
	 			#fix offset
	 			new_lc1 = self.lc.flux[lc1] - diff
	             
	 			new_flux = np.append(new_lc1, self.lc.flux[lc2])                   
	 			self.lc = LightCurve(time = self.lc.time, flux = new_flux, 
					lctype = self.lc.lctype, sector = self.lc.sector)

		elif mid_sector:
			lc3 = (self.lc.time<mid_start)            
			lc4 = (self.lc.time>mid_start) & (self.lc.time<mid_end)
			lc5 = self.lc.time>mid_end

			if lc3.sum()==0:
				if Verbose: print ('  Skipping: No data before flux drop')
				return
			elif lc4.sum()==0:
				if Verbose: print ('  Skipping: No data between flux drops')
				return
			elif lc5.sum()==0:
				if Verbose: print ('  Skipping: No data after flux drop')
				return
			else: 
				mid_diff = self.lc.flux[lc4][-1] - self.lc.flux[lc5][0]
				new_lc4 = self.lc.flux[lc4] - mid_diff

				new_flux = np.concatenate((self.lc.flux[lc3], new_lc4, self.lc.flux[lc5]))          
				self.lc = LightCurve(time = self.lc.time, flux = new_flux, 
				lctype = self.lc.lctype, sector = self.lc.sector) 

		if DoPlot:
			if (offset_loc != None):
				plt.plot(self.lc.time[lc2], self.lc.flux[lc2], 'g.')
				plt.plot(self.lc.time[lc1], self.lc.flux[lc1], '.', label = 'Original offset')
				plt.plot(self.lc.time[lc1], new_lc1, 'r.', label = 'Fixed offset')
			elif mid_sector:
 			    plt.plot(self.lc.time[lc4], self.lc.flux[lc4], 'g.')
 			    plt.plot(self.lc.time[lc5], self.lc.flux[lc5], '.', label = 'Original offset')
 			    plt.plot(self.lc.time[lc4], new_lc4, 'r.', label = 'Fixed offset')
			plt.legend(loc = 3, prop = {'size': 18})
			plt.title('TIC {}'.format(self.tic_id), size = 18)
			plt.xlabel('Time (in days)', size = 18)
			plt.ylabel('Flux', size = 18)
			plt.show()
			plt.close()
        
	def mask_SEM(self, sigma = 5):
        # """
        # Skye Excess Metric mask
        # """
		if sigma == 5:
			sector_10 = sem_mask['sig5_s10'].tolist()
			for time in sector_10:
				self.mask(time - 0.005, time + 0.005, use_detrended_lc = False, DoPlot = False)

			sector_11 = sem_mask['sig5_s11'].tolist()                
			for time in sector_11:
				self.mask(time - 0.005, time + 0.005, use_detrended_lc = False, DoPlot = False)    
                
			sector_10 = sem_mask['sig5_s12'].tolist()
			for time in sector_12:
				self.mask(time - 0.005, time + 0.005, use_detrended_lc = False, DoPlot = False)     

		elif sigma == 3:
			sector_1 = sem_mask['sig3_s1'].tolist()
			for time in sector_1:
				self.mask(time - 0.005, time + 0.005, use_detrended_lc = False, DoPlot = False)

			sector_2 = sem_mask['sig3_s2'].tolist()
			for time in sector_2:
				self.mask(time - 0.005, time + 0.005, use_detrended_lc = False, DoPlot = False)

			sector_3 = sem_mask['sig3_s3'].tolist()
			for time in sector_3:
				self.mask(time - 0.005, time + 0.005, use_detrended_lc = False, DoPlot = False)

			sector_4 = sem_mask['sig3_s4'].tolist()
			for time in sector_4:
				self.mask(time - 0.005, time + 0.005, use_detrended_lc = False, DoPlot = False)

			sector_10 = sem_mask['sig3_s10'].tolist()
			for time in sector_10:
				self.mask(time - 0.005, time + 0.005, use_detrended_lc = False, DoPlot = False)
                
			sector_11 = sem_mask['sig3_s11'].tolist()
			for time in sector_11:
				self.mask(time - 0.005, time + 0.005, use_detrended_lc = False, DoPlot = False)    
                
			sector_12 = sem_mask['sig3_s12'].tolist()
			for time in sector_12:
				self.mask(time - 0.005, time + 0.005, use_detrended_lc = False, DoPlot = False) 

	def clean_lc(self, Verbose = True, DoPlot = True):
  #       """
  #       Function to get rid of systematics in the light curve using the mask, disco_check and fix_offset functions
		# Args: 
		# Verbose (bool): Prints out warnings and any other relevant statements; default: True 
		# DoPlot (bool): Plots fixed region of the light curve for visualization; default: True
  #       """

		for sector in self.sectors:     
			nans = np.count_nonzero(~np.isfinite(self.lc.flux))
			self.ratio_nan = nans/len(self.lc.flux)
			if nans != 0:
				print('Warning! TIC ' + str(self.tic_id) + ' has ' + str(nans) + ' non-finite and bad data points in sector ' + 
              str(sector) + ' data (' + str(np.round(self.ratio_nan*100, decimals = 2)) + '% of the light curve)')
				print('Removing non-finite and bad data points...')
			finite = np.where(np.isfinite(self.lc.flux) & (self.lc.flux > 0))
			self.lc.flux = self.lc.flux[finite]
			self.lc.time = self.lc.time[finite]
			self.lc= LightCurve(time=self.lc.time, flux=self.lc.flux, 
			lctype=self.lc.lctype, sector= self.lc.sector)

		if self.tic_id in disco_fix['tics'].values:
			min_win = disco_fix.loc[disco_fix['tics'] == self.tic_id,'min_win']
			max_win = disco_fix.loc[disco_fix['tics'] == self.tic_id,'max_win']
			for i, j in zip(min_win, max_win):
				if Verbose:				
					print('Offset Detected at {} and {} in TIC {}'.format(i,j,self.tic_id))
				loc2 = self.disco_check(i, j)                 
				self.fix_offset(offset_loc = loc2, DoPlot = DoPlot) 

		elif self.tic_id in mask_fix['tics'].values:
			time_min = mask_fix.loc[mask_fix['tics'] == self.tic_id,'time_min']
			time_max = mask_fix.loc[mask_fix['tics'] == self.tic_id,'time_max']	
			for i,j in zip(time_min,time_max):				
				if Verbose:
					print('Masking TIC {} light curve between {} and {}'.format(self.tic_id,i,j)) 
				self.mask(time_min = i, time_max = j, use_detrended_lc = False, DoPlot = DoPlot)

		elif (sector == 1) & (self.tic_id == 410214986):
			if Verbose:
					print('Fixed Offset in lightcurve of TIC {}'.format(self.tic_id))  
			self.fix_offset(offset_loc = 1327.38911, DoPlot = DoPlot)           
			self.fix_offset(mid_start = 1332.78497, mid_end = 1333.11830, mid_sector = True, DoPlot = DoPlot)  
			
		elif sector == 4:    
			if Verbose:
					print('Fixed offset due to an error in the uploaded guidestar tables at 1413.26 days in the Sector 4 lightcurve of TIC {}'.format(self.tic_id))
			self.fix_offset(offset_loc = 1413.26468, DoPlot = DoPlot)  

		elif sector == 12 :  
			if Verbose:
					print('Masking Sector 12 of TIC {}'.format(self.tic_id)) 
			self.mask(time_min = 1624.94979, time_max = 1628.5, use_detrended_lc = False, DoPlot = DoPlot) 

		elif sector == 19:
			if Verbose:
					print('Fixed Offset in lightcurve of TIC {}'.format(self.tic_id))  
			loc1 = self.disco_check(1821, 1822)                 
			self.fix_offset(offset_loc = loc1, DoPlot = DoPlot) 
			loc2 = self.disco_check(1826, 1827)                 
			self.fix_offset(offset_loc = loc2, DoPlot = DoPlot)  

		elif sector == 21:
			if Verbose:
					print('Fixed Offset in lightcurve of TIC {}'.format(self.tic_id))  
			loc2 = self.disco_check(1882.5, 1884)                 
			self.fix_offset(offset_loc = loc2, DoPlot = DoPlot) 
			loc1 = self.disco_check(1876.5, 1877)       
			self.fix_offset(offset_loc = loc1, DoPlot = DoPlot)                             

		elif sector == 22:
			if Verbose:
				print('Fixed Offset in lightcurve of TIC {}'.format(self.tic_id))  
			loc2 = self.disco_check(1905, 1906.5)                 
			self.fix_offset(offset_loc = loc2, DoPlot = DoPlot) 

		else:
			self.lc= LightCurve(time=self.lc.time, flux=self.lc.flux, 
			lctype=self.lc.lctype, sector= self.lc.sector)	

		return self.ratio_nan

	def search(self, n_transits_min = 2, snr_cut = 7, sde_cut = 7, search_injection=False, Verbose=True, doPlot = True):

  #       """
  #       Function to search the lightcurve with TLS based on user defined snr and SDE cuts: https://transitleastsquares.readthedocs.io/en/latest/
		# Args: 
		# n_transits_min (int): ; default: 2
		# snr_cut (float): SNR cut that defines what a Threshold Crossing Event is; default: 7
		# sde_cut (float): SDE cut that defines what a Threshold Crossing Event is; default: 7
		# search_injection (bool): searches the injected light curve instead of the detrended light curve; default: False
		# Verbose (bool): Prints out warnings and any other relevant statements; default: True 
		# doPlot (bool): Plots fixed region of the light curve for visualization; default: True
  #       """


		if search_injection:
			tls = transitleastsquares(self.lc.time, self.lc.injection)
		else:
			tls = transitleastsquares(self.lc.time_detrend, self.lc.flat)
			print ('Searching in TIC {}'.format(self.tic_id) )


		rad_m = [255475779, 318027415, 396137563, 144609763, 390842111, 323245745, 462080292, 162208214, 75540884,
             99039842, 374542527, 185053184, 374812449]
		if self.tic_id in rad_m:
			radius_max = 10
			radius_min = 0.1
			self.radius = 1   
			
		results = tls.power(u = self.ab, R_star = self.radius,
			duration_grid_step = 1.1, period_min = 0.5, period_max = 25, n_transits_min = n_transits_min, oversampling_factor = 3, T0_fit_margin = 0.01)

		# is this a Threshold crossing event?
		self.TCE= (results.SDE > sde_cut) & (results.snr > snr_cut)   
		if not self.Injection:
			if self.TCE:
				print ('\nTLS found a TCE in TIC {}'.format(self.tic_id))
			else:
				print ('Nothing found in TIC {}'.format(self.tic_id) )
				print('*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*')

		# store strongest signal
		p= self.planet = {}
		p['T0'] = results.T0 #(float) Mid-transit time of the first transit within the time series
		p['rp_rs'] = results.rp_rs #(float) Radius ratio of planet and star using the analytic equations from Heller 2019
		p['Rp'] = results.rp_rs * 110. * self.radius #Radius of Planet
		p['P'] = results.period # Period of the best-fit signal
		p['P_err'] = results.period_uncertainty #  Uncertainty of the best-fit period (half width at half maximum)
		p['P0'] = results.T0 % results.period # Period at which first transit is detected??
		p['depth'] = results.depth #(float) Best-fit transit depth (measured at the transit bottom)
		p['tdur_obs'] = results.duration #(float) Best-fit transit duration
		p['periods'] = results.periods #(array) The period grid used in the search
		p['SDE'] = results.SDE
		p['transit_times'] = results.transit_times
		p['snr'] = results.snr
		p['transit_depths'] = results.transit_depths
		p['transit_depths_uncertainties'] = results.transit_depths_uncertainties
		p['power'] = results.power
		p['model_folded_phase'] = results.model_folded_phase
		p['model_folded_model'] = results.model_folded_model
		p['model_lightcurve_time'] = results.model_lightcurve_time
		p['model_lightcurve_model'] = results.model_lightcurve_model     
		p['folded_phase'] = results.folded_phase
		p['folded_y'] = results.folded_y #(array) Data flux of each point
		p['folded_dy'] = results.folded_dy #(array) Data uncertainty of each point
		p['per_transit_count'] = results.per_transit_count
		p['distinct_transit_count'] = results.distinct_transit_count
		p['empty_transit_count'] = results.empty_transit_count    
		p['snr_per_transit'] = results.snr_per_transit	
		p['odd_even_mismatch'] = results.odd_even_mismatch  
		p['FAP'] = results.FAP    
		p['chi2'] = results.chi2 
		p['chi2red'] = results.chi2red
		p['chi2_min'] = results.chi2_min  
		p['chi2red_min'] = results.chi2red_min       
        
		tdur = t14(R_s=self.radius, M_s=self.mass, P=self.planet['P'], small_planet=True)
		p['tdur_exp'] = tdur

		self.tls_results= results
        
		if not self.Injection:
			if self.TCE:
				if not os.path.exists(self.candidate_path): os.makedirs(self.candidate_path)
				if doPlot:
				    lightcurve_tce(self, save_plot=True)
				    planet_periodogram(self, save_plot=True)
				    planet_phasefold(self, save_plot=True)							
