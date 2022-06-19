import numpy as np
from .helpers import calculate_rms

class LightCurve:
	"""Holds a lightcurve

	Attributes
	----------
	lctype : str
		the type of eleanor lightcurve (cor, raw, flux, psf)
	sector : list
		the list of sectors with data
	time : np.arr
		time
	flux : np.arr
		flux at each `time` (from eleanor)
	flux_n : np.arr
		median normalized `flux`
	flux_ppm : np.arr
		median normalized `flux` in ppm
	rms : float
		rms of `flux`
	injection : np.arr, optional
		injected signal
	time_detrend : np.arr
		`time` series for detrended lightcurve with nans removed
	flat : np.arr
		flattened light curve on `time_detrend`
	flat_ppm : np.arr
		median normalized `flat` in ppm
	flat_rms : float
		rms of detrended light curve `flat`
	trend : np.arr
		trend light curve on `time`
	trend_nz : np.arr
		trend light curve on `time_detrend`
	trend_n : np.arr
		median normalized `trend`
	injection_detrend : np.arr, optional
		injected signal on `time_detrend`
	rot_fit : np.arr
		rotation trend on `time`
	rot_flat : np.arr
		light curve on `time` with rotation divided out
	"""
	def __init__(self, time, flux, lctype='unknown', sector=None, injection=None):

		self.lctype = lctype
		self.sector= sector
		self.time = time
		self.flux = flux
		self.flux_n = self.flux/np.median(self.flux)
		self.flux_ppm = (self.flux_n-1)*1e6
		self.rms = calculate_rms(self.flux_n)

		if injection is not None:
			self.injection = injection

	def add_detrend(self, flat, trend, NonZero=True):

		#look for nans
		nz = np.isfinite(flat) & np.isfinite(trend)

		self.time_detrend = self.time[nz]
		
		self.flat_raw = flat
		self.flat = flat[nz]
		self.flat_ppm = (1-self.flat)*1e6
		self.flat_rms = calculate_rms(flat[nz])

		self.trend = trend # trend_raw?
		self.trend_nz = trend[nz]
		#self.trend_n = trend[nz]/np.median(trend[nz])
		self.trend_n = self.trend/np.nanmedian(self.trend)

		if hasattr(self, 'injection'):
			self.injection_detrend = self.injection[nz]

	def add_rotation(self, rotation):

		self.rot_fit = rotation
		self.rot_flat = self.lc/self.rotation