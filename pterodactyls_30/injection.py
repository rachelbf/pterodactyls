import numpy as np
import batman
from .pterodactyls import dactyl
from .lightcurve import LightCurve
from .vetting import *

cgs_au= 1.49597870691e13
cgs_Rsun= 6.960e10
cgs_Rearth= 6.37814e8

class Injection:
	"""dactyl with injected planet in lightcurve"""
	def __init__(self, dactyl, P, Rp, P0, ecc= 0, b=0.5, 
					Contamination=False, Verbose=True, NoRstar=False):
		self.P= P
		# self.Rp comes later, in case it needs the stellar radius
		self.P0= P0
		self.Verbose=Verbose
		self.Injection=True

		# copy paramters from dactyl object
		assert type(dactyl) is dactyl
		self.tic_id= dactyl.tic_id
		self.radius= dactyl.radius
		self.mass= dactyl.mass
		self.ab= dactyl.ab
		#self.lc= dactyl.lc
		self.lctype= dactyl.lctype
		self.sectors= dactyl.sectors
		
		''' Now set Rp'''
		if NoRstar:
			# Planet size units are Rp/Rstar
			self.rp_rs= Rp
			self.Rp= self.rp_rs * (self.radius*cgs_Rsun) / cgs_Rearth
		else:
			# Planet size in earth radii
			self.Rp= Rp
			self.rp_rs= (Rp*cgs_Rearth)/(self.radius*cgs_Rsun)

		if hasattr(dactyl, 'Tmag'):
			self.Tmag = dactyl.Tmag

		# copy rotation rate and window lengths, if present
		if hasattr(dactyl, 'rot_rate'):
			self.rot_rate = dactyl.rot_rate
			self.rotation = 1/self.rot_rate
		if hasattr(dactyl, 'max_splines'):
			self.max_splines = dactyl.max_splines		
		if hasattr(dactyl, 'window_length'):
			self.window_length = dactyl.window_length

		if hasattr(dactyl, 'fluxratio'):
			self.fluxratio = dactyl.fluxratio
		else:
			self.fluxratio = 1.
			if Contamination:
				print ('Warning: no flux ratio to calculate contamination')	

		# copy detected planet parameters from search
		if hasattr(dactyl, 'vetting_tests_passed'):
			self.known_signal= dactyl.vetting_tests_passed
			if hasattr(dactyl, 'vet_flags'):
				self.known_signal_flags= dactyl.vet_flags
			if hasattr(dactyl, 'planet'):
				self.known_signal_period= dactyl.planet['P']

		else:
			self.known_signal=None

		# Calculate inclination from impact parameter 
		# inc= ...

		# injection signal with batman
		ma = batman.TransitParams()
		ma.t0 = P0  # time of inferior conjunction; first transit is X days after start
		ma.per = P  # orbital period
		ma.rp = self.rp_rs # planet radius (in units of stellar radii)
		ma.a = self.mass**(1./3.)* (ma.per/365.)**(2./3.) * cgs_au/(self.radius*cgs_Rsun)  # semi-major axis (in units of stellar radii)
		ma.inc = 90  # orbital inclination (in degrees)
		ma.ecc = 0  # eccentricity
		ma.w = 90  # longitude of periastron (in degrees)
		ma.u = self.ab  # limb darkening coefficients
		ma.limb_dark = "quadratic"  # limb darkening model

		# inject into all sectors
		m = batman.TransitModel(ma, dactyl.lc.time, 
			supersample_factor=21, exp_time=0.020417)  # initializes model
		injection= m.light_curve(ma)

		if Contamination:
			contam_injection= injection*self.fluxratio + (1-self.fluxratio)
			flux = dactyl.lc.flux* contam_injection
		else:
			flux = dactyl.lc.flux* injection

		self.lc= LightCurve(dactyl.lc.time, flux, injection=injection,
			lctype=dactyl.lc.lctype, sector=dactyl.lc.sector)

	def rotation_rate(self, **kwargs):
		if not hasattr(self,'rot_rate'):
			# call detrending routine from main class
			print ('Warning: you are recalculating rotation rates for each injection')
			dactyl.rotation_rate(self, **kwargs)

	def detrend_lightcurve(self, **kwargs):

		# do not recalculate window length / max splines
		if hasattr(self, 'max_splines') and not 'max_splines' in kwargs:
			kwargs['max_splines']= self.max_splines
		if hasattr(self, 'window_length') and not 'window_length' in kwargs:
			kwargs['window_length']= self.window_length
			
		# call detrending routine from main class
		dactyl.detrend_lightcurve(self, **kwargs)

	def search(self, **kwargs):
		# call search routine from main class
		dactyl.search(self, **kwargs)
		
	def vetting(self, **kwargs):
		# call pre vetting routine from vetting.py
		vetting(self, **kwargs)

	def recover(self, Verbose=True):
		# is the signal recovered?
		if not hasattr(self,'vetting_tests_passed'):
			raise ValueError('Run planet search first')

		# check if an injected signal is recovered at same period
		self.rec_P= np.isclose(self.planet['P'], self.P, rtol=0.01)

		if self.rec_P:
			if self.vetting_tests_passed:	
				self.rec=True
				self.status= 'Recovered'
				if Verbose:
					print ('Recovered injection!')

			elif self.TCE:
				self.rec=False
				self.status= 'Failed pre-vet'
				if Verbose:
					print ('Signal Detected but failed pre-vetting:')
					print (self.vet_flags)
			else:
				self.rec=False
				self.status= 'Failed TCE'
				if Verbose:
					print ('Signal detected but below threshold:')
					print ('  SDE= {:.1f}'.format(self.planet['SDE']))
					print ('  SNR= {:.1f}'.format(self.planet['snr']))
		else:
			if self.vetting_tests_passed:
				self.rec=False
				self.status= 'Other Signal'
				if Verbose:
					print ('Stronger signal detected at other period')
					print ('P_recover= {:.2f}'.format(self.planet['P']) )
					print ('P_inject= {:.2f}'.format(self.P) )


				if self.known_signal is not None:
					print ('Known signal present in lightcurve')
					print (self.known_signal)
					if hasattr(self, 'self.known_signal_period'):
						print ('P_known= {:.2f}'.format(self.known_signal_period) )
					if hasattr(self, 'self.known_signal_flags'):
						print (self.known_signal_flags)
						
			else:
				self.rec=False
				self.status= 'No Period Match'
				if Verbose:
					print ('No signal detected at any period')
	
		return self.rec

