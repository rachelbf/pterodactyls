import numpy as np
import os
from datetime import datetime, timedelta
from IPython.utils import io

from .pterodactyls import dactyl
from .injection import Injection
from .vetting import vetting, flux_contamination
from . import helpers, plotter

class IRgrid:
	''' Injection grid using TRICERATOPS '''	
	def __init__(self, ticID, 
				 preset=None, name=None, 
				 Contamination=True, NoRstar=False,
				 Pbins=None, Rpbins=None, n_per_bin=None):
		
		self.ticID= ticID

		self.Contamination= Contamination
		self.NoRstar= NoRstar
		
		if preset == 'small':
			self.name='small' if name is None else name
			self.Pbins=np.geomspace(1,27, 4)
			if NoRstar:
				self.Rpbins= np.geomspace(0.01, 0.3, 5)
			else:
				self.Rpbins= np.geomspace(1,16, 5)
			self.n_per_bin=1
		elif preset == 'large':
			self.name='large' if name is None else name
			self.Pbins=np.geomspace(1,27, 4)
			if NoRstar:
				self.Rpbins= np.geomspace(0.01, 0.3, 5)
			else:
				self.Rpbins= np.geomspace(1,16, 5)
			self.n_per_bin=100
		else:
			self.name= 'custom' if name is None else name
			self.Pbins= np.asarray(Pbins)
			self.Rpbins= np.asarray(Rpbins)
			self.n_per_bin= 1 if n_per_bin is None else n_per_bin
		
		self.nP= self.Pbins.size-1
		self.nRp= self.Rpbins.size-1
		self.n= self.nP*self.nRp*self.n_per_bin
		
		self.path= 'Injections/{}'.format(self.name)
		if not os.path.exists(self.path): 
			os.makedirs(self.path)
		
		if NoRstar:
			print('Injecting Rp_Rs')
			self.fgrid= '{}/RpRs-{}.npz'.format(self.path, self.ticID)
		else:
			self.fgrid= '{}/{}.npz'.format(self.path, self.ticID)
		
		print ('TIC {}: '.format(self.ticID), end='')
		if os.path.isfile(self.fgrid):
			print ('Grid found in {}'.format(self.fgrid))
		else:
			print ('No grid found')
		
	def run(self, Overwrite=False, Plot=True, Load_from_file=True, Bugfix=False):
		if os.path.isfile(self.fgrid) and not Overwrite:
			#print ('Skipping, grid already in {}'.format(self.fgrid))
			if Load_from_file:
				self.load()
			else:
				print ('Skipperino\n')
			return
		
		print ('Running {}x{}x{} grid'.format(self.nP, 
											  self.nRp, 
											  self.n_per_bin))
		self.P_injs, self.Rp_injs, self.P0_injs= \
			helpers.PeriodRadius_grid(self.Pbins, self.Rpbins, 
									  self.n_per_bin)

		
		print ('Initial Search...')
		try:
			tstart= datetime.now()
			with io.capture_output() as captured:
				self.tic_id= search(self.ticID)
			self.dt_search= (datetime.now()-tstart).total_seconds()
			dt_round= timedelta(seconds=round(self.dt_search))
			print ('  Done in {}'.format(dt_round))
		except:
			print ('Initial Searched Failed')
			if Bugfix:
				raise
			else:
				return
		
		injs= []
		print ('Injections:')
		tstart= datetime.now()
		for i,(P_inj, Rp_inj, P0_inj) in enumerate(zip(self.P_injs, 
													   self.Rp_injs, 
													   self.P0_injs)):
			#print ('{}/{}... '.format(i+1, len(self.P_injs), end=''))
			# replace with tqdm?
			print ('  {}/{} @ {}'.format(i+1, self.n, 
								datetime.now().replace(microsecond=0) ))
			try:
				#tloop= datetime.now()
				with io.capture_output() as captured:
					inj= injection(self.tic_id, P_inj, Rp_inj, P0_inj, 
						Contamination=self.Contamination, NoRstar=self.NoRstar)
				#dtloop= (datetime.now()-tloop).total_seconds()
				#print ('{}'.format(timedelta(seconds=round(dt_loop))))
				inj.Crash= False
				
			except:
				print ('Injection {} Failed'.format(i))
				if Bugfix:
					raise				
				else:
					# deal with this, make empty injection
					inj= Injection(self.tic_id, 
					  P=P_inj, Rp=Rp_inj, P0=P0_inj, 
					  Contamination=Contamination, NoRstar=self.NoRstar,
					  Verbose=False)

					#if not hasattr(inj, 'TCE'): inj.TCE= None
					inj.TCE= None
					inj.rec= False
					inj.status='Crash'
					inj.Crash= True
		
			injs.append(inj)
		
		self.dt_injections= (datetime.now()-tstart).total_seconds()
		
		self.P_rec= np.array([None if inj.Crash else inj.planet['P'] for inj in injs])
		self.Rp_rec= np.array([None if inj.Crash else inj.planet['Rp'] for inj in injs])
		
		self.rec= np.array([inj.rec for inj in injs])
		#self.rec_P= np.array([inj.rec_P for inj in injs])
		self.status= np.array([inj.status for inj in injs])
		
		#self.TCE= np.array([inj.TCE for inj in injs])
		self.snr= np.array([None if inj.Crash else inj.planet['snr'] for inj in injs])
		self.SDE= np.array([None if inj.Crash else inj.planet['SDE'] for inj in injs])
		
		#self.rec_fraction= self.rec.sum()/self.rec.size
		
		# save this too
		self.rp_rs_injs = np.array([inj.rp_rs for inj in injs])

		if self.NoRstar:
			# make sure this is Rp and not Rp_Rstar
			self.Rp_injs = np.array([inj.Rp for inj in injs])

		np.savez(self.fgrid, P_injs=self.P_injs, Rp_injs=self.Rp_injs,
				P0_injs=self.P0_injs, rp_rs_injs=self.rp_rs_injs,
				P_rec=self.P_rec, Rp_rec=self.Rp_rec,
				rec=self.rec, status=self.status,
				snr=self.snr, SDE=self.SDE,
				dt_search=self.dt_search, 
				dt_injections=self.dt_injections,
				)
		
		print ('Finished in {}'.format(
			timedelta(seconds=round(self.dt_injections))))
		print ('  {} per injection'.format(
			timedelta(seconds=round(self.dt_injections/self.nP)) ))
		
		self.rec_fraction= self.rec.sum()/self.rec.size
		print ('  {:.0%} recovered'.format(self.rec_fraction))
		
		if Plot:
			print ('Making plot')
			plot.injection_per_star(self, Save=True, Show=False)
		
		print ()
	
	def load(self):
		npz=np.load(self.fgrid, allow_pickle=True) # allow for None entries
		
		self.P_injs= npz['P_injs']
		self.Rp_injs= npz['Rp_injs']
		self.P0_injs= npz['P0_injs']
		self.rp_rs_injs= npz['rp_rs_injs']
		
		self.P_rec= npz['P_rec']
		self.Rp_rec= npz['Rp_rec']
		self.rec= npz['rec']
		self.status= npz['status']
		self.snr= npz['snr']
		self.SDE= npz['SDE']

		self.dt_search= npz['dt_search'][()]
		self.dt_injections= npz['dt_injections'][()]
		
#		 self.= npz['']
#		 self.= npz['']

		npz.close()
	
		self.rec_fraction= self.rec.sum()/self.rec.size
		
		print ('Total runtime was {}, {:.0%} recovered'.format(
			timedelta(seconds=round(self.dt_injections)),
			self.rec_fraction))
		
		return self
		
	def get_Tmag(self, path='pap1_obs', Verbose=False):
		
		npz=np.load('{}/{}/sectors.npz'.format(path, self.ticID) )
		if 'Tmag' in npz:
			self.Tmag= npz['Tmag'][()]
			if Verbose:
				print ('{}: Tmag {}'.format(self.ticID, self.Tmag))
		else:
			self.Tmag= None
			if Verbose:
				print ('{}: no Tmag'.format(self.ticID))
		npz.close()
		
		return self.Tmag

	def get_fluxratio(self, path='pap1_obs', Verbose=False):
		
		npz=np.load('{}/{}/fratio.npz'.format(path, self.ticID) )
		if 'fratio' in npz:
			self.fluxratio= npz['fratio'][()]
			if Verbose:
				print ('{}: flux ratio {}'.format(self.ticID, self.fluxratio))
		else:
			self.fluxratio= None
			if Verbose:
				print ('{}: no flux ratio'.format(self.ticID))
		npz.close()
		
		return self.fluxratio

def search(ticID):
	tic= TIC(ticID, Verbose=False, save_on_disk=True)
	
	tic.extract_lightcurve(save_on_disk=True, Verbose=False)
	flux_contamination(tic, Verbose = False, doPlot = False)
	tic.clean_lc(Verbose = False, DoPlot = False)
	tic.detrend_lightcurve(method = 'penalized', Verbose = False)
	tic.search(doPlot=False)
	
	vetting(tic, Verbose=False)
	
	return tic

def injection(tic, P, Rp, P0, Contamination=True, NoRstar=False):
	inject= Injection(tic, 
					  P=P, Rp=Rp, P0=P0, 
					  Contamination=Contamination, NoRstar=NoRstar,
					  Verbose=False)
	inject.detrend_lightcurve(method = 'penalized', Verbose=False)
	inject.search(doPlot=False, Verbose=False)
	vetting(inject, Verbose=True)
	inject.recover(Verbose=False)
	
	return inject