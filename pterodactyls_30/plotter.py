import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from astropy.stats import sigma_clip
import transitleastsquares as tls
import scipy
import math
from transitleastsquares import transit_mask
import math, os


"""
Plots all available eleanor light curves
"""	
def lightcurves_eleanor(dactyl):
	fig, ax = plt.subplots(2, 2, sharex=True, sharey = True)

	for lc in dactyl.lcs['raw']:
		ax[0, 0].plot(lc.time, lc.flux_n, '.', label = 'S{}'.format(lc.sector))
	ax[0, 0].set_title('Raw Flux')
	ax[0, 0].set_ylabel('Normalized Flux')
	ax[0, 0].set_xlabel('Time')

	for lc in dactyl.lcs['cor']:
		ax[1, 0].plot(lc.time, lc.flux_n, '.', label = 'S{}'.format(lc.sector))
	ax[1, 0].set_title('Corr Flux')
	ax[1, 0].set_ylabel('Normalized Flux')
	ax[1, 0].set_xlabel('Time')

	if 'pca' in dactyl.lcs:
		for lc in dactyl.lcs['pca']:
			ax[0, 1].plot(lc.time, lc.flux_n, '.', label = 'S{}'.format(lc.sector))
	ax[0, 1].set_title('PCA Flux')
	ax[0, 1].set_ylabel('Normalized Flux')
	ax[0, 1].set_xlabel('Time')

	if 'psf' in dactyl.lcs:
		for lc in dactyl.lcs['psf']:
			ax[1, 1].plot(lc.time, lc.flux_n, '.', label = 'S{}'.format(lc.sector))
	ax[1, 1].set_title('PSF Flux')
	ax[1, 1].set_ylabel('Normalized Flux')
	ax[1, 1].set_xlabel('Time')

	plt.legend(loc='upper center', bbox_to_anchor=(-0.1, -0.25), fancybox=True, shadow=True, ncol=5)
	plt.subplots_adjust(top=0.92, bottom=0.08, left=0.10, right=0.95, hspace=0.35, wspace=0.35)
	
	plt.show()

	plt.show()
	plt.close()

"""
Plots psf light curve for a given target, with adjustable axes
"""
def plot_psf(dactyl, xmin = None, xmax = None, ymin = None, ymax = None):
	plt.figure(figsize=(10,5))
	dactyl.sector_list = ','.join(str(v) for v in dactyl.sectors)
	if 'psf' in dactyl.lcs:
		for lc in dactyl.lcs['psf']:
 			plt.title('TIC {}'.format(dactyl.tic_id), size = 20)  
 			plt.plot(lc.time, lc.flux_n, 'b.')
 			plt.ylabel('PSF Flux', size = 18)
 			plt.xlabel('Time (in days)', size = 18)
 			
 			text = 'Sectors: ' + dactyl.sector_list
 			# plt.text(0.01, 0.9, text, fontweight="normal", size = 18)
 			if xmin != None: plt.xlim(left = xmin)
 			if xmax != None: plt.xlim(right = xmax)
             
 			if ymin != None: plt.ylim(bottom = ymin)
 			if ymax != None: plt.ylim(top = ymax)             

	plt.show()
	plt.close()

"""
Plots psf and pca light curve for a given target, with adjustable axes
"""
def plot_lc(dactyl, xmin = None, xmax = None, ymin = None, ymax = None):
	fig, (ax1, ax2) = plt.subplots(2, 1, figsize = (10, 7), sharex = True)
	dactyl.sector_list = ','.join(str(v) for v in dactyl.sectors)
	if 'psf' in dactyl.lcs:
		for lc in dactyl.lcs['psf']:
 			ax1.set_title('TIC {}'.format(dactyl.tic_id), size = 20)  
 			ax1.plot(lc.time, lc.flux_n, 'b.')
 			ax1.set_ylabel('PSF Flux', size = 18)
 			
 			text = 'Sectors: ' + dactyl.sector_list
 			ax1.text(0.01, 0.9, text, fontweight="normal", size = 18, transform=ax1.transAxes)
 			if xmin != None: ax1.set_xlim(left = xmin)
 			if xmax != None: ax1.set_xlim(right = xmax)
				 
 			if ymin != None: ax1.set_ylim(bottom = ymin)
 			if ymax != None: ax1.set_ylim(top = ymax)
			
	dactyl.pca_time = []
	dactyl.pca_flux = []
	if 'pca' in dactyl.lcs:
		for lc in dactyl.lcs['pca']:
			ax2.plot(lc.time, lc.flux_n, 'r.')
			ax2.set_ylabel('PCA Flux', size = 18)
			ax2.set_xlabel('Time (in days)', size = 18)
			if xmin != None: ax2.set_xlim(left = xmin)
			if xmax != None: ax2.set_xlim(right = xmax)			

			if ymin != None: ax2.set_ylim(bottom = ymin)
			if ymax != None: ax2.set_ylim(top = ymax)
	
	fig.tight_layout()

	plt.show()
	plt.close()
	

"""
Plots rotation figures for a given target
"""	
def plot_rotation(dactyl):
	fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize = (10, 10))
	plt.subplots_adjust(hspace = 0.5)
	text = 'Rotation Rate: ' + str(format(dactyl.rotation, '.2f')) + ' days'
	st = plt.suptitle(text)
	st.set_y(1.01)  
	ax1.plot(1/dactyl.freq, dactyl.power, "k")
	ax1.set_xscale('log')
	ax1.set_ylabel("Power")
	ax1.set_title('TIC {}'.format(dactyl.tic_id))
	ax1.set_xlabel('Period (days)')

	ax2.plot(dactyl.lc.time, dactyl.lc.flux, '.')
	ax2.plot(dactyl.lc.time, dactyl.lc.rot_fit, '-r')	
	ax2.set_title('Stellar Rotation Model')
	ax2.set_xlabel('Time (in days)')
	ax2.set_yticks([])
	ax2.set_ylabel('Flux')
		
	ax3.plot(dactyl.lc.time, dactyl.lc.rot_flat, '.g')
	ax3.set_title('Light Curve after removing stellar rotation')
	ax3.set_xlabel('Time (in days)')
	ax3.set_ylabel('Flux')

	fig.tight_layout()

	plt.show()
	plt.close()

"""
Plots detrended light curve for a given target, with adjustable axes
"""
def lightcurve_detrended(dactyl, xmin = None, xmax = None, ymin = None, ymax = None):
	dactyl.sector_list = ','.join(str(v) for v in dactyl.sectors)
	fig, (ax1, ax2) = plt.subplots(2, 1, figsize = (10,7), sharex = True)
	plt.subplots_adjust(hspace = 0.01)
	
	ax1.scatter(dactyl.lc.time, dactyl.lc.flux, color = 'blue', marker = '.', alpha = 1, label = 'PSF Flux')
	ax1.plot(dactyl.lc.time, dactyl.lc.trend, color = 'red', linewidth = 2.0, linestyle = 'solid', label = 'Trend Line')	
	ax1.legend(loc=1, prop={'size': 18})
	ax1.set_title('TIC {}'.format(dactyl.tic_id))
	ax1.set_ylabel('Flux', size = 18)
	if xmin != None: ax1.set_xlim(left = xmin)
	if xmax != None: ax1.set_xlim(right = xmax)
	if hasattr(dactyl, 'rot_rate'):
		text= '  Rotation Rate: ' + str(format(dactyl.rotation, '.2f')) + ' days'
		ax1.text(0.55, 0.05, text, fontweight="normal", size = 20, transform=ax1.transAxes)

	text2 = 'Sectors: ' + dactyl.sector_list
	ax1.text(0.01, 0.9, text2, fontweight="normal", size = 20, transform=ax1.transAxes)
	ax1.text(0.01, 0.05, "(A)", fontweight="normal", size = 20, transform=ax1.transAxes)
		
	ax2.plot(dactyl.lc.time_detrend, dactyl.lc.flat, '.g', label ='Detrended Flux')
	ax2.legend(loc=1, prop={'size': 18})
	ax2.set_xlabel('Time (in days)', size = 18)
	ax2.set_ylabel('Flux', size = 18)
	if xmin != None: ax2.set_xlim(left = xmin)
	if xmax != None: ax2.set_xlim(right = xmax)
	if ymin != None: ax2.set_ylim(bottom = ymin)
	if ymax != None: ax2.set_ylim(top = ymax)    
    
	ax2.text(0.01, 0.05, "(B)", fontweight="normal", size = 20, transform=ax2.transAxes)
	text3= 'RMS: ' + str(format(dactyl.lc.flat_rms, '.2f')) + 'ppm'
	ax2.text(0.01, 0.9, text3, fontweight="normal", size = 20, transform=ax2.transAxes)	
								  
	fig.tight_layout()
	plt.show()
	plt.close()

"""
Plots light curve with found transits for a given target, with the option to save plot to disk
"""	
def lightcurve_tce(dactyl, save_plot  = False):
	dactyl.sector_list = ','.join(str(v) for v in dactyl.sectors)
	fig, (ax1, ax2) = plt.subplots(2, 1, figsize = (10,7), sharex = True)
	plt.subplots_adjust(hspace = 0.01)
	
	ax1.scatter(dactyl.lc.time, dactyl.lc.flux_n, color = 'blue', marker = '.', alpha = 1, label = 'PSF Flux')
	ax1.plot(dactyl.lc.time, dactyl.lc.trend_n, color = 'red', linewidth = 2.0, linestyle = 'solid', label = 'Trend Line')	
	ax1.legend(loc=1, prop={'size': 18})
	ax1.set_title('TIC {}'.format(dactyl.tic_id))
	ax1.set_ylabel('Flux', size = 18)
	if hasattr(dactyl, 'rot_rate'):
		text= '  Rotation Rate: ' + str(format(dactyl.rotation, '.2f')) + ' days'
		ax1.text(0.55, 0.05, text, fontweight="normal", size = 20, transform=ax1.transAxes)

	text2 = 'Sectors: ' + dactyl.sector_list
	ax1.text(0.01, 0.9, text2, fontweight="normal", size = 20, transform=ax1.transAxes)
	ax1.text(0.01, 0.05, "(A)", fontweight="normal", size = 20, transform=ax1.transAxes)
		
	ax2.plot(dactyl.lc.time_detrend, dactyl.lc.flat, '.g', label ='Detrended Flux')
	for i in dactyl.planet['transit_times']:
		ax2.axvline(x=i, color='k', linestyle='--', alpha = 0.5)
	ax2.legend(loc=1, prop={'size': 18})
	ax2.set_xlabel('Time (in days)', size = 18)
	ax2.set_ylabel('Flux', size = 18)
	ax2.text(0.01, 0.05, "(B)", fontweight="normal", size = 20, transform=ax2.transAxes)
	text3= 'RMS: ' + str(format(dactyl.lc.flat_rms, '.2f')) + 'ppm'	
	text3+= '\nSDE: ' + str(format(dactyl.planet['SDE'], '.2f'))
	text3+= '\nSNR: ' + str(format(dactyl.planet['snr'], '.2f'))
	ax2.text(0.01, 0.7, text3, fontweight="normal", size = 20, transform=ax2.transAxes)
	fig.tight_layout()
	if save_plot:
		plt.savefig(dactyl.candidate_path+'lightcurve_tce', bbox_inches = 'tight')
	plt.show()
	plt.close()

"""
Same as lightcurve_tce but prettier
"""
def paper_plot(dactyl, xmin = None, xmax = None, ymin = None, ymax = None, loc = 0, bbox = [0,1], save_plot = True):
	fig, (ax1, ax2) = plt.subplots(2, 1, figsize = (10,7), sharex = True)    
	plt.subplots_adjust(hspace = 0.01)    
	l1 = ax1.scatter(dactyl.lc.time, dactyl.lc.flux_n, color = 'dimgrey', marker = '.', s = 100, alpha = 0.7, edgecolors = None, label = 'Normalized Flux')
	l2, = ax1.plot(dactyl.lc.time, dactyl.lc.trend_n, color = 'deepskyblue',  alpha=1.0, linewidth = 2.0, linestyle = 'solid', label = 'Trend Line')	
# 	ax1.set_title('TIC {}'.format(dactyl.tic_id), size = 18)
	ax1.set_ylabel('Flux', size = 18)


	ax2.scatter(dactyl.lc.time_detrend, dactyl.lc.flat, color = 'dimgrey', marker = '.', alpha = 0.75, edgecolors = None, label = 'Detrended Flux')
	l4, = ax2.plot(dactyl.planet['model_lightcurve_time'],dactyl.planet['model_lightcurve_model'], alpha=1, color='red', zorder=1, label = 'Transit Model')
	ax2.set_xlabel('TESS Julian Date (in days)', size = 18)
	ax2.set_ylabel('Flux', size = 18)
	if xmin != None: ax2.set_xlim(left = xmin)
	if xmax != None: ax2.set_xlim(right = xmax)    
	if ymin != None: ax2.set_ylim(bottom = ymin)
	if ymax != None: ax2.set_ylim(top = ymax)     
	plt.legend(handles=[l1, l2, l4], loc=loc, bbox_to_anchor=bbox,
           ncol=4, prop={'size': 20}, markerscale=2, shadow=True, fancybox=True)

	if save_plot:
		plt.savefig(dactyl.candidate_path+str(dactyl.tic_id)+'.pdf', bbox_inches = 'tight')

	plt.show()
	plt.close()

"""
Plots periodogram of the planet signal, with option to save plot to disk
"""    
def planet_periodogram(dactyl, save_plot = False):
    plt.figure()
    ax = plt.gca()
    ax.axvline(dactyl.planet['P'], alpha=0.4, lw=3)
    plt.xlim(np.min(dactyl.planet['periods']), np.max(dactyl.planet['periods']))
    for n in range(2, 10):
        ax.axvline(n*dactyl.planet['P'], alpha=0.4, lw=1, linestyle="dashed")
        ax.axvline(dactyl.planet['P'] / n, alpha=0.4, lw=1, linestyle="dashed")
    plt.ylabel(r'SDE')
    plt.xlabel('Period (days)')
    plt.plot(dactyl.planet['periods'], dactyl.planet['power'], color='black', lw=0.5)
    plt.title('Periodogram of TIC {}'.format(dactyl.tic_id))
    plt.xlim(0, max(dactyl.planet['periods']));
    if save_plot:
        plt.savefig(dactyl.candidate_path+'planet_periodogram', bbox_inches = 'tight')
    plt.show()
    plt.close()

"""
Plots zoomed in version of individual transits for a given target, with adjustable axes and option to save plot to disk
"""	    
def transit_zoom(dactyl, xmin = None, xmax = None, save_plot = False):
    limits = 4*dactyl.planet['tdur_obs']
    if len(dactyl.planet['transit_times']) % 3 == 0:
        nrows = len(dactyl.planet['transit_times']) // 3
    else:
        nrows = (len(dactyl.planet['transit_times']) // 3) + 1
    fig, axes = plt.subplots(nrows=nrows, ncols=3, figsize = (12, len(dactyl.planet['transit_times'])//3*4), constrained_layout=True)
    plt.subplots_adjust(hspace = 0.75, wspace = 0.5)
    
    st = plt.suptitle('TIC {}'.format(dactyl.tic_id), fontsize = 18)
    colors = iter(cm.plasma(np.linspace(0, 1, len(dactyl.planet['transit_times']))))
    axes = axes.ravel()
    in_transit = transit_mask(dactyl.lc.time_detrend,dactyl.planet['P'],dactyl.planet['tdur_obs'],dactyl.planet['T0'])
    for i in range(len(axes)):
        try:
            axes[i].scatter(dactyl.lc.time_detrend, dactyl.lc.flat, color=next(colors));
            axes[i].scatter(dactyl.lc.time_detrend[in_transit], dactyl.lc.flat[in_transit], color='red')    
            axes[i].plot(dactyl.planet['model_lightcurve_time'],dactyl.planet['model_lightcurve_model'], alpha=0.5, color='red', zorder=1)
            axes[i].set_xlim([dactyl.planet['transit_times'][i]-limits, dactyl.planet['transit_times'][i]+limits])
            axes[i].set_xlabel('TJD', fontsize=14)
            axes[i].set_ylabel('Flux', fontsize=14)  
        except StopIteration:
            pass
        
    if save_plot:
        plt.savefig(dactyl.candidate_path+'transit_zoom', bbox_inches = 'tight')
    plt.tight_layout()
    plt.show()
    plt.close()
 
"""
Plots phase-folded light curve for a given target, with adjustable axes and option to save plot to disk
"""   
def planet_phasefold(dactyl, xmin = 0.4, xmax = 0.6, save_plot = False):
    plt.plot(dactyl.planet['model_folded_phase'], dactyl.planet['model_folded_model'], color='red')
    plt.scatter(dactyl.planet['folded_phase'], dactyl.planet['folded_y'], color='blue', s=10, alpha=0.5, zorder=2)
    plt.xlabel('Phase')
    plt.ylabel('Relative flux')
    plt.title('Phase Folded Light Curve of TIC {}'.format(dactyl.tic_id))
    if xmin != None: plt.xlim(left = xmin)
    if xmax != None: plt.xlim(right = xmax)
    if save_plot:
        plt.savefig(dactyl.candidate_path+'planet_phasefold', bbox_inches = 'tight')
    plt.show()
    plt.close()

"""
Plot to depict how individual transit depths compare to the mean transit deapth for a given target, with option to save plot to disk
"""
def tdepth(dactyl, save_plot = False):	
	fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize = (10,10), sharex = True)
	plt.subplots_adjust(hspace = 0.1)
	ax1.set_title('TIC {}'.format(dactyl.tic_id))
	ax1.errorbar(dactyl.planet['transit_times'], dactyl.planet['transit_depths'], yerr = dactyl.planet['transit_depths_uncertainties'], fmt = 'o', color = 'red')
	ax1.axhline(y=dactyl.planet['depth'], color='k', linestyle='--', alpha = 0.5)
	depths_og = 1-dactyl.planet['transit_depths']
	depths = [value for value in depths_og if not math.isnan(value)]
	text1 = 'std/avg tdepth = ' + str(np.std(depths)/np.median(depths))
	dactyl.tdepth_metric = np.std(depths)/np.median(depths)
	ax1.text(0.01, 0.9, text1, fontweight="normal", size = 18, transform=ax1.transAxes)
	ax1.set_ylabel('1 - Transit Depths', size = 18)
	ax1.set_ylim(top = 1.0005, bottom = 0.992)
	
	ax2.plot(dactyl.planet['transit_times'], dactyl.planet['snr_per_transit'], 'o', color = 'orange')	
	ax2.axhline(y=dactyl.planet['snr'], color='blue', linestyle='--', alpha = 0.5)	
	ax2.set_ylabel('SNR', size = 18)
	
	ax3.plot(dactyl.planet['transit_times'], dactyl.planet['per_transit_count'], 'o', color = 'green')		
	ax3.set_ylabel('Data Points in Transit', size = 18)	
	ax3.set_xlabel('Time (in days)', size = 18)
	if save_plot:
		plt.savefig(dactyl.candidate_path+'tdepth', bbox_inches = 'tight')	
	plt.show()
	plt.close()

"""
Plots injected light curve
"""	
def lightcurve_injection(inject):	
	plt.plot(inject.lc.time, (inject.lc.injection-1)*1e6)
	plt.title('Injected planet')
	plt.ylabel('Signal [ppm]')
	plt.xlabel('Time')

	text= 'P = {:.1} days'.format(inject.P)
	text+= '\nR= {:.1} Rearth'.format(inject.Rp)
	plt.text(0.02, 0.97, text, ha='left', va='top', transform = plt.gca().transAxes)

	plt.show()
	plt.close()

"""
Plot the Injection-Recovery grid from recovery.py
"""
def injection_per_star(INJ, Save=False, Grid=False, Show=True, DoRpRs=None):	
	rec= INJ.rec
	fail= INJ.status=='Failed'
	FP= INJ.status=='Other Signal'

	if DoRpRs is None:
		PlotRpRs = INJ.NoRstar
		PlotGridAxes=True
	else:
		PlotRpRs= DoRpRs
		PlotGridAxes= (DoRpRs == INJ.NoRstar)
	#print (PlotGridAxes)

	plt.title('TIC {}'.format(INJ.ticID))
	plt.xlabel('Orbital Period [days]')
	if PlotRpRs:
		plt.ylabel('Planet Radius [R$_\\bigstar$]')
	else:
		plt.ylabel('Planet Radius [earth]')
	
	if PlotGridAxes:
		plt.xlim(INJ.Pbins[0]*0.9, INJ.Pbins[-1]*1.1)
		plt.ylim(INJ.Rpbins[0]*0.9, INJ.Rpbins[-1]*1.1)
	
	if Grid:
		for xP in INJ.Pbins:
			plt.vlines(xP, ymin= INJ.Rpbins[0], ymax= INJ.Rpbins[-1], 
						color='0.75')
		for yR in INJ.Rpbins:
			plt.hlines(yR, xmin= INJ.Pbins[0], xmax= INJ.Pbins[-1], 
						color='0.75')
	
	# plot Rp or Rp_Rstar?
	Rp_plot= INJ.rp_rs_injs if PlotRpRs else INJ.Rp_injs
	#print (Rp_plot)

	plt.loglog(INJ.P_injs, Rp_plot, 'y.', 
			   label='Injected')
	plt.loglog(INJ.P_injs[rec], Rp_plot[rec], 'x', color='purple', ls='', 
			   label='Recovered')
	plt.loglog(INJ.P_injs[fail], Rp_plot[fail], 'v', color='blue', ls='',
			  label='False Negative')
	plt.loglog(INJ.P_injs[FP], Rp_plot[FP], 'o', color='r', ls='',
			  label='False Positive')

	plt.legend(bbox_to_anchor=(1, 1),)
	
	if Save:
		os.makedirs('png/injections/{}'.format(INJ.name), exist_ok=True)
		if PlotRpRs:  fplot= 'png/injections/{}/{}-RpRs'.format(INJ.name, INJ.ticID)
		else:       fplot= 'png/injections/{}/{}'.format(INJ.name, INJ.ticID)
		plt.savefig(fplot, bbox_inches='tight',dpi=150)
		if Show: plt.show()
		plt.close()