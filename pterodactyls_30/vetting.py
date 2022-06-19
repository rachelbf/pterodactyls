import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import triceratops.triceratops as tr
from lightkurve import TessLightCurve
from tabulate import tabulate
import pandas as pd
from .plotter import *
import math
import os
import EDIunplugged as EDI
import lightkurve as lk
from vetting import centroid_test

"""
Vets the planet signal using tests inspired from Kepler's Robovetter and K2's EDI Vetter; See section 3.4 in pterodactyls paper for details
Args:
    tdur_cutoff (float): Cutoff value for the transit duration vetting test; default: 0.75
    tdepth_cutoff (float): Cutoff value for the transit depth vetting test; default: 1.5
    doplot (bool): Plots diagnostic plots; default: False   
    Verbose (bool): Prints out warnings and any other relevant statements; default: False
""" 
def vetting(dactyl, tdur_cutoff = 0.75, tdepth_cutoff = 1.5, doplot = False, Verbose = False):   
    if hasattr(dactyl, 'TCE'):
        if dactyl.TCE:         
            dactyl.vet_flags = []
            print ('Running Vetting routine...')
            print("==========================================") 
            depths_og = 1-dactyl.planet['transit_depths']
            depths = [value for value in depths_og if not math.isnan(value)]
            dactyl.tdepth_metric = np.std(depths)/np.median(depths)
            dactyl.tdur_metric = dactyl.planet['tdur_obs']/dactyl.planet['tdur_exp']
            
            if dactyl.planet['distinct_transit_count'] == 1:
                  dactyl.single_transit_flag = True
                  dactyl.vet_flags.append('single_transit_flag')                  
            else:     
                  dactyl.single_transit_flag = False            
                  
            if np.abs(1 - dactyl.planet['P']/dactyl.rotation) < 0.01:
                  dactyl.per_rot_flag = True
                  dactyl.vet_flags.append('per_rot_flag')
            else:
                  dactyl.per_rot_flag = False                
                  
            if np.abs(0.5 - dactyl.planet['P']/dactyl.rotation) < 0.01:
                  dactyl.per_half_rot_flag = True
                  dactyl.vet_flags.append('per_half_rot_flag')
            else:
                  dactyl.per_half_rot_flag = False  
                    
            if np.abs(dactyl.tdepth_metric) >= tdepth_cutoff:
                  dactyl.tdepth_flag = True
                  dactyl.vet_flags.append('tdepth_flag')
            else:
                  dactyl.tdepth_flag = False  
 
            
            if (dactyl.tdur_metric <= (1-tdur_cutoff)) | (dactyl.tdur_metric >= 1/(1-tdur_cutoff)):
                  dactyl.tdur_flag = True
                  dactyl.vet_flags.append('tdur_flag')
            else:
                  dactyl.tdur_flag = False  

            ''' Run Eddie Vetter '''
            dactyl.params=EDI.parameters(dactyl.tls_results, limbDark=dactyl.ab, impact=0, snrThreshold=7, minTransit=2)
            EDI.Go(dactyl.params, telescope='TESS', print=False)
            if Verbose:
                print('Single Transit: ', dactyl.single_transit_flag)
                print('Orbital Period Similar to Stellar Rotation Period: ', dactyl.per_rot_flag)
                print('Orbital Period Similar to half Stellar Rotation Period: ', dactyl.per_half_rot_flag)
                print('Inconsistent Transit Depths: ', dactyl.tdepth_flag)   
                print('Observed Transit Duration Inconsistent with Expected transit Duration: ', dactyl.tdur_flag)
                print('EDIVetter Individual Transit Test: ', dactyl.params.transMaskFP) 
                print('EDIVetter Odd/Even Transit Variation: ', dactyl.params.evenOddFP)
                print("==========================================") 

                
            if (len(dactyl.vet_flags)==0) & (dactyl.params.transMaskFP == False) & (dactyl.params.evenOddFP == False):
                dactyl.vetting_tests_passed = True
                if Verbose:
                    print ('\nTCE passed all vetting tests! \nPoI found in TIC {}'.format(dactyl.tic_id))
                    print ('  SDE = {:.1f}'.format(dactyl.planet['SDE']) )
                    print ('  SNR = {:.1f}'.format(dactyl.planet['snr']) )
                    print ('  Period  = {:.2f} days'.format(dactyl.planet['P']) )
                    print ('  Period error = {:.5f} days'.format(dactyl.planet['P_err']) )
                    print ('  Radius = {:.2f} earth radii'.format(dactyl.planet['rp_rs'] * 110.*dactyl.radius) ) 
                    print ('  T0 = {:.2f}'.format(dactyl.planet['T0']) )
                    print(len(dactyl.planet['transit_times']), 'transit times in time series:', \
                    ['{0:0.5f}'.format(i) for i in dactyl.planet['transit_times']])
                    print('Number of data points during each unique transit', dactyl.planet['per_transit_count'])
                    print('The number of transits with intransit data points', dactyl.planet['distinct_transit_count'])
                    print('The number of transits with no intransit data points', dactyl.planet['empty_transit_count'])
                    print('Transit depth:', format(dactyl.planet['depth'], '.5f'))
                    print('Transit duration:', format(dactyl.planet['tdur_obs'], '.5f'), 'days')
                    print('Transit depths (mean)', dactyl.planet['transit_depths'])
                    print('Transit depth uncertainties', dactyl.planet['transit_depths_uncertainties'])
                    print('Radius ratio of planet and star:', format(dactyl.planet['rp_rs'], '.5f'))
                    print('Stellar Radius:', format(dactyl.radius, '.5f'), 'Rsun')
                
                    if doplot: 
                        lightcurve_tce(dactyl, save_plot=True)
                        planet_periodogram(dactyl, save_plot=True)
                        planet_phasefold(dactyl, save_plot=True)
                        tdepth(dactyl, save_plot=True)
                    print ('*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*')
            else:
                dactyl.vetting_tests_passed=False
                if dactyl.params.transMaskFP:
                    dactyl.vet_flags.append('EDI: transMaskFP')
                if dactyl.params.evenOddFP:
                    dactyl.vet_flags.append('EDI: evenOddFP')
        else:
            dactyl.vetting_tests_passed=False

"""
Calculates flux contribution of a given target to a 5x5 pixel aperture using triceratops
Args:
    doplot (bool): Plots diagnostic plots; default: False   
    Verbose (bool): Prints out warnings and any other relevant statements; default: False
    save_on_disk (bool): Saves flux ratio and top 5 flux contributors to disk for optimization; default: False 
    Overwrite (bool): Overwrites the flux ratio and top 5 flux contributors that were saved to disk; default: False
""" 
def flux_contamination(dactyl, doPlot = True, Verbose = False, save_on_disk = True, Overwrite= False):      
    if save_on_disk:
        if not os.path.exists(dactyl.fratio_path): os.makedirs(dactyl.fratio_path)
        flux_file = dactyl.fratio_path+'fratio.npz'
    if save_on_disk and os.path.isfile(flux_file) and not Overwrite:
        print('(load {})'.format(flux_file))
        npz = np.load(flux_file)
        dactyl.fluxratio = npz['fratio'][()] # variable not array
        npz.close()
    else:
        
        if save_on_disk:
            print('Calculating Flux Ratio of TIC {}'.format(dactyl.tic_id))
            sectors = np.array(dactyl.sectors)
            target = tr.target(ID = dactyl.tic_id, sectors = sectors, search_radius=6.5)
                            
            for i,sector in enumerate(sectors):
                print("Sector: ", sector)
                if doPlot:
                    target.plot_field(sector=sector, ap_pixels=None)                    
            target.calc_depths(tdepth=0, all_ap_pixels=None)            
            if Verbose:
                print(tabulate(target.stars, headers="keys", floatfmt=".4f", numalign="left", tablefmt="simple"))
            
            dactyl.fluxratio = target.stars.iloc[0]['fluxratio']

            # save 5 brightest stars if TCE
            if os.path.isdir(dactyl.candidate_path):
                contam_file = dactyl.candidate_path+'contamination.txt'              
                contam_csv = target.stars.nlargest(5, ['fluxratio'])
                contam_csv.to_csv(contam_file, sep = ',', index = False)
          
        #    print(target.stars.iloc[0])    
            np.savez_compressed(flux_file, fratio = dactyl.fluxratio)
            print ('(saved {})'.format(flux_file))

"""
Vets phase folded light curve using triceratops: https://github.com/stevengiacalone/triceratops
Args:
    phase_min (float): how much pre-transit baseline to consider for triceratops analysis; default: 0.4
    phase_max (float): how much post-transit baseline to consider for triceratops analysis; default: 0.6
    flatpriors (bool): When True, turns off Gyr-old occurrence rates priors used by triceratops; default: False 
    doPlot (bool): Plots diagnostic plots; default: False   
    Verbose (bool): Prints out warnings and any other relevant statements; default: False
    save_on_disk (bool): Saves flux ratio and top 5 flux contributors to disk for optimization; default: False 
"""         
def vet_phase_fold_lc(dactyl, phase_min = 0.4, phase_max = 0.6, flatpriors = False, Verbose = True, doPlot = True, save_on_disk = True):     
    if hasattr(dactyl, 'TCE'):
        if dactyl.TCE:
            sectors = np.array(dactyl.sectors)
            target = tr.target(ID = dactyl.tic_id, sectors = sectors, search_radius=6.5)
            
            ''' Alter phase folded light curve for better looking comparison with triceratops'''
            ph_lc = pd.DataFrame(
                {'folded_phase': dactyl.planet['folded_phase'],
                  'folded_y': dactyl.planet['folded_y'],
                  'folded_dy': dactyl.planet['folded_dy']
                })
            mask = (ph_lc['folded_phase'] >= phase_min) & (ph_lc['folded_phase'] <=phase_max)
            folded_phase = ph_lc['folded_phase'][mask] - 0.5 #shift from TLS axes to triceratops axes
            folded_y = ph_lc['folded_y'][mask]
            folded_dy = ph_lc['folded_dy'][mask]
            print('Vetting the phase-folded light curve of TIC {}'.format(dactyl.tic_id))
            
            target.calc_depths(tdepth=1-dactyl.planet['depth'], all_ap_pixels=None)            
                
            target.stars['mass'] = target.stars['mass'].replace(np.nan, 1.0)
            target.stars['rad'] = target.stars['rad'].replace(np.nan, 1.0)
            target.stars['Teff'] = target.stars['Teff'].replace(np.nan, 5777.0)
            
            if dactyl.cluster == '32Ori':
                target.stars['plx'] = target.stars['plx'].replace(np.nan, 10.41666667)
            elif dactyl.cluster == 'ABDMG':
                target.stars['plx'] = target.stars['plx'].replace(np.nan, 33.33333333)
            elif dactyl.cluster == 'AlphaPer':
                target.stars['plx'] = target.stars['plx'].replace(np.nan, 5)                
            elif dactyl.cluster == 'Argus':
                target.stars['plx'] = target.stars['plx'].replace(np.nan, 8.333333333)
            elif dactyl.cluster == 'BetaPic':
                target.stars['plx'] = target.stars['plx'].replace(np.nan, 33.33333333)                
            elif dactyl.cluster == 'Blanco1':
                target.stars['plx'] = target.stars['plx'].replace(np.nan, 3.95256917)
            elif dactyl.cluster == 'Carina':
                target.stars['plx'] = target.stars['plx'].replace(np.nan, 16.66666667)
            elif dactyl.cluster == 'CarinaNear':
                target.stars['plx'] = target.stars['plx'].replace(np.nan, 33.33333333)                
            elif dactyl.cluster == 'Columba':
                target.stars['plx'] = target.stars['plx'].replace(np.nan, 20)                
            elif dactyl.cluster == 'ComaBer':
                target.stars['plx'] = target.stars['plx'].replace(np.nan, 11.76470588)                
            elif dactyl.cluster == 'EtaCha':
                target.stars['plx'] = target.stars['plx'].replace(np.nan, 10.52631579)                
            elif dactyl.cluster == 'Hyades':
                target.stars['plx'] = target.stars['plx'].replace(np.nan, 22.22222222)                
            elif dactyl.cluster == 'IC2391':
                target.stars['plx'] = target.stars['plx'].replace(np.nan, 6.711409396) 
            elif dactyl.cluster == 'IC2602':
                target.stars['plx'] = target.stars['plx'].replace(np.nan, 6.711409396)  
            elif dactyl.cluster == 'LCC':
                target.stars['plx'] = target.stars['plx'].replace(np.nan, 9.090909091)  
            elif dactyl.cluster == 'MUTA':
                target.stars['plx'] = target.stars['plx'].replace(np.nan, 6.666666667)  
            elif dactyl.cluster == 'NGC2451':
                target.stars['plx'] = target.stars['plx'].replace(np.nan, 3.703703704)  
            elif dactyl.cluster == 'Octans':
                target.stars['plx'] = target.stars['plx'].replace(np.nan, 7.692307692)  
            elif dactyl.cluster == 'PiEri':
                target.stars['plx'] = target.stars['plx'].replace(np.nan, 6.535947712)  
            elif dactyl.cluster == 'Platais':
                target.stars['plx'] = target.stars['plx'].replace(np.nan, 7.692307692)  
            elif dactyl.cluster == 'Pleiades':
                target.stars['plx'] = target.stars['plx'].replace(np.nan, 7.462686567)  
            elif dactyl.cluster == 'THA':
                target.stars['plx'] = target.stars['plx'].replace(np.nan, 21.73913043)                 
            elif dactyl.cluster == 'TWHya':
                target.stars['plx'] = target.stars['plx'].replace(np.nan, 16.66666667)                 
            elif dactyl.cluster == 'UCL':
                target.stars['plx'] = target.stars['plx'].replace(np.nan, 7.692307692)                 
            elif dactyl.cluster == 'UMa':
                target.stars['plx'] = target.stars['plx'].replace(np.nan, 40)                 
            elif dactyl.cluster == 'UCRA':
                target.stars['plx'] = target.stars['plx'].replace(np.nan, 6.802721088)                 
            elif dactyl.cluster == 'UppSco':
                target.stars['plx'] = target.stars['plx'].replace(np.nan, 7.692307692)                 
            elif dactyl.cluster == 'VCA':
                target.stars['plx'] = target.stars['plx'].replace(np.nan, 11.76470588) 
            elif dactyl.cluster == 'XFOR':
                target.stars['plx'] = target.stars['plx'].replace(np.nan, 10)                 
               
                
            target.stars['mass'][target.stars['ID'] == '466035000'] = 3.1788               
            target.stars['rad'][target.stars['ID'] == '466035000'] = 6.47  
            target.stars['Teff'][target.stars['ID'] == '466035000'] = 6000             
            target.stars['plx'][target.stars['ID'] == '466035000'] = 0.00119721623
                
            #vetting phase folded light curve
            time, flux, flux_err = folded_phase*dactyl.planet['P'], folded_y, folded_dy
            P_orb = dactyl.planet['P']
            
            lc_binsize = (time.max()-time.min())/100
            lc = TessLightCurve(time=time, flux=flux, flux_err=flux_err)

            target.calc_probs(time=lc.time.value, flux_0=lc.flux.value, flux_err_0=np.mean(lc.flux_err.value), P_orb=P_orb, parallel=True, flatpriors = flatpriors)
            # target.calc_probs(time=lc.time.value, flux_0=lc.flux.value, flux_err_0=np.mean(lc.flux_err.value), P_orb=P_orb, drop_scenario = ["PTP", "PEB", "STP", "SEB", "DTP", "DEB", "BTP", "BEB"], parallel=True, flatpriors = flatpriors)
            df_results = target.probs
            if flatpriors:
                prob_file = dactyl.candidate_path+'triceratops_results_flat.txt'
                fname = dactyl.candidate_path+'vet_lc_flat'
            else:
                prob_file = dactyl.candidate_path+'triceratops_results_prior.txt'
                fname = dactyl.candidate_path+'vet_lc_prior'
            if Verbose:
                print("False Positive Probability: ", np.round(target.FPP, 4))
                print("Nearby False Positive Probability: ", np.round(target.NFPP, 4))
                print(tabulate(df_results, headers="keys", floatfmt=".4f", numalign="left", tablefmt="simple"))
                df_results.to_csv(prob_file, sep = ',', index = False)
            if doPlot:
                target.plot_fits(time=lc.time.value, flux_0=lc.flux.value, flux_err_0=np.mean(lc.flux_err.value), save = True, fname = fname)
            
"""
Does centroid test to determine if transit signal is indeed coming from a given target, from Christina Hedges' vetting package: https://iopscience.iop.org/article/10.3847/2515-5172/ac376a
Args:
    periods (float): orbital period of the transit signal
    T0s (float): T0s of the transit signal
    tdurs (float): transit duration of the transit signal
    aperture_mask (str): An aperture mask type to perform the test on. See `lk.TargetPixelFile._parse_aperture_mask`. 'all' will use all pixels,`pipeline` will use the TESS pipeline mask, and `threshold` 
                    will use all contiguous pixels above a 3 sigma threshold; default: 'pipeline'
    doPlot (bool): Plots centroid plot; default: False   
    savePlot (bool): Saves plot to disk for optimization; default: True
"""  
def TESS_centroid_test(dactyl, periods, T0s, tdurs, aperture_mask = 'pipeline', doPlot = True, savePlot = True): 
    for s in dactyl.sectors:
        tpf = lk.search_targetpixelfile('TIC'+str(dactyl.tic_id), mission='TESS', sector=s).download()
        period = periods
        T0 = T0s
        dur = tdurs
        r = centroid_test(tpf, period, T0, dur, aperture_mask=aperture_mask, plot=doPlot)
        if savePlot:
            plt.savefig(dactyl.candidate_path+'TIC{}_s{}_centroidtest.pdf'.format(dactyl.tic_id,s), bbox_inches = 'tight')