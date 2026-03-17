#!/bin/python

import numpy as np
from matplotlib import pyplot as plt
import matplotlib.ticker
from matplotlib import rcParams
import matplotlib.ticker as tck
from matplotlib.patches import Circle
from matplotlib.patches import Ellipse
from matplotlib.patches import FancyArrowPatch
from matplotlib.ticker import FuncFormatter
from astropy.io import fits
from astropy import wcs
from astropy.coordinates import SkyCoord, Longitude, Latitude, Angle
from astropy.visualization import ZScaleInterval
import astropy.units as u
from ztfquery import lightcurve
import requests
import swifttools.ukssdc.data.SXPS as uds
import mastcasjobs
from astroquery.simbad import Simbad
from astroquery.hips2fits import hips2fits
from astropy.wcs.utils import proj_plane_pixel_scales
import io
import os
import shutil
import time
from scipy.stats import binned_statistic

from periodicity_search import LombScargle_search, LombScargle_FAP, LombScargle_FAP_level, PDM_search, PDM_FAP, PDM_FAP_level, PDM_bisearch, PDM_trisearch

path=os.getcwd()

#We get as reference epoch for period folding the time (in MJD) of the first image in the gband
#t_0 = float(np.genfromtxt(path+'/gfilter/usedImages.txt',dtype='str')[0].split('_')[2])

analysis = input("Enter 'STELLA', 'LCO' or 'ZTF' for the type of data (default value 'LCO'):\n") or 'LCO'
while analysis!='STELLA' and analysis!='LCO' and analysis!='ZTF':
   analysis = input("Invalid value, please insert 'STELLA' or 'LCO' or 'ZTF':\n") or 'LCO'
if analysis == 'ZTF':
   FoVname = input("Enter the name of the target (4FGL source name without spaces):\n")
lc_plot = input("Enter 'yes' if you want to produce light curves plots (default value 'yes'):\n")  or "yes"
while lc_plot!='yes' and lc_plot!='no':
   lc_plot = input("Invalid value, please insert 'yes' or 'no':\n") or "yes"
if lc_plot == 'yes':
   ncand = int(input("Enter the number corresponding to the periodic source for which you want to publish plots:\n"))

   per_search_sub = input("Enter 'LS1', 'PDM', 'multiPDM' or 'LS2' for the periodicity search method (default value 'PDM'):\n") or 'PDM'
   while per_search_sub!='LS1' and per_search_sub!='PDM' and per_search_sub!='multiPDM' and per_search_sub!='LS2':
      per_search_sub = input("Invalid value, please insert 'LS1', 'PDM', 'multiPDM' or 'LS2':\n") or 'PDM'
   n_periods_sub = float(input("Enter an integer positive number for number of period steps (10-10000) to use (default value 2000):\n")  or "2000")
   while n_periods_sub<10 or n_periods_sub>10000.:
      n_periods_sub = float(input("Invalid value, please insert a integer number between 10 and 10000:\n") or "2000")
   n_periods_sub = int(n_periods_sub)
   P_low_sub  = float(input("Enter the start period (in days) used (default value 0.02):\n") or 0.02)
   while P_low_sub<=0.:
      P_low_sub = float(input("Invalid value, please insert a positive start period:\n") or 0.02)
   P_up_sub  = float(input("Enter the final period (in days) used (default value 2.5):\n") or 2.5)
   while P_up_sub<=0.:
      P_up_sub = float(input("Invalid value, please insert a positive final period:\n") or 2.5)
   if per_search_sub == 'PDM' or per_search_sub == 'multiPDM':
      n_bin_sub = float(input("Enter an integer positive number for number of phase bins to use in the PDM (or multiband-PDM) method (default value 10):\n")  or "10")
      while n_bin_sub<0. or n_bin_sub%1!=0.:
         n_bin_sub = float(input("Invalid value, please insert a positive integer number (default value 10):\n") or "10")
      n_bin_sub = int(n_bin_sub)
   period_error = input("Enter 'yes' if you want to estimate the period error (default value 'no', it can strongly increase the computation times):\n")  or "no"
   while period_error!='yes' and period_error!='no':
      period_error = input("Invalid value, please insert 'yes' or 'no':\n") or "no"
   if period_error == 'yes':
      n_steps_sub = float(input("Enter an integer positive number for number of steps to use in the simulation for the error estimation on the period (default value 1000):\n")  or "1000")
      while n_steps_sub<10 or n_steps_sub>10000.:
         n_steps_sub = float(input("Invalid value, please insert a integer number between 10 and 10000:\n") or "1000")
      n_steps_sub = int(n_steps_sub)
   double_period = input("Enter 'yes' if you want to fold light curves at P=2*P_best (assumes ellipsoidal modulation, default value 'no'):\n")  or "no"
   while double_period!='yes' and double_period!='no':
      double_period = input("Invalid value, please insert 'yes' or 'no':\n") or "no"
   PS_plot = input("Enter 'yes' if you want to plot also power spectrum (default value 'no'):\n")  or "no"
   while PS_plot!='yes' and PS_plot!='no':
      PS_plot = input("Invalid value, please insert 'yes' or 'no':\n") or "no"
   if PS_plot == 'yes':
      period_FAPl = input("Enter 'yes' if you want to estimate the False-Alarm-Probability level at 0.1%, better than 3sigma (default value 'no', it can increase the computation times):\n")  or "no"
      while period_FAPl!='yes' and period_FAPl!='no':
         period_FAPl = input("Invalid value, please insert 'yes' or 'no':\n") or "no"
      if period_FAPl == 'yes':
         #Probability of periodic signal in the data
         conf_level = 99.9
         FAPl_factor = float(input("Enter an integer positive number for the reduction factor of the period grid (default value 1):\n")  or "1")
         while FAPl_factor<1 or FAPl_factor>100:
            FAPl_factor = float(input("Invalid value, please insert a integer number between 1 and 100:\n") or "1")
         FAPl_factor = int(FAPl_factor)

   if analysis == 'ZTF':
      lc_rebin = input("Enter 'yes' if you want to phase-rebin ZTF data (default value 'no'):\n")  or "no"
      while lc_rebin!='yes' and lc_rebin!='no':
         lc_rebin = input("Invalid value, please insert 'yes' or 'no':\n") or "no"
      if lc_rebin == 'yes':
         rebin_factor = int(input("Enter the number of phase bins that you want to use on ZTF data (default value 10):\n")  or "10")
         while rebin_factor<1:
            rebin_factor = int(input("Invalid value, please insert a number major than 1:\n") or "10")

#Built-in,minimum number of data points per phase bin (if lower bin is not used in PDM calculation)
min_bin = 2

#Defining scale factor (in pixel) to use as squared region around the target, with side=2*bound
bound = 1.05
#Defining radius size (in arcmin) to use for plotting 4FGL sources nearby to our candidate target
radius_gamma = 11
#Defining radius size (in arcsec) to consider association of an optical candidate with an X-ray counterpart
radius_X = 6.5
#Color range factor to be multiplied to data stdev in pyplot.imshow
colscale = 0.33
#colscale = 0.12
#Time spacing factor for colour light curve
t_space = 2.0
#Space between axes and labels
labelpad = 1.7
markersize=4
markersize_main=5
alpha=1.0
tsize = 18
subsize = 12
size = 15
axlabelsize = 16
axticksize = 16
#figxsize1 = 25
#figysize1 = 35
#figxsize2 = 30
#figysize2 = 42
figxsize = 7
figysize = 7
#Title y position
title_height = 0.95

if not os.path.isdir(path+'/finalperiodplots_paper'):
   os.system('mkdir '+path+'/finalperiodplots_paper')
   
if analysis == 'LCO' or analysis == 'STELLA':
   varSources_fin1=np.genfromtxt(path+'/rfilter/varSources_1.csv', delimiter=',')
   if os.path.exists(path+'/final_periodic_can_1.txt'):
      list_targ = np.int_(np.genfromtxt(path+'/final_periodic_can_1.txt'))
      nvar=np.size(list_targ)

      if nvar!=1:
         perSources_fin = varSources_fin1[list_targ-1]
      else:
         perSources_fin = varSources_fin1[list_targ-1,:]
else:
   varSources_fin1=np.genfromtxt(path+'/ZTF_search/ZTF_coord.dat')
   if os.path.exists(path+'/final_periodic_can_ZTF.txt'):
      list_targ = np.int_(np.genfromtxt(path+'/final_periodic_can_ZTF.txt'))
      nvar=np.size(list_targ)
      if nvar!=1:
         perSources_fin = varSources_fin1[list_targ-1]
      else:
         perSources_fin = varSources_fin1[list_targ-1,:]

#Load eventual X-ray counterparts for the Field of View from Chandra (with csc2.py script), XMM-Newton (with XMM4d13s.py script) and Swift (with swift2SXPS.py script) and radio counterparts from Karri's table if we are dealing with STELLA 2015 dataset
#if os.path.exists(path+'/2sxps.cat'):
#   os.remove(path+'/2sxps.cat')
#if os.path.exists(path+'/2sxps.dat'):
#   os.remove(path+'/2sxps.dat')
#if os.path.exists(path+'/2sxps.html'):
#   os.remove(path+'/2sxps.html')
#if os.path.exists(path+'/2sxps.log'):
#   os.remove(path+'/2sxps.log')
#if os.path.exists(path+'/2sxps.reg'):
#   os.remove(path+'/2sxps.reg')
#if os.path.exists(path+'/2sxps.tex'):
#   os.remove(path+'/2sxps.tex')
#if os.path.exists(path+'/csc2.cat'):
#   os.remove(path+'/csc2.cat')
#if os.path.exists(path+'/csc2.reg'):
#   os.remove(path+'/csc2.reg')
#if os.path.exists(path+'/XMM4d13s.cat'):
#   os.remove(path+'/XMM4d13s.cat')
#if os.path.exists(path+'/XMM4d13s.reg'):
#   os.remove(path+'/XMM4d13s.reg')
if os.path.exists(path+'/median_r_after.fits'):
   images = np.genfromtxt(path+'/rfilter.txt',dtype='str')
   #hdul = fits.open(path+'/rfilter/lsc1m004-fa03-20220609-0155-e91.fits')
   hdul = fits.open(path+'/median_r_after.fits')
   #hdul = fits.open(path+'/rfilter/median_r.fits')
   data = hdul[0].data
   hdr = hdul[0].header
   w = wcs.WCS(hdr)
   platescale = float(hdr['PIXSCALE'])
   hdul.close()
   if analysis=='STELLA':
      FoVname = hdr['OBJNAME']
      deltat = float(hdr['EXPT'])
      telescope = 'STELLA/WiFSIP'
   elif analysis=='LCO':
      FoVname = hdr['OBJECT']
      deltat = float(hdr['EXPTIME'])
      telescope = 'LCO/Sinistro'
   else:
      telescope = 'P48/ZTF'

   #FoVname = '4FGLJ2117.9+3729'
   #platescale = 0.333

#Load the corresponding 4FGL Fermi coordinates, 2sigma ellipses (in degrees) and 1 sigma semi-major axis (in degrees) of the target
#hdul = fits.open('/export/work/marcotu/gll_psc_v16.fit')
if FoVname[4] == '_':
   FoVname = FoVname.split('_')[0]+FoVname.split('_')[1]
hdul = fits.open('/export/work/marcotu/gll_psc_v28.fit')
[[ratar, dectar, ell_a_68, ell_a_95, ell_b_95, ell_theta_95]] = [[hdul[1].data['RAJ2000'][i],hdul[1].data['DEJ2000'][i],hdul[1].data['Conf_68_SemiMajor'][i],hdul[1].data['Conf_95_SemiMajor'][i],hdul[1].data['Conf_95_SemiMinor'][i],hdul[1].data['Conf_95_PosAng'][i]] for i in range(0,len(hdul[1].data)) if FoVname.replace("L","L ")==hdul[1].data['Source_Name'][i]]
hdul.close()
if not os.path.exists(path+'/median_r_after.fits'):
   #Query HIPS fits image of field from PanSTARRS around the target coordinates
   hdul = hips2fits.query(hips='CDS/P/PanSTARRS/DR1/r', width=300, height=300, projection='AZP', fov=(3 * ell_a_95* u.deg), ra=Longitude(ratar * u.deg), dec=Latitude(dectar * u.deg), format='fits', min_cut=0.5, max_cut=99.5)
   data = hdul[0].data
   hdr = hdul[0].header
   w = wcs.WCS(hdr)
   pixel_scales = proj_plane_pixel_scales(w)  # Returns scales in degrees per pixel
   platescale = 3600*pixel_scales[0]  # Assuming square pixels, take first value
   hdul.close()
   m, s = np.mean(data), np.std(data)

xtar, ytar = w.all_world2pix(ratar, dectar, 1)
ell_a_95 = ell_a_95*3600/platescale
ell_b_95 = ell_b_95*3600/platescale

#os.system('swift2SXPS.py '+str(ratar)+' '+str(dectar)+' '+str(4*ell_a_68*60.)+' -1 0')
if os.path.exists(path+'/2sxps.cat'):
   ra_Xcparts_Swift, dec_Xcparts_Swift, err_Xcparts_Swift = np.genfromtxt(path+'/2sxps.cat',usecols=(2,3,4),unpack=True)
   x_Xcparts_Swift, y_Xcparts_Swift = w.all_world2pix(ra_Xcparts_Swift, dec_Xcparts_Swift, 1)
   err_Xcparts_Swift = err_Xcparts_Swift/float(platescale)

#os.system('csc2.py '+str(ratar)+' '+str(dectar)+' '+str(4*ell_a_68*60.))
if os.path.exists(path+'/csc2.cat'):
   ra_Xcparts_Chandra, dec_Xcparts_Chandra, ell_a_Chandra, ell_b_Chandra, ell_theta_Chandra = np.genfromtxt(path+'/csc2.cat',usecols=(1,2,3,4,5),unpack=True)
   x_Xcparts_Chandra, y_Xcparts_Chandra = w.all_world2pix(ra_Xcparts_Chandra, dec_Xcparts_Chandra, 1)
   ell_a_Chandra = ell_a_Chandra/float(platescale)
   ell_b_Chandra = ell_b_Chandra/float(platescale)

#os.system('XMM4d13s.py '+str(ratar)+' '+str(dectar)+' '+str(4*ell_a_68*60.))
if os.path.exists(path+'/XMM4d13s.cat'):
   ra_Xcparts_XMM, dec_Xcparts_XMM, err_Xcparts_XMM = np.genfromtxt(path+'/XMM4d13s.cat',usecols=(1,2,3),unpack=True)
   x_Xcparts_XMM, y_Xcparts_XMM = w.all_world2pix(ra_Xcparts_XMM, dec_Xcparts_XMM, 1)
   err_Xcparts_XMM = err_Xcparts_XMM/float(platescale)

if os.path.exists(path+'/eRASS1d1.cat'):
   ra_Xcparts_eRASS, dec_Xcparts_eRASS, err_Xcparts_eRASS = np.genfromtxt(path+'/eRASS1d1.cat',usecols=(1,2,3),unpack=True)
   x_Xcparts_eRASS, y_Xcparts_eRASS = w.all_world2pix(ra_Xcparts_eRASS, dec_Xcparts_eRASS, 1)
   err_Xcparts_eRASS = err_Xcparts_eRASS/float(platescale)

data_crop = data[max(0,round(xtar-bound*max(ell_a_95,ell_b_95))):min(np.shape(data)[1]-1,round(xtar+bound*max(ell_a_95,ell_b_95))),max(0,round(ytar-bound*max(ell_a_95,ell_b_95))):min(np.shape(data)[0]-1,round(ytar+bound*max(ell_a_95,ell_b_95)))]

m, s = np.mean(data_crop), np.std(data_crop)

if lc_plot == 'yes':
   #Produce multi-band light curves plot and the corresponding field of view with optical and X-ray match (if any) for the periodic #ncand
   if nvar!=1:
      ra=perSources_fin[ncand-1,0]
      dec=perSources_fin[ncand-1,1]
      index=int(list_targ[ncand-1])
   else:
      ra=perSources_fin[0]
      dec=perSources_fin[1]
      index=int(list_targ)

   if period_error == 'yes':
      best_period_temp = np.zeros(n_steps_sub)
      
   if analysis == 'LCO' or analysis == 'STELLA':
      if os.path.exists(path+'/gfilter/outputcats/doerPhot_V'+str(index)+'.csv') and os.path.exists(path+'/rfilter/outputcats/doerPhot_V'+str(index)+'.csv') and os.path.exists(path+'/ifilter/outputcats/doerPhot_V'+str(index)+'.csv'):
         mjd_g, m_g, sigmam_g = np.genfromtxt(path+'/gfilter/outputcats/doerPhot_V'+str(index)+'.csv', delimiter=',', usecols=(8,-2,-1),unpack=True)
         #mjd_g, m_g, sigmam_g = np.genfromtxt(path+'/ZTF_g_J2117B.dat', usecols=(0,1,2),unpack=True)
         mjd_r, m_r, sigmam_r = np.genfromtxt(path+'/rfilter/outputcats/doerPhot_V'+str(index)+'.csv', delimiter=',', usecols=(8,-2,-1),unpack=True)
         #mjd_r, m_r, sigmam_r = np.genfromtxt(path+'/ZTF_r_J2117B.dat', usecols=(0,1,2),unpack=True)
         #ztf_mjd_r, ztf_m_r, ztf_sigmam_r = np.genfromtxt(path+'/ZTF_search/ZTF_r_23.dat', usecols=(0,1,2),unpack=True)
         mjd_i, m_i, sigmam_i = np.genfromtxt(path+'/ifilter/outputcats/doerPhot_V'+str(index)+'.csv', delimiter=',', usecols=(8,-2,-1),unpack=True)
   else:
      if os.path.exists(path+'/ZTF_search/ZTF_g_'+str(index)+'.dat'):
         mjd_g, m_g, sigmam_g = np.genfromtxt(path+'/ZTF_search/ZTF_g_'+str(index)+'.dat', usecols=(0,1,2),unpack=True)
      if os.path.exists(path+'/ZTF_search/ZTF_r_'+str(index)+'.dat'):
         mjd_r, m_r, sigmam_r = np.genfromtxt(path+'/ZTF_search/ZTF_r_'+str(index)+'.dat', usecols=(0,1,2),unpack=True)
      #m_g = m_g[mjd_g<57343]
      #sigmam_g = sigmam_g[mjd_g<57343]
      #mjd_g = mjd_g[mjd_g<57343]

      #m_r = m_r[mjd_r<57343]
      #sigmam_r = sigmam_r[mjd_r<57343]
      #mjd_r = mjd_r[mjd_r<57343]

      #m_i = m_i[mjd_i<57343]
      #sigmam_i = sigmam_i[mjd_i<57343]
      #mjd_i = mjd_i[mjd_i<57343]
      
      #m_g = np.array([m_g[j] for j in range(0,len(mjd_g)) if mjd_g[j] >= 58646.6170972 and mjd_g[j] <= 59062.74206480])
      #sigmam_g = np.array([sigmam_g[j] for j in range(0,len(mjd_g)) if mjd_g[j] >= 58646.6170972 and mjd_g[j] <= 59062.74206480])
      #mjd_g = np.array([mjd_g[j] for j in range(0,len(mjd_g)) if mjd_g[j] >= 58646.6170972 and mjd_g[j] <= 59062.74206480])
      #m_r = np.array([m_r[j] for j in range(0,len(mjd_r)) if mjd_r[j] >= 58632.8208449 and mjd_r[j] <= 59047.1511574])
      #sigmam_r = np.array([sigmam_r[j] for j in range(0,len(mjd_r)) if mjd_r[j] >= 58632.8208449 and mjd_r[j] <= 59047.1511574])
      #mjd_r = np.array([mjd_r[j] for j in range(0,len(mjd_r)) if mjd_r[j] >= 58632.8208449 and mjd_r[j] <= 59047.1511574])

   if per_search_sub == 'LS1':
      periods, periodogram = LombScargle_search(mjd_r, m_r, sigmam_r, P_low_sub, P_up_sub, n_periods_sub, 1)
      best_period = periods[~np.isnan(periodogram)][np.argmax(periodogram[~np.isnan(periodogram)])]
      #best_period = periods[~np.isnan(periodogram)][np.argsort(periodogram[~np.isnan(periodogram)])][-2]
      #best_period = 0.3054

      if period_error == 'yes':
         #Light curve simulation following Boldin et al. (2012) with a gaussian distribution (m'_i=np.random.normal(m_i,sigmam_i)) to estimate the error on the period
         for i in range(0,n_steps_sub):
            #m_r_temp = m_r + sigmam_r*np.random.uniform(-1.,1.,np.size(m_r))
            m_r_temp = [np.random.normal(m_r[j],sigmam_r[j]) for j in range(0,np.size(m_r))]
            m_r_temp = np.array(m_r_temp)
            periods_temp, periodogram_temp = LombScargle_search(mjd_r, m_r_temp, sigmam_r, P_low_sub, P_up_sub, n_periods_sub, 1)
            best_period_temp[i] = periods_temp[~np.isnan(periodogram_temp)][np.argmax(periodogram_temp[~np.isnan(periodogram_temp)])]
         sigma_best_period = np.std(best_period_temp,ddof=1)
   elif per_search_sub == 'PDM':
      periods, periodogram = PDM_search(mjd_r, m_r, P_low_sub, P_up_sub, 0., n_periods_sub, n_bin_sub, min_bin)
      best_period = periods[~np.isnan(periodogram)][np.argmin(periodogram[~np.isnan(periodogram)])]
      best_period = 0.16508

      if period_error == 'yes':
         for i in range(0,n_steps_sub):
            m_r_temp = [np.random.normal(m_r[j],sigmam_r[j]) for j in range(0,np.size(m_r))]
            m_r_temp = np.array(m_r_temp)
            periods_temp, periodogram_temp = PDM_search(mjd_r, m_r_temp, P_low_sub, P_up_sub, 0., n_periods_sub, n_bin_sub, min_bin)
            best_period_temp[i] = periods_temp[~np.isnan(periodogram_temp)][np.argmin(periodogram_temp[~np.isnan(periodogram_temp)])]
         sigma_best_period = np.std(best_period_temp,ddof=1)
   elif per_search_sub == 'LS2':
      periods, periodogram = LombScargle_search(mjd_r, m_r, sigmam_r, P_low_sub, P_up_sub, n_periods_sub, 2)
      best_period = periods[~np.isnan(periodogram)][np.argmax(periodogram[~np.isnan(periodogram)])]

      if period_error == 'yes':
         #Light curve simulation following Boldin et al. (2012) with a gaussian distribution (m'_i=np.random.normal(m_i,sigmam_i)) to estimate the error on the period
         for i in range(0,n_steps_sub):
            #m_r_temp = m_r + sigmam_r*np.random.uniform(-1.,1.,np.size(m_r))
            m_r_temp = [np.random.normal(m_r[j],sigmam_r[j]) for j in range(0,np.size(m_r))]
            m_r_temp = np.array(m_r_temp)
            periods_temp, periodogram_temp = LombScargle_search(mjd_r, m_r_temp, sigmam_r, P_low_sub, P_up_sub, n_periods_sub, 2)
            best_period_temp[i] = periods_temp[~np.isnan(periodogram_temp)][np.argmax(periodogram_temp[~np.isnan(periodogram_temp)])]
         sigma_best_period = np.std(best_period_temp,ddof=1)
   else:
      periods, periodogram = PDM_trisearch(mjd_g, m_g, mjd_r, m_r, mjd_i, m_i, P_low_sub, P_up_sub, 0., n_periods_sub, n_bin_sub, min_bin)
      best_period = periods[~np.isnan(periodogram)][np.argmin(periodogram[~np.isnan(periodogram)])]

      if period_error == 'yes':
         for i in range(0,n_steps_sub):
            m_g_temp = [np.random.normal(m_g[j],sigmam_g[j]) for j in range(0,np.size(m_g))]
            m_r_temp = [np.random.normal(m_r[j],sigmam_r[j]) for j in range(0,np.size(m_r))]
            m_i_temp = [np.random.normal(m_i[j],sigmam_i[j]) for j in range(0,np.size(m_i))]
            m_g_temp = np.array(m_g_temp)
            m_r_temp = np.array(m_r_temp)
            m_i_temp = np.array(m_i_temp)
            periods_temp, periodogram_temp = PDM_trisearch(mjd_g, m_g_temp, mjd_r, m_r_temp, mjd_i, m_i_temp, P_low_sub, P_up_sub, 0., n_periods_sub, n_bin_sub, min_bin)
            best_period_temp[i] = periods_temp[~np.isnan(periodogram_temp)][np.argmin(periodogram_temp[~np.isnan(periodogram_temp)])]
         sigma_best_period = np.std(best_period_temp,ddof=1)

   #t_0 = mjd_r[np.argmax(m_r)]
   t_0 = mjd_r[np.argsort(m_r)][-3]
   #t_0 = mjd_r[np.argmin(m_r)]-0.5*best_period
   #t_0 = 57343.939
   print("Epoch", t_0)
   print("Best period", best_period)
   if period_error == 'yes':
      print("Uncertainty on the best period", sigma_best_period)
   if 'mjd_g' in locals():
      if double_period == 'yes':
         phase_g = ((mjd_g-t_0)/(2*best_period)) % 1
      else:
         phase_g = ((mjd_g-t_0)/(best_period)) % 1

      #m_g_bin, bin_edges, binnumber = binned_statistic(phase_g[~np.isnan(m_g)],m_g[~np.isnan(m_g)],statistic='mean',bins=50)
      #phase_g_bin = bin_edges[1:] - (bin_edges[1] - bin_edges[0])/2
      #counts, bin_edges, binnumber = binned_statistic(phase_g[~np.isnan(m_g)],m_g[~np.isnan(m_g)],statistic='count',bins=50)
      #sigmam_g_bin, bin_edges, binnumber = binned_statistic(phase_g[~np.isnan(m_g)],m_g[~np.isnan(m_g)],statistic='std',bins=50)
      #sigmam_g_bin = sigmam_g_bin/np.sqrt(counts-1)
      
      if analysis == 'ZTF':
         if lc_rebin == 'yes':
            m_g_bin, bin_edges, binnumber = binned_statistic(phase_g[~np.isnan(m_g)],m_g[~np.isnan(m_g)],statistic='mean',bins=rebin_factor)
            phase_g_bin = bin_edges[1:] - (bin_edges[1] - bin_edges[0])/2
            counts, bin_edges, binnumber = binned_statistic(phase_g[~np.isnan(m_g)],m_g[~np.isnan(m_g)],statistic='count',bins=rebin_factor)
            sigmam_g_bin, bin_edges, binnumber = binned_statistic(phase_g[~np.isnan(m_g)],m_g[~np.isnan(m_g)],statistic='std',bins=rebin_factor)
            sigmam_g_bin = sigmam_g_bin/np.sqrt(counts-1)
   if 'mjd_r' in locals():
      if double_period == 'yes':
         phase_r = ((mjd_r-t_0)/(2*best_period)) % 1
      else:
         phase_r = ((mjd_r-t_0)/(best_period)) % 1

      #m_r_bin, bin_edges, binnumber = binned_statistic(phase_r[~np.isnan(m_r)],m_r[~np.isnan(m_r)],statistic='mean',bins=50)
      #phase_r_bin = bin_edges[1:] - (bin_edges[1] - bin_edges[0])/2
      #counts, bin_edges, binnumber = binned_statistic(phase_r[~np.isnan(m_r)],m_r[~np.isnan(m_r)],statistic='count',bins=50)
      #sigmam_r_bin, bin_edges, binnumber = binned_statistic(phase_r[~np.isnan(m_r)],m_r[~np.isnan(m_r)],statistic='std',bins=50)
      #sigmam_r_bin = sigmam_r_bin/np.sqrt(counts-1)
      if analysis == 'ZTF':
         if lc_rebin == 'yes':
            m_r_bin, bin_edges, binnumber = binned_statistic(phase_r[~np.isnan(m_r)],m_r[~np.isnan(m_r)],statistic='mean',bins=rebin_factor)
            phase_r_bin = bin_edges[1:] - (bin_edges[1] - bin_edges[0])/2
            counts, bin_edges, binnumber = binned_statistic(phase_r[~np.isnan(m_r)],m_r[~np.isnan(m_r)],statistic='count',bins=rebin_factor)
            sigmam_r_bin, bin_edges, binnumber = binned_statistic(phase_r[~np.isnan(m_r)],m_r[~np.isnan(m_r)],statistic='std',bins=rebin_factor)
            sigmam_r_bin = sigmam_r_bin/np.sqrt(counts-1)
   if 'mjd_i' in locals():
      if double_period == 'yes':
         phase_i = ((mjd_i-t_0)/(2*best_period)) % 1
      else:
         phase_i = ((mjd_i-t_0)/(best_period)) % 1

      #m_i_bin, bin_edges, binnumber = binned_statistic(phase_i[~np.isnan(m_i)],m_i[~np.isnan(m_i)],statistic='mean',bins=15)
      #phase_i_bin = bin_edges[1:] - (bin_edges[1] - bin_edges[0])/2
      #counts, bin_edges, binnumber = binned_statistic(phase_i[~np.isnan(m_i)],m_i[~np.isnan(m_i)],statistic='count',bins=15)
      #sigmam_i_bin, bin_edges, binnumber = binned_statistic(phase_i[~np.isnan(m_i)],m_i[~np.isnan(m_i)],statistic='std',bins=15)
      #sigmam_i_bin = sigmam_i_bin/np.sqrt(counts-1)
   if 'mjd_g' in locals() and 'mjd_r' in locals():
      if analysis == 'LCO' or analysis == 'STELLA':
         mjd_gr = []
         g_r = []
         sigma_gr = []
         for i in range(0,np.size(mjd_g)):
            for j in range(0,np.size(mjd_r)):
               if (mjd_r[j]-mjd_g[i]) < t_space*deltat/86400. and (mjd_r[j]-mjd_g[i])>0:
                  mjd_gr.append((mjd_g[i]+mjd_r[j])/2.)
                  g_r.append(m_g[i]-m_r[j])
                  sigma_gr.append(np.sqrt(sigmam_g[i]**2+sigmam_r[j]**2))
         mjd_gr = np.asarray(mjd_gr)
         g_r = np.asarray(g_r)
         sigma_gr = np.asarray(sigma_gr)
         if double_period == 'yes':
            phase_gr = ((mjd_gr-t_0)/(2*best_period)) % 1
         else:
            phase_gr = ((mjd_gr-t_0)/(best_period)) % 1
         #phase_gr = []
         #g_r = []
         #sigma_gr = []
         #for i in range(0,np.size(phase_g)-1):
         #   for j in range(0,np.size(phase_r)):
         #      if phase_r[np.argsort(phase_r)][j] > phase_g[np.argsort(phase_g)][i] and phase_r[np.argsort(phase_r)][j] < phase_g[np.argsort(phase_g)][i+1]:
         #         if (phase_r[np.argsort(phase_r)][j]-phase_g[np.argsort(phase_g)][i])<(phase_g[np.argsort(phase_g)][i+1]-phase_r[np.argsort(phase_r)][j]):
         #            phase_gr.append((phase_g[np.argsort(phase_g)][i]+phase_r[np.argsort(phase_r)][j])/2.)
         #            g_r.append(m_g[np.argsort(phase_g)][i]-m_r[np.argsort(phase_r)][j])
         #            sigma_gr.append(np.sqrt(sigmam_g[np.argsort(phase_g)][i]**2+sigmam_r[np.argsort(phase_r)][j]**2))
         #         else:
         #            phase_gr.append((phase_g[np.argsort(phase_g)][i+1]+phase_r[np.argsort(phase_r)][j])/2.)
         #            g_r.append(m_g[np.argsort(phase_g)][i+1]-m_r[np.argsort(phase_r)][j])
         #            sigma_gr.append(np.sqrt(sigmam_g[np.argsort(phase_g)][i+1]**2+sigmam_r[np.argsort(phase_r)][j]**2))
         #phase_gr = np.asarray(phase_gr)
         #g_r = np.asarray(g_r)
         #sigma_gr = np.asarray(sigma_gr)

         #g_r_bin, bin_edges, binnumber = binned_statistic(phase_gr[~np.isnan(g_r)],g_r[~np.isnan(g_r)],statistic='mean',bins=70)
         #phase_gr_bin = bin_edges[1:] - (bin_edges[1] - bin_edges[0])/2
         #counts, bin_edges, binnumber = binned_statistic(phase_gr[~np.isnan(g_r)],g_r[~np.isnan(g_r)],statistic='count',bins=70)
         #sigma_gr_bin, bin_edges, binnumber = binned_statistic(phase_gr[~np.isnan(g_r)],g_r[~np.isnan(g_r)],statistic='std',bins=70)
         #sigma_gr_bin = sigma_gr_bin/np.sqrt(counts-1)
         #phase_gr = phase_gr_bin
         #g_r = g_r_bin
         #sigma_gr = sigma_gr_bin
         #phase_g = phase_g_bin
         #m_g = m_g_bin
         #sigmam_g = sigmam_g_bin
         #phase_r = phase_r_bin
         #m_r = m_r_bin
         #sigmam_r = sigmam_r_bin
      else:
         phase_gr = []
         g_r = []
         sigma_gr = []
         for i in range(0,np.size(phase_g)-1):
            for j in range(0,np.size(phase_r)):
               if phase_r[np.argsort(phase_r)][j] > phase_g[np.argsort(phase_g)][i] and phase_r[np.argsort(phase_r)][j] < phase_g[np.argsort(phase_g)][i+1]:
                  if (phase_r[np.argsort(phase_r)][j]-phase_g[np.argsort(phase_g)][i])<(phase_g[np.argsort(phase_g)][i+1]-phase_r[np.argsort(phase_r)][j]):
                     phase_gr.append((phase_g[np.argsort(phase_g)][i]+phase_r[np.argsort(phase_r)][j])/2.)
                     g_r.append(m_g[np.argsort(phase_g)][i]-m_r[np.argsort(phase_r)][j])
                     sigma_gr.append(np.sqrt(sigmam_g[np.argsort(phase_g)][i]**2+sigmam_r[np.argsort(phase_r)][j]**2))
                  else:
                     phase_gr.append((phase_g[np.argsort(phase_g)][i+1]+phase_r[np.argsort(phase_r)][j])/2.)
                     g_r.append(m_g[np.argsort(phase_g)][i+1]-m_r[np.argsort(phase_r)][j])
                     sigma_gr.append(np.sqrt(sigmam_g[np.argsort(phase_g)][i+1]**2+sigmam_r[np.argsort(phase_r)][j]**2))
         phase_gr = np.asarray(phase_gr)
         g_r = np.asarray(g_r)
         sigma_gr = np.asarray(sigma_gr)
                  
      if analysis == 'ZTF':
         if lc_rebin == 'yes':
            g_r_bin, bin_edges, binnumber = binned_statistic(phase_gr[~np.isnan(g_r)],g_r[~np.isnan(g_r)],statistic='mean',bins=rebin_factor)
            phase_gr_bin = bin_edges[1:] - (bin_edges[1] - bin_edges[0])/2
            counts, bin_edges, binnumber = binned_statistic(phase_gr[~np.isnan(g_r)],g_r[~np.isnan(g_r)],statistic='count',bins=rebin_factor)
            sigma_gr_bin, bin_edges, binnumber = binned_statistic(phase_gr[~np.isnan(g_r)],g_r[~np.isnan(g_r)],statistic='std',bins=rebin_factor)
            sigma_gr_bin = sigma_gr_bin/np.sqrt(counts-1)
            phase_gr = phase_gr_bin
            g_r = g_r_bin
            sigma_gr = sigma_gr_bin
            phase_g = phase_g_bin
            m_g = m_g_bin
            sigmam_g = sigmam_g_bin
            phase_r = phase_r_bin
            m_r = m_r_bin
            sigmam_r = sigmam_r_bin
         
   if 'mjd_r' in locals() and 'mjd_i' in locals():
      mjd_ri = []
      r_i = []
      sigma_ri = []
      for i in range(0,np.size(mjd_r)):
         for j in range(0,np.size(mjd_i)):
            if (mjd_i[j]-mjd_r[i]) < t_space*deltat/86400. and (mjd_i[j]-mjd_r[i])>0:
               mjd_ri.append((mjd_r[i]+mjd_i[j])/2.)
               r_i.append(m_r[i]-m_i[j])
               sigma_ri.append(np.sqrt(sigmam_r[i]**2+sigmam_i[j]**2))
      mjd_ri = np.asarray(mjd_ri)
      r_i = np.asarray(r_i)
      sigma_ri = np.asarray(sigma_ri)
      if double_period == 'yes':
         phase_ri = ((mjd_ri-t_0)/(2*best_period)) % 1
      else:
         phase_ri = ((mjd_ri-t_0)/(best_period)) % 1

      #r_i_bin, bin_edges, binnumber = binned_statistic(phase_ri[~np.isnan(r_i)],r_i[~np.isnan(r_i)],statistic='mean',bins=15)
      #phase_ri_bin = bin_edges[1:] - (bin_edges[1] - bin_edges[0])/2
      #counts, bin_edges, binnumber = binned_statistic(phase_ri[~np.isnan(r_i)],r_i[~np.isnan(r_i)],statistic='count',bins=15)
      #sigma_ri_bin, bin_edges, binnumber = binned_statistic(phase_ri[~np.isnan(r_i)],r_i[~np.isnan(r_i)],statistic='std',bins=15)
      #sigma_ri_bin = sigma_ri_bin/np.sqrt(counts-1)

   ##ztf_t_0 = ztf_mjd_r[np.argmax(ztf_m_r)]
   #ztf_t_0 = ztf_mjd_r[np.argsort(ztf_m_r)][-10]
   #print("ZTF epoch", ztf_t_0)
   #if double_period == 'yes':
   #   ztf_phase_r = ((ztf_mjd_r-ztf_t_0)/(2*best_period)) % 1
   #else:
   #   ztf_phase_r = ((ztf_mjd_r-ztf_t_0)/(best_period)) % 1
   #ztf_m_r_bin, bin_edges, binnumber = binned_statistic(ztf_phase_r[~np.isnan(ztf_m_r)],ztf_m_r[~np.isnan(ztf_m_r)],statistic='mean',bins=int(np.round(np.size(ztf_m_r)/10)))
   #ztf_phase_r_bin = bin_edges[1:] - (bin_edges[1] - bin_edges[0])/2
   #counts, bin_edges, binnumber = binned_statistic(ztf_phase_r[~np.isnan(ztf_m_r)],ztf_m_r[~np.isnan(ztf_m_r)],statistic='count',bins=int(np.round(np.size(ztf_m_r)/10)))
   #ztf_sigmam_r_bin, bin_edges, binnumber = binned_statistic(ztf_phase_r[~np.isnan(ztf_m_r)],ztf_m_r[~np.isnan(ztf_m_r)],statistic='std',bins=int(np.round(np.size(ztf_m_r)/10)))
   #ztf_sigmam_r_bin = ztf_sigmam_r_bin/np.sqrt(counts-1)
   #ztf_phase_r = ztf_phase_r_bin
   #ztf_m_r = ztf_m_r_bin
   #ztf_sigmam_r = ztf_sigmam_r_bin
   
   figure, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 7), gridspec_kw={'height_ratios': [3, 1.5]}, sharex=True, dpi=300)
   #figure, (ax1, ax2) = plt.subplots(2, 1, figsize=(4, 3), gridspec_kw={'height_ratios': [4, 1.5]}, sharex=True, dpi=300)
   #figure, ax1 = plt.subplots(figsize=(5, 4), dpi=150)
   #plt.title(r""+FoVname.split('L')[0]+"L "+FoVname.split('L')[1]+"-A, "+telescope+" data",fontsize=12)
   plt.subplot(2, 1, 1)
   plt.ylabel(r'Apparent magnitude',fontsize=axlabelsize)
   #plt.ylabel(r'Apparent magnitude',fontsize=12)
   plt.gca().invert_yaxis()
   if 'mjd_g' in locals():
      plt.errorbar(phase_g, m_g, sigmam_g, fmt='o', markerfacecolor='cyan', markeredgecolor='cyan', ecolor='cyan', capsize=0, markersize=4.0, label=r"$\mathit{g'}$ band")
      #plt.errorbar(phase_g+1., m_g, sigmam_g, fmt='o', markerfacecolor='cyan', markeredgecolor='cyan', ecolor='cyan', capsize=0, markersize=2.5, label='_nolegend_')
      #plt.errorbar(phase_g_bin, m_g_bin, sigmam_g_bin, fmt='o', markerfacecolor='cyan', markeredgecolor='cyan', ecolor='cyan', capsize=0, markersize=4.0, label=r"$\mathit{g'}$ band")
      #plt.errorbar(phase_g_bin+1., m_g_bin, sigmam_g_bin, fmt='o', markerfacecolor='cyan', markeredgecolor='cyan', ecolor='cyan', capsize=0, markersize=4.0, label='_nolegend_')
      #plt.errorbar(phase_g, m_g, sigmam_g, fmt='o', markerfacecolor='cyan', markeredgecolor='cyan', ecolor='cyan', capsize=0, markersize=2.5, label=r"g' band")
   if 'mjd_r' in locals():
      plt.errorbar(phase_r, m_r, sigmam_r, fmt='o', markerfacecolor='lime', markeredgecolor='lime', ecolor='lime', capsize=0, markersize=4.0, label=r"$\mathit{r'}$ band")
      #plt.errorbar(phase_r+1., m_r, sigmam_r, fmt='o', markerfacecolor='lime', markeredgecolor='lime', ecolor='lime', capsize=0, markersize=2.5, label='_nolegend_')
      #plt.errorbar(phase_r_bin, m_r_bin, sigmam_r_bin, fmt='o', markerfacecolor='lime', markeredgecolor='lime', ecolor='lime', capsize=0, markersize=4.0, label=r"$\mathit{r'}$ band")
      #plt.errorbar(phase_r_bin+1., m_r_bin, sigmam_r_bin, fmt='o', markerfacecolor='lime', markeredgecolor='lime', ecolor='lime', capsize=0, markersize=4.0, label='_nolegend_')
      #plt.errorbar(phase_r, m_r, sigmam_r, fmt='o', markerfacecolor='lime', markeredgecolor='lime', ecolor='lime', capsize=0, markersize=2.5, label=r"r' band")
      #plt.errorbar(ztf_phase_r, ztf_m_r, ztf_sigmam_r, fmt='o', markerfacecolor='none', markeredgecolor='green', ecolor='green', capsize=0, markersize=4.0, alpha=0.6, label=r"ZTF $\mathit{r'}$ band")
      #plt.errorbar(ztf_phase_r+1., ztf_m_r, ztf_sigmam_r, fmt='o', markerfacecolor='none', markeredgecolor='green', ecolor='green', capsize=0, markersize=4.0, alpha=0.6, label='_nolegend_')
   if 'mjd_i' in locals():
      plt.errorbar(phase_i, m_i, sigmam_i, fmt='o', markerfacecolor='fuchsia', markeredgecolor='fuchsia', ecolor='fuchsia', capsize=0, markersize=4.0, label=r"$\mathit{i'}$ band")
      #plt.errorbar(phase_i+1., m_i, sigmam_i, fmt='o', markerfacecolor='fuchsia', markeredgecolor='fuchsia', ecolor='fuchsia', capsize=0, markersize=2.5, label='_nolegend_')
      #plt.errorbar(phase_i_bin, m_i_bin, sigmam_i_bin, fmt='o', markerfacecolor='fuchsia', markeredgecolor='fuchsia', ecolor='fuchsia', capsize=0, markersize=4.0, label=r"$\mathit{i'}$ band")
      #plt.errorbar(phase_i_bin+1., m_i_bin, sigmam_i_bin, fmt='o', markerfacecolor='fuchsia', markeredgecolor='fuchsia', ecolor='fuchsia', capsize=0, markersize=4.0, label='_nolegend_')
      #plt.errorbar(phase_i, m_i, sigmam_i, fmt='o', markerfacecolor='fuchsia', markeredgecolor='fuchsia', ecolor='fuchsia', capsize=0, markersize=2.5, label=r"i' band")
   ax1.xaxis.set_minor_locator(tck.AutoMinorLocator())
   ax1.yaxis.set_minor_locator(tck.AutoMinorLocator())
   #ax1.set_ylim([14.75,13.5])
   ax1.tick_params(axis='both', which='both', direction='in', top=True, right=True,labelsize=axticksize)
   #ax1.tick_params(axis='both', which='both', direction='in', top=True, right=True,labelsize=12)
   #plt.xlabel(r'Phase ($P_{\mathrm{best}}$=%.4f'%np.round(best_period,4)+'(5) d)',fontsize=12)
   ax1.legend(loc=(0.,0.22),frameon=False,fontsize=subsize)
   #ax1.legend(loc=(-0.03,0.0),frameon=False,fontsize=12)
      

   plt.subplot(2, 1, 2)
   ##plt.xlabel(r'Orbital phase',fontsize=axlabelsize)
   plt.ylabel(r'Color',fontsize=axlabelsize)
   ##plt.xlabel(r'Orbital phase')
   ##plt.ylabel(r'Color')
   plt.gca().invert_yaxis()
   if 'mjd_g' in locals() and 'mjd_r' in locals():
      plt.errorbar(phase_gr, g_r, sigma_gr, fmt='o', color='blue', capsize=0, markersize=4.0, label=r"$\mathit{g'}$-$\mathit{r'}$")
      #plt.errorbar(phase_gr+1., g_r, sigma_gr, fmt='o', color='blue', capsize=0, markersize=4.0, label='_nolegend_')
      #plt.errorbar(phase_gr_bin, g_r_bin, sigma_gr_bin, fmt='o', color='blue', capsize=0, markersize=4.0, label=r"$\mathit{g'}$-$\mathit{r'}$")
      #plt.errorbar(phase_gr_bin+1., g_r_bin, sigma_gr_bin, fmt='o', color='blue', capsize=0, markersize=4.0, label='_nolegend_')
   #   #plt.errorbar(phase_gr, g_r, sigma_gr, fmt='o', color='blue', capsize=0, markersize=2.5, label=r"(g'-r')")
   if 'mjd_r' in locals() and 'mjd_i' in locals():
      plt.errorbar(phase_ri, r_i, sigma_ri, fmt='o', color='red', capsize=0, markersize=4.0, label=r"$\mathit{r'}$-$\mathit{i'}$")
      #plt.errorbar(phase_ri+1., r_i, sigma_ri, fmt='o', color='red', capsize=0, markersize=4.0, label='_nolegend_')
      #plt.errorbar(phase_ri_bin, r_i_bin, sigma_ri_bin, fmt='o', color='red', capsize=0, markersize=4.0, label=r"$\mathit{r'}$-$\mathit{i'}$")
   #   #plt.errorbar(phase_ri_bin+1., r_i_bin, sigma_ri_bin, fmt='o', color='red', capsize=0, markersize=4.0, label='_nolegend_')
   #   #plt.errorbar(phase_ri, r_i, sigma_ri, fmt='o', color='red', capsize=0, markersize=2.5, label=r"(r'-i')")
   ax2.xaxis.set_minor_locator(tck.AutoMinorLocator())
   ax2.yaxis.set_minor_locator(tck.AutoMinorLocator())
   ax2.tick_params(axis='both', which='both', direction='in', top=True, right=True,labelsize=axticksize)
   ax2.legend(loc=(0.,0.23),frameon=False,fontsize=subsize)
   ##ax2.legend(loc=(0.,0.),frameon=False,fontsize=9)
   if double_period == 'yes':
      if not np.isnan(max(periodogram)) and not np.isinf(max(periodogram)) and not np.isnan(min(periodogram)) and not np.isinf(min(periodogram)) and not np.isnan(best_period) and not np.isinf(best_period):
         if period_error == 'yes':
            error_pow = int("{:e}".format(2*sigma_best_period).split('e')[1][2])
            if int("{:e}".format(2*sigma_best_period).split('e')[0][2])>5:
               error_digit = int("{:e}".format(2*sigma_best_period).split('e')[0][0])+1
               if error_digit == 10:
                  error_pow = error_pow-1
                  error_digit = int(1)
            else:
               error_digit = int("{:e}".format(sigma_best_period).split('e')[0][0])
            plt.xlabel(r'Phase ($2\times P_{\mathrm{best}}$='+str(np.round(2*best_period,error_pow))+'('+str(error_digit)+') d)',fontsize=axlabelsize)
         else:
            plt.xlabel(r'Phase ($2\times P_{\mathrm{best}}$=%.5f'%np.round(2*best_period,5)+'(6) d)',fontsize=axlabelsize)
   else:
      if period_error == 'yes':
         error_pow = int("{:e}".format(sigma_best_period).split('e')[1][2])
         if int("{:e}".format(sigma_best_period).split('e')[0][2])>5:
            error_digit = int("{:e}".format(sigma_best_period).split('e')[0][0])+1
            if error_digit == 10:
               error_pow = error_pow-1
               error_digit = int(1)
         else:
            error_digit = int("{:e}".format(sigma_best_period).split('e')[0][0])
         plt.xlabel(r'Phase ($P_{\mathrm{best}}$='+str(np.round(best_period,error_pow))+'('+str(error_digit)+') d)',fontsize=axlabelsize)
      else:
         #plt.xlabel(r'Phase ($P_{\mathrm{best}}$=%.3f'%np.round(best_period,3)+'(7) d)',fontsize=axlabelsize)
         plt.xlabel(r'Orbital phase',fontsize=axlabelsize)

   #broken axis
   #plt.subplots_adjust(wspace=0, hspace=0)
   #ax1.spines.bottom.set_visible(False)
   #ax1.xaxis.tick_top()
   #ax1.tick_params(labeltop=True)
   #ax1.tick_params(top=True)
   #ax1.tick_params(axis='x', which='minor', top=True)
   #ax1.xaxis.tick_bottom()
   #ax1.tick_params(labelbottom=True)
   #ax1.tick_params(bottom=True)
   #ax1.tick_params(axis='x', which='minor', bottom=True)

   plt.subplots_adjust(wspace=0, hspace=0)
   ax1.spines.bottom.set_visible(False)
   ax1.xaxis.tick_top()
   ax1.tick_params(labeltop=True)
   ax2.xaxis.tick_bottom()
   ax2.tick_params(top=True)
   ax2.tick_params(axis='x', which='minor', top=True)
   
   if analysis == 'LCO' or analysis == 'STELLA':
      if double_period == 'yes':
         plt.savefig(path+'/finalperiodplots_paper/'+FoVname+'_'+str(ncand)+'_lc_doubleP.pdf',bbox_inches='tight',pad_inches=0)
      else:
         plt.savefig(path+'/finalperiodplots_paper/'+FoVname+'_'+str(ncand)+'_lc_defence.png',bbox_inches='tight',pad_inches=0)
   else:
      if double_period == 'yes':
         plt.savefig(path+'/finalperiodplots_paper/'+FoVname+'_'+str(ncand)+'_lc_doubleP_ZTF.pdf',bbox_inches='tight',pad_inches=0)
      else:
         plt.savefig(path+'/finalperiodplots_paper/'+FoVname+'_'+str(ncand)+'_lc_ZTF_defence.png',bbox_inches='tight',pad_inches=0)
   plt.close(figure)

   print("##### OUR DATA #####")
   if 'mjd_g' in locals():
      print("Peak-to-peak amplitudes in g", max(m_g)-min(m_g))
      print("Mid and range uncertainty of mag g:", (max(m_g)+min(m_g))/2, (max(m_g)-min(m_g))/2)
   if 'mjd_r' in locals():
      print("Peak-to-peak amplitudes in r", max(m_r)-min(m_r))
      print("Mid and range uncertainty of mag r:", (max(m_r)+min(m_r))/2, (max(m_r)-min(m_r))/2)
   if 'mjd_i' in locals():
      print("Peak-to-peak amplitudes in i", max(m_i)-min(m_i))
      print("Mid and range uncertainty of mag i:", (max(m_i)+min(m_i))/2, (max(m_i)-min(m_i))/2)
   print("Phase corresponding to r brightness peak:", phase_r[np.argmin(m_r)], "Phase corresponding to r brightness minimum:", phase_r[np.argmax(m_r)])
   if 'mjd_g' in locals() and 'mjd_r' in locals():
      print("Color (g-r) observed correspondingly to r brightness peak with error:", g_r[np.argmin(np.fabs(phase_gr-phase_r[np.argmin(m_r)]))], sigma_gr[np.argmin(np.fabs(phase_gr-phase_r[np.argmin(m_r)]))])
      print("Color (g-r) observed correspondingly to r brightness minimum with error:", g_r[np.argmin(np.fabs(phase_gr-phase_r[np.argmax(m_r)]))], sigma_gr[np.argmin(np.fabs(phase_gr-phase_r[np.argmax(m_r)]))])
      print("Phase corresponding to peak temperature in (g-r):", phase_gr[np.argmin(g_r)])
      print("Phase corresponding to minimum temperature in (g-r):", phase_gr[np.argmax(g_r)])
      print("Color (g-r) corresponding to the peak temperature with error:", min(g_r), sigma_gr[np.argmin(g_r)])
      print("Color (g-r) corresponding to the minimum temperature with error:", max(g_r), sigma_gr[np.argmax(g_r)])
      print("Mean and standard deviation of the color (g-r):", np.mean(g_r), np.std(g_r,ddof=1))
   if 'mjd_r' in locals() and 'mjd_i' in locals():
      print("Color (r-i) observed correspondingly to r brightness peak with error:", r_i[np.argmin(np.fabs(phase_ri-phase_r[np.argmin(m_r)]))], sigma_ri[np.argmin(np.fabs(phase_ri-phase_r[np.argmin(m_r)]))])
      print("Color (r-i) observed correspondingly to r brightness minimum with error:", r_i[np.argmin(np.fabs(phase_ri-phase_r[np.argmax(m_r)]))], sigma_ri[np.argmin(np.fabs(phase_ri-phase_r[np.argmax(m_r)]))])
      print("Phase corresponding to peak temperature in (r-i):", phase_ri[np.argmin(r_i)])
      print("Phase corresponding to minimum temperature in (r-i):", phase_ri[np.argmax(r_i)])
      print("Color (r-i) corresponding to the peak temperature with error:", min(r_i), sigma_ri[np.argmin(r_i)])
      print("Color (r-i) corresponding to the minimum temperature with error:", max(r_i), sigma_ri[np.argmax(r_i)])
      print("Mean and standard deviation of the color (r-i):", np.mean(r_i), np.std(r_i,ddof=1))
      
   if PS_plot == 'yes':
      figure, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 7), gridspec_kw={'height_ratios': [3, 1.5]}, dpi=300)
      #figure.suptitle(r""+FoVname.split('L')[0]+"L "+FoVname.split('L')[1],y=0.95,fontsize=tsize)
      plt.subplot(2, 1, 1)
      plt.plot(periods,periodogram, color='black', linewidth=1.5)
      if not np.isnan(max(periodogram)) and not np.isinf(max(periodogram)) and not np.isnan(min(periodogram)) and not np.isinf(min(periodogram)) and not np.isnan(best_period) and not np.isinf(best_period):
         if period_error == 'yes':
            error_pow = int("{:e}".format(sigma_best_period).split('e')[1][2])
            if int("{:e}".format(sigma_best_period).split('e')[0][2])>5:
               error_digit = int("{:e}".format(sigma_best_period).split('e')[0][0])+1
               if error_digit == 10:
                  error_pow = error_pow-1
                  error_digit = int(1)
            else:
               error_digit = int("{:e}".format(sigma_best_period).split('e')[0][0])
            plt.axvline(x=best_period, color='r', ls='-', linewidth=1.0, label=r"$P_{\mathrm{best}}$="+str(np.round(best_period,error_pow))+"("+str(error_digit)+") d")
         else:
            plt.axvline(x=best_period, color='r', ls='-', linewidth=1.0, label=r"$P_{\mathrm{best}}$=%.5f"%np.round(best_period,5)+"(6) d")
      #rcParams['xtick.minor.size'] = axticksize
      ax1.set_xscale('log')
      #ax1.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
      ax1.set_xlim(plt.xlim()[0]+0.004,plt.xlim()[1]-0.2)
      if per_search_sub == 'LS1':
         #plt.ylabel("Lomb-Scargle power in "+r"$\mathit{r'}$",fontsize=axlabelsize)
         plt.ylabel("Lomb-Scargle power",fontsize=axlabelsize)
         if period_FAPl == 'yes':
            tic_1 = time.perf_counter()
            FAP_periods, FAP_levels = LombScargle_FAP_level(mjd_r, m_r, sigmam_r, P_low_sub, P_up_sub, n_periods_sub, 1, FAPl_factor, conf_level)
            toc_1 = time.perf_counter()
            print(f"FAP levels computed in {toc_1 - tic_1:0.4f} seconds")
            #FAP_periods = periods
            #FAP_levels = 0.1*np.ones(np.size(periods))
            plt.plot(FAP_periods,FAP_levels,color='orange',linewidth=1.0,ls='--')
            plt.text(plt.xlim()[0]*1.7, np.mean(FAP_levels)*1.2, r'0.1% FAP', color='darkorange', ha='center',fontsize=size)
      elif per_search_sub == 'LS2':
         plt.ylabel("Lomb-Scargle power",fontsize=axlabelsize)
         if period_FAPl == 'yes':
            #tic_1 = time.perf_counter()
            #FAP_periods, FAP_levels = LombScargle_FAP_level(mjd_r, m_r, sigmam_r, P_low_sub, P_up_sub, n_periods_sub, 2, FAPl_factor, conf_level)
            #toc_1 = time.perf_counter()
            #print(f"FAP levels computed in {toc_1 - tic_1:0.4f} seconds")
            FAP_periods = periods
            FAP_levels = 0.02*np.ones(np.size(periods))
            plt.plot(FAP_periods,FAP_levels,color='orange',linewidth=1.0,ls='--')
            plt.text(plt.xlim()[0]*1.7, np.mean(FAP_levels)*1.2, r'0.1% FAP', color='darkorange', ha='center',fontsize=size)
      elif per_search_sub == 'PDM':
         #plt.ylabel("Normalized phase dispersion in "+r"$\mathit{r'}$",fontsize=axlabelsize)
         plt.ylabel("Normalized phase dispersion",fontsize=axlabelsize)
         if period_FAPl == 'yes':
            tic_1 = time.perf_counter()
            FAP_periods, FAP_levels = PDM_FAP_level(mjd_r, m_r, P_low_sub, P_up_sub, 0., n_periods_sub, n_bin_sub, min_bin, FAPl_factor, conf_level)
            toc_1 = time.perf_counter()
            print(f"FAP levels computed in {toc_1 - tic_1:0.4f} seconds")
            #FAP_periods = periods
            #FAP_levels = 0.97*np.ones(np.size(periods))
            plt.plot(FAP_periods,FAP_levels,color='orange',linewidth=1.5,ls='--')
            plt.text(plt.xlim()[0]*1.7, 0.97*np.mean(FAP_levels), r'0.1% FAP', color='darkorange', ha='center',fontsize=size)
      ax1.yaxis.set_minor_locator(tck.AutoMinorLocator())
      ticks = [0.1, 0.9]
      ax1.set_xticks(ticks)
      formatter = FuncFormatter(lambda x, _: f"{x:g}")
      ax1.xaxis.set_major_formatter(formatter)
      ax1.tick_params(axis='both', which='both', direction='in', top=True, right=True,labelsize=axticksize)
      ax1.set_xlabel(r"Period (d)",fontsize=axlabelsize)
      ax1.legend(loc=(0.,0.45),frameon=False,fontsize=size)

      plt.subplot(2, 1, 2)
      #print(best_period,"{:e}".format(sigma_best_period))
      ax2.set_xlim([0.1,0.5])
      ax2.set_ylim([0.72,1.])
      plt.plot(periods,periodogram, color='black', linewidth=1.5)
      plt.axvline(x=best_period, color='r', ls='-', linewidth=1.0)
      #ax2.set_xscale('log')
      #ax2.tick_params(axis='x', labelsize=axticksize)
      #ax2.tick_params(axis='y', labelsize=axticksize)
      #if per_search_sub == 'LS1' or per_search_sub == 'LS2':
      #   plt.ylabel("Power in "+r"$\mathit{r'}$",fontsize=axlabelsize)
      #elif per_search_sub == 'PDM':
      #   plt.ylabel("Normalized variance in "+r"$\mathit{r'}$",fontsize=axlabelsize)
      if period_FAPl == 'yes':
         plt.plot(FAP_periods,FAP_levels,color='orange',linewidth=1.5,ls='--')
         #plt.text(plt.xlim()[0]*1.4, np.mean(FAP_levels)*1.02, r'0.1% FAP', color='orange', ha='center',fontsize=axticksize)
      #else:
      #   #ax1.set_title("PDM multi-band data periodogram",fontsize=subsize)
      #   ax2.set_ylabel("Sum of normalized variances",fontsize=axlabelsize)
      ##if period_FAP == 'yes':
      ##   if per_search_sub == 'LS1':
      ##      FAP = LombScargle_FAP(mjd_r, m_r, sigmam_r, P_low_sub, P_up_sub, n_periods_sub, 1, n_bootstraps_sub)
      ##      print("False-Alarm-Probability:", FAP)
      ##      ax2.text(0.1,0.3,'FAP='+str(np.round(FAP,3)),fontsize=subsize)
      ##   elif per_search_sub == 'PDM':
      ##      FAP = PDM_FAP(mjd_r, m_r, P_low_sub, P_up_sub, 0., n_periods_sub, n_bin_sub, min_bin, n_bootstraps_sub)
      ##      print("False-Alarm-Probability:", FAP)
      ##      ax2.text(0.1,0.3,'FAP='+str(np.round(FAP,3)),fontsize=subsize)
      ##   else:
      ##      ax2.text(0.1,0.3,'FAP=0.001',fontsize=subsize)
      ##else:
      ##   ax2.text(0.105,0.035,r'FAP$<10^{-4}$',color='red',fontsize=subsize)
      ax2.yaxis.set_minor_locator(tck.AutoMinorLocator())
      ax2.xaxis.set_minor_locator(tck.AutoMinorLocator())
      ticks = [0.1,0.2,0.3,0.4,0.5]
      ax2.set_xticks(ticks)
      formatter = FuncFormatter(lambda x, _: f"{x:g}")
      ax2.xaxis.set_major_formatter(formatter)
      ax2.tick_params(axis='both', which='both', direction='in', top=True, right=True,labelsize=axticksize)
      #ax2.set_xlabel(r"Period (d)",fontsize=axlabelsize)
      
      #broken axis
      plt.subplots_adjust(hspace=0.25)
      ax1.tick_params(top=True)
      ax2.xaxis.tick_bottom()
      ax2.tick_params(axis='x', which='minor', bottom=True)

      if analysis == 'LCO' or analysis == 'STELLA':
         plt.savefig(path+'/finalperiodplots_paper/'+FoVname+'_'+str(ncand)+'_PS_'+per_search_sub+'.pdf',bbox_inches='tight',pad_inches=0)
      else:
         plt.savefig(path+'/finalperiodplots_paper/'+FoVname+'_'+str(ncand)+'_PS_'+per_search_sub+'_ZTF.pdf',bbox_inches='tight',pad_inches=0)
      plt.close(figure)

#Plot target's FoV in the r band, with all the variable stars found inside
#STELLA FoV
figure, ax = plt.subplots(figsize=(4., 4.), dpi=300)
#INT FoV
#figure, ax = plt.subplots(figsize=(3.5, 5), dpi=300)
axticksize = 8
axlabelsize = 8 
subsize = 10
tsize = 10

if os.path.exists(path+'/median_r_after.fits') and os.path.exists(path+'/final_periodic_can_1.txt'):
   varSources_fin1=np.genfromtxt(path+'/rfilter/varSources_1.csv', delimiter=',')
   list_targ = np.int_(np.genfromtxt(path+'/final_periodic_can_1.txt'))
   nvar=np.size(list_targ)
   if nvar!=1:
      perSources_fin = varSources_fin1[list_targ-1]
      ra=perSources_fin[:,0]
      dec=perSources_fin[:,1]
   else:
      perSources_fin = varSources_fin1[list_targ-1,:]
      ra=perSources_fin[0]
      dec=perSources_fin[1]
   x, y = w.all_world2pix(ra, dec, 1)
   
elif not os.path.exists(path+'/median_r_after.fits') and not os.path.exists(path+'/final_periodic_can_1.txt') and os.path.exists(path+'/final_periodic_can_ZTF.txt'):
   varSources_fin1=np.genfromtxt(path+'/ZTF_search/ZTF_coord.dat')
   list_targ = np.int_(np.genfromtxt(path+'/final_periodic_can_ZTF.txt'))
   nvar=np.size(list_targ)
   if nvar!=1:
      perSources_fin = varSources_fin1[list_targ-1]
      ra=perSources_fin[:,0]
      dec=perSources_fin[:,1]
   else:
      perSources_fin = varSources_fin1[list_targ-1,:]
      ra=perSources_fin[0]
      dec=perSources_fin[1]
   x, y = w.all_world2pix(ra, dec, 1)

ax = plt.subplot(projection=w)
#plt.title(r""+FoVname.split('L')[0]+"L "+FoVname.split('L')[1]+" "+telescope+" field",fontsize=tsize)
#ax.coords[0].set_ticks(number=4)
right_asc = ax.coords[0]
decl = ax.coords[1]
#right_asc = ax.coords[1]
#decl = ax.coords[0]
#right_asc.set_axislabel('R. A. (J2000)')
#decl.set_axislabel('Declination (J2000)')
plt.xlabel('R. A. (J2000)',fontsize=axlabelsize)
plt.ylabel('Declination (J2000)',fontsize=axlabelsize,labelpad=-0.8)
#plt.ylabel('Declination (J2000)',fontsize=axlabelsize,labelpad=+0.35)
#plt.xlabel('Declination (J2000)',fontsize=axlabelsize)
#plt.ylabel('R. A. (J2000)',fontsize=axlabelsize,labelpad=-0.5)
   
#plt.xlim(max(0,xtar-0.45*bound*max(ell_a_95,ell_b_95)), min(np.shape(data)[1]-1,xtar+0.97*bound*max(ell_a_95,ell_b_95)))
#plt.ylim(max(0,ytar-0.9*bound*max(ell_a_95,ell_b_95)), min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)))
plt.xlim(max(0,xtar-0.42*bound*max(ell_a_95,ell_b_95)), min(np.shape(data)[1]-1,xtar+0.42*bound*max(ell_a_95,ell_b_95)))
plt.ylim(max(0,ytar-0.85*bound*max(ell_a_95,ell_b_95)), min(np.shape(data)[0]-1,ytar+0.42*bound*max(ell_a_95,ell_b_95)))
#plt.xlim(max(0,xtar-bound*max(ell_a_95,ell_b_95)), xtar+117)
#plt.ylim(ytar-200, ytar+250)
#plt.xlim(x-550,x+500)
#plt.ylim(y-800,y+500)
#plt.xlim(xtar-bound, xtar+bound)
#plt.ylim(ytar-bound, ytar+bound)
#plt.imshow(data, interpolation='nearest', cmap='gray', vmin=m-colscale*s, vmax=m+colscale*s, origin='lower')

right_asc.set_ticks_position('lb') 
right_asc.set_ticklabel_position('b')
#right_asc.set_ticklabel_position('l')
decl.set_ticks_position('lb')
decl.set_ticklabel_position('l')
#decl.set_ticklabel_position('b') 
right_asc.tick_params(labelsize=axticksize)
decl.tick_params(labelsize=axticksize)
zscale = ZScaleInterval()
plt.imshow(zscale(data), cmap="gray")
plt.text(xtar+40,ytar,FoVname.split('L')[0]+'L '+FoVname.split('L')[1],horizontalalignment='center',verticalalignment='center',color='gold',fontsize=tsize)
target = Ellipse(xy=(xtar, ytar), width=2*ell_a_95, height=2*ell_b_95, angle=-ell_theta_95)
target.set_facecolor('none')
target.set_edgecolor('gold')
ax.add_artist(target)
#if np.size(sources_3FGL)!=0:
#   for i in range(0,np.shape(sources_3FGL)[0]):
#      if x_3FGL[i]>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_3FGL[i]<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_3FGL[i]>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_3FGL[i]<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
#         plt.text(x_3FGL[i]-15,y_3FGL[i],sources_3FGL[i].split('L')[0]+'L '+sources_3FGL[i].split('L')[1],horizontalalignment='center',verticalalignment='center',color='yellow',fontsize=tsize)
#         a = Ellipse(xy=(x_3FGL[i], y_3FGL[i]), width=2*3600./float(platescale)*coord_3FGL[i,3], height=2*3600./float(platescale)*coord_3FGL[i,4], angle=coord_3FGL[i,5]+90)
#         a.set_facecolor('none')
#         a.set_edgecolor('yellow')
#         ax.add_artist(a)
#if np.size(sources_4FGLe)!=0:
#   for i in range(0,np.shape(sources_4FGLe)[0]):
#      if x_4FGLe[i]>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_4FGLe[i]<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_4FGLe[i]>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_4FGLe[i]<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
#         plt.text(x_4FGLe[i],y_4FGLe[i],'4FGLe',horizontalalignment='center',verticalalignment='center',color='cyan',fontsize=tsize,alpha=0.8)
#         a = Ellipse(xy=(x_4FGLe[i], y_4FGLe[i]), width=2*3600./float(platescale)*np.asarray(coord_4FGLe)[i,2], height=2*3600./float(platescale)*np.asarray(coord_4FGLe)[i,3], angle=np.asarray(coord_4FGLe)[i,4]+90)
#         a.set_facecolor('none')
#         a.set_edgecolor('cyan')
#         ax.add_artist(a)

if 'ra' in locals():
   if nvar == 1:
      if os.path.exists(path+'/2sxps.cat'):
         if np.size(x_Xcparts_Swift)!=1:
            for i in range(0,np.size(x_Xcparts_Swift)):
               if SkyCoord(ra, dec, frame='fk5', unit='deg').separation(SkyCoord(ra_Xcparts_Swift[i],dec_Xcparts_Swift[i],frame='fk5',unit='deg')).arcsecond <= radius_X:
                  a = Circle(xy=(x_Xcparts_Swift[i], y_Xcparts_Swift[i]), radius=err_Xcparts_Swift[i], alpha=alpha)
                  a.set_linewidth(0.5)
                  a.set_facecolor('none')
                  plt.text(x_Xcparts_Swift[i]+35,y_Xcparts_Swift[i]-45,'SWIFT',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=tsize)
                  a.set_edgecolor('fuchsia')
                  ax.add_artist(a)
         else:
            if SkyCoord(ra, dec, frame='fk5', unit='deg').separation(SkyCoord(ra_Xcparts_Swift,dec_Xcparts_Swift,frame='fk5',unit='deg')).arcsecond <= radius_X:
               a = Circle(xy=(x_Xcparts_Swift, y_Xcparts_Swift), radius=err_Xcparts_Swift, alpha=alpha)
               a.set_linewidth(0.5)
               a.set_facecolor('none')
               plt.text(x_Xcparts_Swift,y_Xcparts_Swift-100,'SWIFT',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=subsize)
               a.set_edgecolor('fuchsia')
               ax.add_artist(a)

      if os.path.exists(path+'/csc2.cat'):
         if np.size(x_Xcparts_Chandra)!=1:
            for i in range(0,np.size(x_Xcparts_Chandra)):
               if SkyCoord(ra, dec, frame='fk5', unit='deg').separation(SkyCoord(ra_Xcparts_Chandra[i],dec_Xcparts_Chandra[i],frame='fk5',unit='deg')).arcsecond <= radius_X:
                  a = Ellipse(xy=(x_Xcparts_Chandra[i], y_Xcparts_Chandra[i]), width=2*ell_a_Chandra[i], height=2*ell_b_Chandra[i], angle=ell_theta_Chandra[i], alpha=alpha)
                  a.set_linewidth(0.5)
                  a.set_facecolor('none')
                  plt.text(x_Xcparts_Chandra[i],y_Xcparts_Chandra[i]+50,'CHANDRA',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=subsize)
                  a.set_edgecolor('fuchsia')
                  ax.add_artist(a)
         else:
            if SkyCoord(ra, dec, frame='fk5', unit='deg').separation(SkyCoord(ra_Xcparts_Chandra,dec_Xcparts_Chandra,frame='fk5',unit='deg')).arcsecond <= radius_X:
               a = Ellipse(xy=(x_Xcparts_Chandra, y_Xcparts_Chandra), width=2*ell_a_Chandra, height=2*ell_b_Chandra, angle=ell_theta_Chandra, alpha=alpha)
               a.set_linewidth(0.5)
               a.set_facecolor('none')
               plt.text(x_Xcparts_Chandra,y_Xcparts_Chandra+50,'CHANDRA',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=subsize)
               a.set_edgecolor('fuchsia')
               ax.add_artist(a)

      if os.path.exists(path+'/XMM4d13s.cat'):
         if np.size(x_Xcparts_XMM)!=1:
            for i in range(0,np.size(x_Xcparts_XMM)):
               if SkyCoord(ra, dec, frame='fk5', unit='deg').separation(SkyCoord(ra_Xcparts_XMM[i],dec_Xcparts_XMM[i],frame='fk5',unit='deg')).arcsecond <= radius_X:
                  a = Circle(xy=(x_Xcparts_XMM[i], y_Xcparts_XMM[i]), radius=err_Xcparts_XMM[i], alpha=alpha)
                  a.set_linewidth(0.5)
                  a.set_facecolor('none')
                  plt.text(x_Xcparts_XMM[i]-55,y_Xcparts_XMM[i]+12,'XMM',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.5*tsize)
                  a.set_edgecolor('fuchsia')
                  ax.add_artist(a)
         else:
            if SkyCoord(ra, dec, frame='fk5', unit='deg').separation(SkyCoord(ra_Xcparts_XMM,dec_Xcparts_XMM,frame='fk5',unit='deg')).arcsecond <= radius_X:
               a = Circle(xy=(x_Xcparts_XMM, y_Xcparts_XMM), radius=err_Xcparts_XMM, alpha=alpha)
               a.set_linewidth(0.5)
               a.set_facecolor('none')
               plt.text(x_Xcparts_XMM+50,y_Xcparts_XMM,'XMM',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=subsize)
               a.set_edgecolor('fuchsia')
               ax3.add_artist(a)

      if os.path.exists(path+'/eRASS1d1.cat'):
         if np.size(x_Xcparts_eRASS)!=1:
            for i in range(0,np.size(x_Xcparts_eRASS)):
               if SkyCoord(ra, dec, frame='fk5', unit='deg').separation(SkyCoord(ra_Xcparts_eRASS[i],dec_Xcparts_eRASS[i],frame='fk5',unit='deg')).arcsecond <= radius_X:
                  a = Circle(xy=(x_Xcparts_eRASS[i], y_Xcparts_eRASS[i]), radius=err_Xcparts_eRASS[i], alpha=alpha)
                  a.set_linewidth(0.5)
                  a.set_facecolor('none')
                  plt.text(x_Xcparts_eRASS[i]+15,y_Xcparts_eRASS[i],'eRASS',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=subsize)
                  a.set_edgecolor('fuchsia')
                  ax.add_artist(a)
         else:
            if SkyCoord(ra, dec, frame='fk5', unit='deg').separation(SkyCoord(ra_Xcparts_eRASS,dec_Xcparts_eRASS,frame='fk5',unit='deg')).arcsecond <= radius_X:
               a = Circle(xy=(x_Xcparts_eRASS, y_Xcparts_eRASS), radius=err_Xcparts_eRASS, alpha=alpha)
               a.set_linewidth(0.5)
               a.set_facecolor('none')
               plt.text(x_Xcparts_eRASS,y_Xcparts_eRASS+50,'eRASS',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=subsize)
               a.set_edgecolor('fuchsia')
               ax3.add_artist(a)
         
      plt.text(x+560,y+560,'CANDIDATE',horizontalalignment='center',verticalalignment='center',color='lime',fontsize=subsize)
      #plt.text(x-90,y+75,'OPTICAL',horizontalalignment='center',verticalalignment='center',color='darkorange',fontsize=tsize)
      a = Circle(xy=(x, y), radius=25, alpha=alpha)
      a.set_linewidth(1.0)
      a.set_facecolor('none')
      a.set_edgecolor('lime')
      ax.add_artist(a)
      a = FancyArrowPatch((x+500, y+500), (x+15, y+15), mutation_scale=18/1.25)
      #a = Arrow(x=x+100, y=y+100, dx=-100, dy=-100, width=1.0, alpha=alpha)
      a.set_facecolor('lime')
      a.set_edgecolor('lime')
      ax.add_artist(a)

   else:
      for i in range(0,nvar):
         if os.path.exists(path+'/2sxps.cat'):
            if np.size(x_Xcparts_Swift)!=1:
               for j in range(0,np.size(x_Xcparts_Swift)):
                  if SkyCoord(ra[i], dec[i], frame='fk5', unit='deg').separation(SkyCoord(ra_Xcparts_Swift[j],dec_Xcparts_Swift[j],frame='fk5',unit='deg')).arcsecond <= radius_X:
                     a = Circle(xy=(x_Xcparts_Swift[j], y_Xcparts_Swift[j]), radius=err_Xcparts_Swift[j], alpha=alpha)
                     a.set_linewidth(0.3)
                     a.set_facecolor('none')
                     plt.text(x_Xcparts_Swift[j]+50,y_Xcparts_Swift[j]-50,'SWIFT',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=tsize)
                     a.set_edgecolor('fuchsia')
                     ax.add_artist(a)
            else:
               if SkyCoord(ra[i], dec[i], frame='fk5', unit='deg').separation(SkyCoord(ra_Xcparts_Swift,dec_Xcparts_Swift,frame='fk5',unit='deg')).arcsecond <= radius_X:
                  a = Circle(xy=(x_Xcparts_Swift, y_Xcparts_Swift), radius=err_Xcparts_Swift, alpha=alpha)
                  a.set_linewidth(0.5)
                  a.set_facecolor('none')
                  plt.text(x_Xcparts_Swift,y_Xcparts_Swift-50,'SWIFT',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=tsize)
                  a.set_edgecolor('fuchsia')
                  ax.add_artist(a)

         if os.path.exists(path+'/csc2.cat'):
            if np.size(x_Xcparts_Chandra)!=1:
               for j in range(0,np.size(x_Xcparts_Chandra)):
                  if SkyCoord(ra[i], dec[i], frame='fk5', unit='deg').separation(SkyCoord(ra_Xcparts_Chandra[j],dec_Xcparts_Chandra[j],frame='fk5',unit='deg')).arcsecond <= radius_X:
                     a = Ellipse(xy=(x_Xcparts_Chandra[j], y_Xcparts_Chandra[j]), width=2*ell_a_Chandra[j], height=2*ell_b_Chandra[j], angle=ell_theta_Chandra[j], alpha=alpha)
                     a.set_linewidth(0.5)
                     a.set_facecolor('none')
                     plt.text(x_Xcparts_Chandra[j],y_Xcparts_Chandra[j]+50,'CHANDRA',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=subsize)
                     a.set_edgecolor('fuchsia')
                     ax.add_artist(a)
            else:
               if SkyCoord(ra[i], dec[i], frame='fk5', unit='deg').separation(SkyCoord(ra_Xcparts_Chandra,dec_Xcparts_Chandra,frame='fk5',unit='deg')).arcsecond <= radius_X:
                  a = Ellipse(xy=(x_Xcparts_Chandra, y_Xcparts_Chandra), width=2*ell_a_Chandra, height=2*ell_b_Chandra, angle=ell_theta_Chandra, alpha=alpha)
                  a.set_linewidth(0.5)
                  a.set_facecolor('none')
                  plt.text(x_Xcparts_Chandra,y_Xcparts_Chandra+50,'CHANDRA',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=subsize)
                  a.set_edgecolor('fuchsia')
                  ax.add_artist(a)

         if os.path.exists(path+'/XMM4d13s.cat'):
            if np.size(x_Xcparts_XMM)!=1:
               for j in range(0,np.size(x_Xcparts_XMM)):
                  if SkyCoord(ra[i], dec[i], frame='fk5', unit='deg').separation(SkyCoord(ra_Xcparts_XMM[j],dec_Xcparts_XMM[j],frame='fk5',unit='deg')).arcsecond <= radius_X:
                     a = Circle(xy=(x_Xcparts_XMM[j], y_Xcparts_XMM[j]), radius=err_Xcparts_XMM[j], alpha=alpha)
                     a.set_linewidth(0.5)
                     a.set_facecolor('none')
                     plt.text(x_Xcparts_XMM[j]-50,y_Xcparts_XMM[j]+20,'XMM',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.5*tsize)
                     a.set_edgecolor('fuchsia')
                     ax.add_artist(a)
            else:
               if SkyCoord(ra[i], dec[i], frame='fk5', unit='deg').separation(SkyCoord(ra_Xcparts_XMM,dec_Xcparts_XMM,frame='fk5',unit='deg')).arcsecond <= radius_X:
                  a = Circle(xy=(x_Xcparts_XMM, y_Xcparts_XMM), radius=err_Xcparts_XMM, alpha=alpha)
                  a.set_linewidth(0.5)
                  a.set_facecolor('none')
                  plt.text(x_Xcparts_XMM,y_Xcparts_XMM+50,'XMM',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=subsize)
                  a.set_edgecolor('fuchsia')
                  ax3.add_artist(a)

         if os.path.exists(path+'/eRASS1d1.cat'):
            if np.size(x_Xcparts_eRASS)!=1:
               for j in range(0,np.size(x_Xcparts_eRASS)):
                  if SkyCoord(ra[i], dec[i], frame='fk5', unit='deg').separation(SkyCoord(ra_Xcparts_eRASS[j],dec_Xcparts_eRASS[j],frame='fk5',unit='deg')).arcsecond <= radius_X:
                     a = Circle(xy=(x_Xcparts_eRASS[j], y_Xcparts_eRASS[j]), radius=err_Xcparts_eRASS[j], alpha=alpha)
                     a.set_linewidth(0.5)
                     a.set_facecolor('none')
                     plt.text(x_Xcparts_eRASS[j],y_Xcparts_eRASS[j]+50,'eRASS',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=subsize)
                     a.set_edgecolor('fuchsia')
                     ax.add_artist(a)
            else:
               if SkyCoord(ra[i], dec[i], frame='fk5', unit='deg').separation(SkyCoord(ra_Xcparts_eRASS,dec_Xcparts_eRASS,frame='fk5',unit='deg')).arcsecond <= radius_X:
                  a = Circle(xy=(x_Xcparts_eRASS, y_Xcparts_eRASS), radius=err_Xcparts_eRASS, alpha=alpha)
                  a.set_linewidth(0.5)
                  a.set_facecolor('none')
                  plt.text(x_Xcparts_eRASS,y_Xcparts_eRASS+50,'eRASS',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=subsize)
                  a.set_edgecolor('fuchsia')
                  ax3.add_artist(a)
         if i==0:
            plt.text(x[i]-90,y[i]+210,'CANDIDATE J2117-A',horizontalalignment='center',verticalalignment='center',color='lime',fontsize=tsize)
            a = Circle(xy=(x[i], y[i]), radius=12, alpha=alpha)
            a.set_linewidth(1.0)
            a.set_facecolor('none')
            a.set_edgecolor('lime')
            ax.add_artist(a)
            a = FancyArrowPatch((x[i]+100, y[i]+180), (x[i]+7, y[i]+7), mutation_scale=18/1.25)
             #a = Arrow(x=x+100, y=y+100, dx=-100, dy=-100, width=1.0, alpha=alpha)
            a.set_facecolor('lime')
            a.set_edgecolor('lime')
            ax.add_artist(a)
         elif i==1:
            plt.text(x[i]+40,y[i]+210,'CANDIDATE J2117-B',horizontalalignment='center',verticalalignment='center',color='lime',fontsize=tsize)
            a = Circle(xy=(x[i], y[i]), radius=12, alpha=alpha)
            a.set_linewidth(1.0)
            a.set_facecolor('none')
            a.set_edgecolor('lime')
            ax.add_artist(a)
            a = FancyArrowPatch((x[i]+100, y[i]+180), (x[i]+7, y[i]+7), mutation_scale=18/1.25)
             #a = Arrow(x=x+100, y=y+100, dx=-100, dy=-100, width=1.0, alpha=alpha)
            a.set_facecolor('lime')
            a.set_edgecolor('lime')
            ax.add_artist(a)
         elif i==2:
            plt.text(x[i]+130,y[i]+130,'J1748-C',horizontalalignment='center',verticalalignment='center',color='darkorange',fontsize=tsize)
            a = Circle(xy=(x[i], y[i]), radius=11, alpha=alpha)
            a.set_linewidth(1.0)
            a.set_facecolor('none')
            a.set_edgecolor('darkorange')
            ax.add_artist(a)
            a = FancyArrowPatch((x[i]+120, y[i]+120), (x[i]+5, y[i]+5), mutation_scale=18/1.25)
            #a = Arrow(x=x+100, y=y+100, dx=-100, dy=-100, width=1.0, alpha=alpha)
            a.set_facecolor('darkorange')
            a.set_edgecolor('darkorange')
            ax.add_artist(a)
         elif i==3:
            plt.text(x[i]-65,y[i]-145,'J1748-D',horizontalalignment='center',verticalalignment='center',color='darkorange',fontsize=tsize)
            a = Circle(xy=(x[i], y[i]), radius=11, alpha=alpha)
            a.set_linewidth(1.0)
            a.set_facecolor('none')
            a.set_edgecolor('darkorange')
            ax.add_artist(a)
            a = FancyArrowPatch((x[i]-120, y[i]-120), (x[i]-5, y[i]-5), mutation_scale=18/1.25)
            #a = Arrow(x=x+100, y=y+100, dx=-100, dy=-100, width=1.0, alpha=alpha)
            a.set_facecolor('darkorange')
            a.set_edgecolor('darkorange')
            ax.add_artist(a)
         elif i==4:
            plt.text(x[i]-190,y[i]-75,'CANDIDATE',horizontalalignment='center',verticalalignment='center',color='lime',fontsize=tsize)
            a = Circle(xy=(x[i], y[i]), radius=11, alpha=alpha)
            a.set_linewidth(1.0)
            a.set_facecolor('none')
            a.set_edgecolor('lime')
            ax.add_artist(a)
            a = FancyArrowPatch((x[i]-165, y[i]-50), (x[i]-5, y[i]-5), mutation_scale=18/1.25)
            #a = Arrow(x=x+100, y=y+100, dx=-100, dy=-100, width=1.0, alpha=alpha)
            a.set_facecolor('lime')
            a.set_edgecolor('lime')
            ax.add_artist(a)

plt.savefig(path+'/finalperiodplots_paper/'+FoVname+'_FoV_defence.png',bbox_inches='tight',pad_inches=0)
plt.close(figure)
