#!/bin/python

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib import rcParams
from matplotlib.patches import Circle
from matplotlib.patches import Ellipse
from astropy.io import fits
from astropy import wcs
from astropy.coordinates import SkyCoord
from astropy.coordinates import Angle
from astropy.timeseries import LombScargle
import astropy.units as u
from ztfquery import lightcurve
import requests
import swifttools.ukssdc.data.SXPS as uds
import mastcasjobs
from astroquery.simbad import Simbad
import io
import os
import shutil
import time

from periodicity_search import LombScargle_search, LombScargle_FAP_level, PDM_search, PDM_bisearch, PDM_trisearch, PDM_FAP_level

path=os.getcwd()

#We get as reference epoch for period folding the time (in MJD) of the first image in the gband
#t_0 = float(np.genfromtxt(path+'/gfilter/usedImages.txt',dtype='str')[0].split('_')[2])

analysis = input("Enter 'STELLA', 'INT' or 'LCO' for the type of data (default value 'LCO'):\n") or 'LCO'
while analysis!='STELLA' and analysis!='INT' and analysis!='LCO':
   analysis = input("Invalid value, please insert 'STELLA', 'INT' or 'LCO':\n") or 'LCO'
fext = float(input("Enter the scale factor for the FWHM used for the extraction radius (default value 1.2):\n") or 1.2)
while fext<=0.:
   fext = float(input("Invalid value, please insert a positive extraction factor:\n") or 1.2)
nlist = int(input("Enter the list of candidate variables you want to inspect (1 or 2) (default value 1):\n") or 1)
while nlist!=1 and nlist!=2:
   nlist = int(input("Invalid value, please insert 1 or 2:\n") or 1)
eRASS_search = input("Enter 'yes' if you want to look for eRASS X-ray counterparts (default value 'no', it can take several minutes because of the eRASS catalog size if this is the first time you query for this field):\n")  or "no"
while eRASS_search!='yes' and eRASS_search!='no':
   eRASS_search = input("Invalid string, please insert 'yes' or 'no':\n") or "no"
ATLAS = input("Do you want ATLAS forced photometry? (default 'yes'):\n") or 'yes'
while ATLAS!='yes' and ATLAS!='no':
   ATLAS = input("Invalid string, please insert 'yes' or 'no':\n") or 'yes'

per_search = input("Enter 'LS', 'PDM' or 'multiPDM' for the periodicity search method (default value 'PDM'):\n") or 'PDM'
while per_search!='LS' and per_search!='PDM' and per_search!='multiPDM':
   per_search = input("Invalid value, please insert 'LS', 'PDM' or 'multiPDM':\n") or 'PDM'
n_periods_sub = float(input("Enter an integer positive number for number of period steps (10-10000) to use (default value 2000):\n")  or "2000")
while n_periods_sub<10 or n_periods_sub>10000.:
   n_periods_sub = float(input("Invalid value, please insert a integer number between 10 and 10000:\n") or "2000")
n_periods_sub = int(n_periods_sub)
n_periods = float(input("Enter an integer positive number for number of period steps (10-10000) to use on ZTF data with better sampling (default value 2000):\n")  or "2000")
while n_periods<10 or n_periods>10000.:
   n_periods = float(input("Invalid value, please insert a integer number between 10 and 10000:\n") or "2000")
n_periods = int(n_periods)
P_low  = float(input("Enter the start period (in days) used (default value 0.02):\n") or 0.02)
while P_low<=0.:
   P_low = float(input("Invalid value, please insert a positive start period:\n") or 0.02)
P_up  = float(input("Enter the final period (in days) used (default value 2.5):\n") or 2.5)
while P_up<=0.:
   P_up = float(input("Invalid value, please insert a positive final period:\n") or 2.5)

if per_search == 'PDM' or per_search == 'multiPDM':
   n_bin_sub = float(input("Enter an integer positive number for number of phase bins to use in the PDM (or multiband-PDM) method (default value 10):\n")  or "10")
   while n_bin_sub<0. or n_bin_sub%1!=0.:
      n_bin_sub = float(input("Invalid value, please insert a positive integer number (default value 10):\n") or "10")
   n_bin_sub = int(n_bin_sub)
   n_bin = float(input("Enter an integer positive number for number of phase bins to use on ZTF data in the PDM (or multiband-PDM) method (default value 10):\n")  or "10")
   while n_bin<0. or n_bin%1!=0.:
      n_bin = float(input("Invalid value, please insert a positive integer number:\n") or "10")
   n_bin = int(n_bin)

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

#Defining scale factor (in sigma) to use as squared region around the target, with side=2*bound*semimaj_ax
bound = 2
#Defining radius size (in arcmin) to use for plotting 3FGL sources nearby to our candidate target
radius_gamma = 11
#Defining radius size (in arcsec) to consider association of an optical candidate with an X-ray counterpart
radius_X = 3
#Defining scale factor (in pixel) to use as squared zoom region around the variable source, with side=2*zoom
zoom = 100
#Color range factor to be multiplied to data stdev in pyplot.imshow
colscale = 0.15
#Number of decimal digits to which round the coordinates for not plotting twice the aperture radius of the variable candidate in the last panel (zoomed FoV)
ndig=5
#Built-in,minimum number of data points per phase bin (if lower bin is not used in PDM calculation)
min_bin = 2
#Time spacing factor for colour light curve
#t_space = 1.5
t_space = 2.0
markersize=8
markersize_main=10
alpha=0.25
tsize = 35
subsize = 22
axlabelsize = 18
axticksize = 16
#figxsize1 = 25
#figysize1 = 35
#figxsize2 = 30
#figysize2 = 42
figxsize = 35
figysize = 52
#figxsize = 40
#figysize = 65

tic = time.perf_counter()
list_targ = np.int_(np.genfromtxt(path+'/final_periodic_can_'+str(nlist)+'.txt'))
nvar=np.size(list_targ)
#nvar=1
sourcesr = np.genfromtxt(path+'/rfilter/starVariability.csv', delimiter=',')
varSources_fin1=np.genfromtxt(path+'/rfilter/varSources_1.csv', delimiter=',')
#varSources_fin1=np.genfromtxt(path+'/varSources_1_'+var+'.csv', delimiter=',')
#varSources_fin2=np.genfromtxt(path+'/rfilter/varSources_2.csv', delimiter=',')
#if os.path.exists(path+'/varSources_2_'+var+'.csv'):
#   varSources_fin2=np.genfromtxt(path+'/varSources_2_'+var+'.csv', delimiter=',')
#else:
#   varSources_fin2=np.genfromtxt(path+'/rfilter/varSources_2.csv', delimiter=',')

if not (os.path.isdir(path+'/finalperiodplots'+str(nlist)+'_'+per_search)):
   os.system('mkdir '+path+'/finalperiodplots'+str(nlist)+'_'+per_search)

if nlist==1:
   if nvar!=1:
      perSources_fin = varSources_fin1[list_targ-1]
   else:
      perSources_fin = varSources_fin1[list_targ-1,:]
else:
   if nvar!=1:
     perSources_fin = varSources_fin2[list_targ-1]
   else:
      perSources_fin = varSources_fin2[list_targ-1,:]

images = np.genfromtxt(path+'/rfilter.txt',dtype='str')
hdul = fits.open(path+'/median_r_after.fits')
data = hdul[0].data
hdr = hdul[0].header
rx = hdul[1].data['x']
ry = hdul[1].data['y']
hdul.close()
if analysis=='STELLA':
   FoVname = hdr['OBJNAME']
   deltat = float(hdr['EXPT'])
   platescale = hdr['PIXSCALE']
   telescope = 'STELLA/WiFSIP'
elif analysis=='INT':
   FoVname = hdr['OBJECT']
   deltat = float(hdr['EXPTIME'])
   platescale = hdr['SECPPIX']
   telescope = 'INT/WFC'
elif analysis=='LCO':
   FoVname = hdr['OBJECT'].replace("L_","L")
   deltat = float(hdr['EXPTIME'])
   platescale = hdr['PIXSCALE']
   hdul = fits.open(path+'/rfilter/'+images[0])
   if hdul[0].header['INSTRUME']=='kb84':
      telescope = 'LCO/SBIG'
   else:
      telescope = 'LCO/Sinistro'
   hdul.close()
w = wcs.WCS(hdr)

#Load the corresponding Fermi coordinates, 2sigma ellipses (in degrees) and 1 sigma semi-major axis (in degrees) of the target
hdul = fits.open('/export/work/marcotu/gll_psc_v28.fit')
[[ratar, dectar, ell_a_68, ell_a_95, ell_b_95, ell_theta_95]] = [[hdul[1].data['RAJ2000'][i],hdul[1].data['DEJ2000'][i],hdul[1].data['Conf_68_SemiMajor'][i],3600./float(platescale)*hdul[1].data['Conf_95_SemiMajor'][i],3600./float(platescale)*hdul[1].data['Conf_95_SemiMinor'][i],hdul[1].data['Conf_95_PosAng'][i]] for i in range(0,len(hdul[1].data)) if FoVname.replace("L","L ")==hdul[1].data['Source_Name'][i]]
xtar, ytar = w.all_world2pix(ratar, dectar, 1)
hdul.close()

#Load the 3FGL Fermi updated coordinates for our target and eventual closeby sources (within radius arcmin from our target), also extended sources
hdul = fits.open('/export/work/marcotu/gll_psc_v16.fit')
sources_3FGL = [hdul[1].data['Source_Name'][i] for i in range(0,len(hdul[1].data))]
coord_3FGL = [[hdul[1].data['RAJ2000'][i],hdul[1].data['DEJ2000'][i],hdul[1].data['Conf_68_SemiMajor'][i],hdul[1].data['Conf_95_SemiMajor'][i],hdul[1].data['Conf_95_SemiMinor'][i],hdul[1].data['Conf_95_PosAng'][i]] for i in range(0,len(hdul[1].data))]
sources_3FGL = np.asarray(sources_3FGL)
coord_3FGL = np.asarray(coord_3FGL)

sep = SkyCoord(ratar, dectar, frame='fk5', unit='deg').separation(SkyCoord(coord_3FGL[:,0],coord_3FGL[:,1],frame='fk5',unit='deg')).arcminute
sources_3FGL = [sources_3FGL[i] for i in range(0,len(sources_3FGL)) if sep[i]<=radius_gamma]
sources_3FGL = np.asarray(sources_3FGL)
coord_3FGL = [coord_3FGL[i,:] for i in range(0,len(coord_3FGL)) if sep[i]<=radius_gamma]
coord_3FGL = np.asarray(coord_3FGL)
sep = [sep[i] for i in range(0,len(sep)) if sep[i]<=radius_gamma]

if np.size(sources_3FGL)!=0:
   x_3FGL, y_3FGL = w.all_world2pix(coord_3FGL[:,0], coord_3FGL[:,1], 1)

#Load eventual 4FGL extended closeby sources (within radius_gamma arcmin from our target)
#sources_4FGLe = [hdul[2].data['Source_Name'][i] for i in range(0,len(hdul[2].data))]
#coord_4FGLe = [[hdul[2].data['RAJ2000'][i],hdul[2].data['DEJ2000'][i],hdul[2].data['Model_SemiMajor'][i],hdul[2].data['Model_SemiMinor'][i],hdul[2].data['Model_PosAng'][i]] for i in range(0,len(hdul[2].data))]
#sources_4FGLe = np.asarray(sources_4FGLe)
#coord_4FGLe = np.asarray(coord_4FGLe)
#sep = SkyCoord(ratar, dectar, frame='fk5', unit='deg').separation(SkyCoord(coord_4FGLe[:,0],coord_4FGLe[:,1],frame='fk5',unit='deg')).arcminute
#sources_4FGLe = [sources_4FGLe[i] for i in range(0,len(sources_4FGLe)) if sep[i]<=radius_gamma]
#coord_4FGLe = [coord_4FGLe[i] for i in range(0,len(coord_4FGLe)) if sep[i]<=radius_gamma]
#sep = [sep[i] for i in range(0,len(sep)) if sep[i]<=radius_gamma]
#sources_4FGLe = np.asarray(sources_4FGLe)
#coord_4FGLe = np.asarray(coord_4FGLe)
#if np.size(sources_4FGLe)!=0:
#   x_4FGLe, y_4FGLe = w.all_world2pix(np.asarray(coord_4FGLe)[:,0], np.asarray(coord_4FGLe)[:,1], 1)
hdul.close()

#Load eventual X-ray counterparts for the Field of View from Chandra (with csc2.py script), XMM-Newton (with XMM4d13s.py script) and Swift (with swift2SXPS.py script) and radio counterparts from NVSS, VLASS, FIRST and MSC 4FGL radio surveys
if os.path.exists(path+'/2sxps.cat'):
   os.remove(path+'/2sxps.cat')
if os.path.exists(path+'/2sxps.reg'):
   os.remove(path+'/2sxps.reg')
if os.path.exists(path+'/csc2.cat'):
   os.remove(path+'/csc2.cat')
if os.path.exists(path+'/csc2.reg'):
   os.remove(path+'/csc2.reg')
if os.path.exists(path+'/XMM4d13s.cat'):
   os.remove(path+'/XMM4d13s.cat')
if os.path.exists(path+'/XMM4d13s.reg'):
   os.remove(path+'/XMM4d13s.reg')
if os.path.exists(path+'/nvss.cat'):
   os.remove(path+'/nvss.cat')
if os.path.exists(path+'/nvss.reg'):
   os.remove(path+'/nvss.reg')
if os.path.exists(path+'/vlass.cat'):
   os.remove(path+'/vlass.cat')
if os.path.exists(path+'/vlass.reg'):
   os.remove(path+'/vlass.reg')
if os.path.exists(path+'/first.cat'):
   os.remove(path+'/first.cat')
if os.path.exists(path+'/first.reg'):
   os.remove(path+'/first.reg')
if os.path.exists(path+'/MSC_4FGL.cat'):
   os.remove(path+'/MSC_4FGL.cat')
if os.path.exists(path+'/MSC_4FGL.reg'):
   os.remove(path+'/MSC_4FGL.reg')
if os.path.exists(path+'/FASTgpps.cat'):
   os.remove(path+'/FASTgpps.cat')
if os.path.exists(path+'/FASTgpps.reg'):
   os.remove(path+'/FASTgpps.reg')


#Counts-flux conversion factor for X-ray upper limits, taken from the Heasoft tool Pimms,
#with power law index 1.2, in 0.3-10 keV energy range and computing the galactic nh from the "HI4PI collaboration 2016" 2D HI map
#counts_flux is the unabs. flux (in erg/s/cm^2) corresponding to 1 c/s acquired from Swift/XRT/PC
counts_flux = 5.929E-11

#Flux-unabsorbed flux conversion factor for X-ray Chandra/ACIS camera sources, taken from the Heasoft tool Pimms,
#with power law index 1.2, in 0.5-7 keV energy range and computing the galactic nh from the "HI4PI collaboration 2016" 2D HI map
#unabs_flux_Chandra is the unabs. flux (in erg/s/cm^2) corresponding to 1 erg/s/cm^2
unabs_flux_Chandra = 1.168
#Flux-unabsorbed flux conversion factor for X-ray XMM/EPIC/PN camera sources, taken from the Heasoft tool Pimms,
#with power law index 1.2, in 0.2-12 keV energy range and computing the galactic nh from the "HI4PI collaboration 2016" 2D HI map
#unabs_flux_XMM is the unabs. flux (in erg/s/cm^2) corresponding to 1 erg/s/cm^2
unabs_flux_XMM = 1.038
#Flux-unabsorbed flux conversion factor for X-ray sources of the eROSITA western main catalogue, taken from the Heasoft tool Pimms,
#with power law index 1.2, in 0.2-2.3 keV energy range and computing the galactic nh from the "HI4PI collaboration 2016" 2D HI map
#unabs_flux_eRASS is the unabs. flux (in erg/s/cm^2) corresponding to 1 erg/s/cm^2
unabs_flux_eRASS = 1.182

hdul = fits.open(path+'/rfilter/'+images[0])
hdr = hdul[0].header
hdul.close()
#if hdr['INSTRUME']=='WIFSIP V2.0':
#   df = pd.read_excel('/export/work/marcotu/COBIPULSE-North/STELLA/2015B/DATA/reduced/Radio_Xray_counterparts.xlsx',usecols=(0,2,4,6,7),dtype={'3FGL Field':str, 'RA':float, 'DEC':float, 'Type':str, 'Observatory':str})
#   Fields_cparts = df['3FGL Field'].to_numpy()
#   ra_cparts = df['RA'].to_numpy()
#   dec_cparts = df['DEC'].to_numpy()
#   type_cparts = df['Type'].to_numpy()
#   obs_Rcparts = df['Observatory'].to_numpy()
#   ra_cparts = [ra_cparts[i] for i in range(0,len(Fields_cparts)) if Fields_cparts[i]==FoVname.replace("3FGL","") and type_cparts[i]=='R']
#   dec_cparts = [dec_cparts[i] for i in range(0,len(Fields_cparts)) if Fields_cparts[i]==FoVname.replace("3FGL","") and type_cparts[i]=='R']
#   #type_cparts = [type_cparts[i] for i in range(0,len(Fields_cparts)) if Fields_cparts[i]==FoVname.replace("3FGL","")]
#   obs_Rcparts = [obs_Rcparts[i] for i in range(0,len(Fields_cparts)) if Fields_cparts[i]==FoVname.replace("3FGL","") and type_cparts[i]=='R']
#   #Fields_cparts = [Fields_cparts[i] for i in range(0,len(Fields_cparts)) if Fields_cparts[i]==FoVname.replace("3FGL","")]
#   if np.size(obs_Rcparts)!=0:
#      x_Rcparts, y_Rcparts = w.all_world2pix(ra_cparts, dec_cparts, 1)
#   del df, Fields_cparts, ra_cparts, dec_cparts, type_cparts

os.system('swift2SXPS.py '+str(ratar)+' '+str(dectar)+' '+str(4*ell_a_68*60.)+' -1 0')
if os.path.exists(path+'/2sxps.cat'):
   ra_Xcparts_Swift, dec_Xcparts_Swift, err_Xcparts_Swift, flux_X_Swift = np.genfromtxt(path+'/2sxps.cat',usecols=(2,3,4,8),unpack=True)
   x_Xcparts_Swift, y_Xcparts_Swift = w.all_world2pix(ra_Xcparts_Swift, dec_Xcparts_Swift, 1)
   err_Xcparts_Swift = err_Xcparts_Swift/float(platescale)

os.system('csc2.py '+str(ratar)+' '+str(dectar)+' '+str(4*ell_a_68*60.))
if os.path.exists(path+'/csc2.cat'):
   ra_Xcparts_Chandra, dec_Xcparts_Chandra, ell_a_Chandra, ell_b_Chandra, ell_theta_Chandra, flux_X_Chandra = np.genfromtxt(path+'/csc2.cat',usecols=(1,2,3,4,5,6),unpack=True)
   x_Xcparts_Chandra, y_Xcparts_Chandra = w.all_world2pix(ra_Xcparts_Chandra, dec_Xcparts_Chandra, 1)
   ell_a_Chandra = ell_a_Chandra/float(platescale)
   ell_b_Chandra = ell_b_Chandra/float(platescale)
   flux_X_Chandra = flux_X_Chandra*unabs_flux_Chandra

os.system('XMM4d13s.py '+str(ratar)+' '+str(dectar)+' '+str(4*ell_a_68*60.))
if os.path.exists(path+'/XMM4d13s.cat'):
   ra_Xcparts_XMM, dec_Xcparts_XMM, err_Xcparts_XMM, flux_X_XMM = np.genfromtxt(path+'/XMM4d13s.cat',usecols=(1,2,3,15),unpack=True)
   x_Xcparts_XMM, y_Xcparts_XMM = w.all_world2pix(ra_Xcparts_XMM, dec_Xcparts_XMM, 1)
   err_Xcparts_XMM = err_Xcparts_XMM/float(platescale)
   flux_X_XMM = flux_X_XMM*unabs_flux_XMM

if eRASS_search == 'yes':
   #Perform the eROSITA query just if the catalogue in this field wasn't already produced before, the query takes several minutes!
   if not os.path.exists(path+'/eRASS1d1.cat'):
      os.system('eRASS1d1.py '+str(ratar)+' '+str(dectar)+' '+str(4*ell_a_68*60.))
if os.path.exists(path+'/eRASS1d1.cat'):
   ra_Xcparts_eRASS, dec_Xcparts_eRASS, err_Xcparts_eRASS, flux_X_eRASS = np.genfromtxt(path+'/eRASS1d1.cat',usecols=(1,2,3,5),unpack=True)
   x_Xcparts_eRASS, y_Xcparts_eRASS = w.all_world2pix(ra_Xcparts_eRASS, dec_Xcparts_eRASS, 1)
   err_Xcparts_eRASS = err_Xcparts_eRASS/float(platescale)
   flux_X_eRASS = flux_X_eRASS*unabs_flux_eRASS

cat_R = []
x_Rcparts = []
y_Rcparts = []
errmaj_Rcparts = []
errmin_Rcparts = []
errangle_Rcparts = []
flux_R = []
os.system('nvss.py '+str(ratar)+' '+str(dectar)+' '+str(4*ell_a_68*60.))
if os.path.exists(path+'/nvss.cat'):
   ra_Rcparts_nvss, dec_Rcparts_nvss = np.genfromtxt(path+'/nvss.cat',usecols=(1,2),dtype='U',unpack=True)
   ra_Rcparts_nvss = Angle(ra_Rcparts_nvss,unit="hourangle").degree
   dec_Rcparts_nvss = Angle(dec_Rcparts_nvss,unit="deg").degree
   errmaj_Rcparts_nvss, errmin_Rcparts_nvss, errangle_Rcparts_nvss, flux_R_nvss = np.genfromtxt(path+'/nvss.cat',usecols=(3,4,5,6),unpack=True)
   x_Rcparts_nvss, y_Rcparts_nvss = w.all_world2pix(ra_Rcparts_nvss, dec_Rcparts_nvss, 1)
   errmaj_Rcparts_nvss = errmaj_Rcparts_nvss/float(platescale)
   errmin_Rcparts_nvss = errmin_Rcparts_nvss/float(platescale)
   if np.size(ra_Rcparts_nvss)>1:
      cat_R.extend(np.full(np.size(ra_Rcparts_nvss),'nvss'))
      x_Rcparts.extend(x_Rcparts_nvss)
      y_Rcparts.extend(y_Rcparts_nvss)
      errmaj_Rcparts.extend(errmaj_Rcparts_nvss)
      errmin_Rcparts.extend(errmin_Rcparts_nvss)
      errangle_Rcparts.extend(errangle_Rcparts_nvss)
   else:
      cat_R.append('nvss')
      x_Rcparts.append(x_Rcparts_nvss)
      y_Rcparts.append(y_Rcparts_nvss)
      errmaj_Rcparts.append(errmaj_Rcparts_nvss)
      errmin_Rcparts.append(errmin_Rcparts_nvss)
      errangle_Rcparts.append(errangle_Rcparts_nvss)

os.system('vlass.py '+str(ratar)+' '+str(dectar)+' '+str(4*ell_a_68*60.))
if os.path.exists(path+'/vlass.cat'):
   ra_Rcparts_vlass, dec_Rcparts_vlass, errmaj_Rcparts_vlass, errmin_Rcparts_vlass, errangle_Rcparts_vlass, flux_R_vlass = np.genfromtxt(path+'/vlass.cat',usecols=(2,3,4,5,6,7),unpack=True)
   x_Rcparts_vlass, y_Rcparts_vlass = w.all_world2pix(ra_Rcparts_vlass, dec_Rcparts_vlass, 1)
   errmaj_Rcparts_vlass = errmaj_Rcparts_vlass/float(platescale)
   errmin_Rcparts_vlass = errmin_Rcparts_vlass/float(platescale)
   if np.size(ra_Rcparts_vlass)>1:
      cat_R.extend(np.full(np.size(ra_Rcparts_vlass),'vlass'))
      x_Rcparts.extend(x_Rcparts_vlass)
      y_Rcparts.extend(y_Rcparts_vlass)
      errmaj_Rcparts.extend(errmaj_Rcparts_vlass)
      errmin_Rcparts.extend(errmin_Rcparts_vlass)
      errangle_Rcparts.extend(errangle_Rcparts_vlass)
   else:
      cat_R.append('vlass')
      x_Rcparts.append(x_Rcparts_vlass)
      y_Rcparts.append(y_Rcparts_vlass)
      errmaj_Rcparts.append(errmaj_Rcparts_vlass)
      errmin_Rcparts.append(errmin_Rcparts_vlass)
      errangle_Rcparts.append(errangle_Rcparts_vlass)

os.system('first.py '+str(ratar)+' '+str(dectar)+' '+str(4*ell_a_68*60.))
if os.path.exists(path+'/first.cat'):
   ra_Rcparts_first, dec_Rcparts_first = np.genfromtxt(path+'/first.cat',usecols=(1,2),dtype='U',unpack=True)
   ra_Rcparts_first = Angle(ra_Rcparts_first,unit="hourangle").degree
   dec_Rcparts_first = Angle(dec_Rcparts_first,unit="deg").degree
   errmaj_Rcparts_first, errmin_Rcparts_first, errangle_Rcparts_first, flux_R_first = np.genfromtxt(path+'/first.cat',usecols=(3,4,5,6),unpack=True)
   x_Rcparts_first, y_Rcparts_first = w.all_world2pix(ra_Rcparts_first, dec_Rcparts_first, 1)
   errmaj_Rcparts_first = errmaj_Rcparts_first/float(platescale)
   errmin_Rcparts_first = errmin_Rcparts_first/float(platescale)
   if np.size(ra_Rcparts_first)>1:
      cat_R.extend(np.full(np.size(ra_Rcparts_first),'first'))
      x_Rcparts.extend(x_Rcparts_first)
      y_Rcparts.extend(y_Rcparts_first)
      errmaj_Rcparts.extend(errmaj_Rcparts_first)
      errmin_Rcparts.extend(errmin_Rcparts_first)
      errangle_Rcparts.extend(errangle_Rcparts_first)
   else:
      cat_R.append('first')
      x_Rcparts.append(x_Rcparts_first)
      y_Rcparts.append(y_Rcparts_first)
      errmaj_Rcparts.append(errmaj_Rcparts_first)
      errmin_Rcparts.append(errmin_Rcparts_first)
      errangle_Rcparts.append(errangle_Rcparts_first)

os.system('MSC_4FGL.py '+str(ratar)+' '+str(dectar)+' '+str(4*ell_a_68*60.))
if os.path.exists(path+'/MSC_4FGL.cat'):
   ra_Rcparts_MSC, dec_Rcparts_MSC, errmaj_Rcparts_MSC, errmin_Rcparts_MSC, flux_R_MSC = np.genfromtxt(path+'/MSC_4FGL.cat',usecols=(1,2,3,4,6),unpack=True)
   x_Rcparts_MSC, y_Rcparts_MSC = w.all_world2pix(ra_Rcparts_MSC, dec_Rcparts_MSC, 1)
   errmaj_Rcparts_MSC = 2*errmaj_Rcparts_MSC*3600./float(platescale)
   errmin_Rcparts_MSC = 2*errmin_Rcparts_MSC*3600./float(platescale)
   if np.size(ra_Rcparts_MSC)>1:
      cat_R.extend(np.full(np.size(ra_Rcparts_MSC),'MSC'))
      x_Rcparts.extend(x_Rcparts_MSC)
      y_Rcparts.extend(y_Rcparts_MSC)
      errmaj_Rcparts.extend(errmaj_Rcparts_MSC)
      errmin_Rcparts.extend(errmin_Rcparts_MSC)
      errangle_Rcparts.extend(np.full(np.size(ra_Rcparts_MSC),np.nan))
   else:
      cat_R.append('MSC')
      x_Rcparts.append(x_Rcparts_MSC)
      y_Rcparts.append(y_Rcparts_MSC)
      errmaj_Rcparts.append(errmaj_Rcparts_MSC)
      errmin_Rcparts.append(errmin_Rcparts_MSC)
      errangle_Rcparts.append(np.nan)

os.system('FASTgpps.py '+str(ratar)+' '+str(dectar)+' '+str(4*ell_a_68*60.))
if os.path.exists(path+'/FASTgpps.cat'):
   ra_Rcparts_FAST, dec_Rcparts_FAST = np.genfromtxt(path+'/FASTgpps.cat',usecols=(2,3),dtype='U',unpack=True)
   ra_Rcparts_FAST = Angle(ra_Rcparts_FAST,unit="hourangle").degree
   dec_Rcparts_FAST = Angle(dec_Rcparts_FAST,unit="deg").degree
   errmaj_Rcparts_FAST, errmin_Rcparts_FAST, flux_R_FAST = np.genfromtxt(path+'/FASTgpps.cat',usecols=(4,5,6),unpack=True)
   x_Rcparts_FAST, y_Rcparts_FAST = w.all_world2pix(ra_Rcparts_FAST, dec_Rcparts_FAST, 1)
   errmaj_Rcparts_FAST = max(2*errmaj_Rcparts_FAST/float(platescale),4)
   errmin_Rcparts_FAST = max(2*errmin_Rcparts_FAST/float(platescale),4)
   if np.size(ra_Rcparts_FAST)>1:
      cat_R.extend(np.full(np.size(ra_Rcparts_FAST),'FAST'))
      x_Rcparts.extend(x_Rcparts_FAST)
      y_Rcparts.extend(y_Rcparts_FAST)
      errmaj_Rcparts.extend(errmaj_Rcparts_FAST)
      errmin_Rcparts.extend(errmin_Rcparts_FAST)
      errangle_Rcparts.extend(np.full(np.size(ra_Rcparts_FAST),np.nan))
   else:
      cat_R.append('FAST')
      x_Rcparts.append(x_Rcparts_FAST)
      y_Rcparts.append(y_Rcparts_FAST)
      errmaj_Rcparts.append(errmaj_Rcparts_FAST)
      errmin_Rcparts.append(errmin_Rcparts_FAST)
      errangle_Rcparts.append(np.nan)

cat_R = np.asarray(cat_R)
x_Rcparts = np.asarray(x_Rcparts)
y_Rcparts = np.asarray(y_Rcparts)
errmaj_Rcparts = np.asarray(errmaj_Rcparts)
errmin_Rcparts = np.asarray(errmin_Rcparts)
errangle_Rcparts = np.asarray(errangle_Rcparts)
flux_R = np.asarray(flux_R)

data_crop = data[max(0,round(xtar-bound*max(ell_a_95,ell_b_95))):min(np.shape(data)[1]-1,round(xtar+bound*max(ell_a_95,ell_b_95))),max(0,round(ytar-bound*max(ell_a_95,ell_b_95))):min(np.shape(data)[0]-1,round(ytar+bound*max(ell_a_95,ell_b_95)))]

#Checking nights dates
if analysis=='STELLA':
   dates = [s.split("-")[0].split("science")[1].split("A")[0] for s in images]
elif analysis=='INT':
   dates = []
   for i in range(0,len(images)):
      hdul = fits.open(path+'/rfilter/'+images[i])
      hdr = hdul[0].header
      if hdr['UTOBS'][0]=='0':
         dates.append(hdr['DATE-OBS'][:-1]+str(int(hdr['DATE-OBS'][-1])-1))
         if dates[i][-1]=='0':
            dates[i] = dates[i][:-2]+'31'
      else:
         dates.append(hdr['DATE-OBS'])
      hdul.close()
elif analysis=='LCO':
   dates = [s.split("-")[2] for s in images]
nights = []
if analysis=='INT':
   for i in range(0,np.size(dates)):
      if i==0:
         nights.append(dates[i])
      else:
         if dates[i]!=dates[i-1]:
            nights.append(dates[i])
else:
   for i in range(0,np.size(dates)):
      if i==0:
         nights.append(dates[i][:4]+"-"+dates[i][4:6]+"-"+dates[i][6:])
      else:
         if dates[i]!=dates[i-1]:
            nights.append(dates[i][:4]+"-"+dates[i][4:6]+"-"+dates[i][6:])
N_nights = np.size(nights)

#Checking nights mjd_ref
images = np.genfromtxt(path+'/rfilter/usedImages.txt',dtype='str')
mjd = [float(s.split('_')[2]) for s in images if 'median' not in s]
mjd_ref = np.zeros(N_nights)
counter = 0
for i in range(0,np.size(mjd)):
   if i!=0:
      if (mjd[i]-mjd[i-1]) > 0.8:
         mjd_ref[counter] = mjd[i-1]
         counter+=1
      if i==np.size(mjd)-1:
         mjd_ref[counter] = mjd[i]

#We get as reference epoch for the light curves the integer part of the time (in MJD) of the first image in the corresponding filter
ref_ep = np.trunc(float(np.genfromtxt(path+'/rfilter/usedImages.txt',dtype='str')[0].split('_')[2]))

if nvar!=1:
   best_periods = np.zeros(np.shape(perSources_fin)[0])
else:
   best_periods = np.zeros(1)

m, s = np.mean(data_crop), np.std(data_crop)

#For loop in the candidate variables
for i in range(0,nvar):
   if nvar!=1:
      sigma=perSources_fin[i,3]
      dm=perSources_fin[i,2]
      ra=perSources_fin[i,0]
      dec=perSources_fin[i,1]
      index=int(list_targ[i])
      #index=int(perSources_fin[i,4])
   else:
      sigma=perSources_fin[3]
      dm=perSources_fin[2]
      ra=perSources_fin[0]
      dec=perSources_fin[1]
      index=int(list_targ)
      #index=int(perSources_fin[4])

   sep_4FGL = SkyCoord(ra, dec, frame='fk5', unit='deg').separation(SkyCoord(ratar,dectar,frame='fk5',unit='deg')).degree
   if np.size(sources_3FGL)!=0:
      sep_3FGL = SkyCoord(ra, dec, frame='fk5', unit='deg').separation(SkyCoord(coord_3FGL[:,0],coord_3FGL[:,1],frame='fk5',unit='deg')).degree
      ell_a_68_3FGL = coord_3FGL[np.argmin(sep_3FGL),2]
      sep_3FGL = sep_3FGL[np.argmin(sep_3FGL)]

   #Check for association between optical candidate and X-ray counterparts within radius_X arcsec
   if os.path.exists(path+'/2sxps.cat'):
      sep_X = SkyCoord(ra, dec, frame='fk5', unit='deg').separation(SkyCoord(ra_Xcparts_Swift,dec_Xcparts_Swift,frame='fk5',unit='deg')).arcsecond
      if np.size(sep_X)!=1:
         if min(sep_X)<=radius_X:
            flux_Xcpart_Swift = flux_X_Swift[np.argmin(sep_X)]
      else:
         if sep_X<=radius_X:
            flux_Xcpart_Swift = flux_X_Swift
            
   if os.path.exists(path+'/csc2.cat'):
      sep_X = SkyCoord(ra, dec, frame='fk5', unit='deg').separation(SkyCoord(ra_Xcparts_Chandra,dec_Xcparts_Chandra,frame='fk5',unit='deg')).arcsecond
      if np.size(sep_X)!=1:
         if min(sep_X)<=radius_X:
            flux_Xcpart_Chandra = flux_X_Chandra[np.argmin(sep_X)]
      else:
         if sep_X<=radius_X:
            flux_Xcpart_Chandra = flux_X_Chandra

   if os.path.exists(path+'/XMM4d13s.cat'):
      sep_X = SkyCoord(ra, dec, frame='fk5', unit='deg').separation(SkyCoord(ra_Xcparts_XMM,dec_Xcparts_XMM,frame='fk5',unit='deg')).arcsecond
      if np.size(sep_X)!=1:
         if min(sep_X)<=radius_X:
            flux_Xcpart_XMM = flux_X_XMM[np.argmin(sep_X)]
      else:
         if sep_X<=radius_X:
            flux_Xcpart_XMM = flux_X_XMM

   if os.path.exists(path+'/eRASS1d1.cat'):
      sep_X = SkyCoord(ra, dec, frame='fk5', unit='deg').separation(SkyCoord(ra_Xcparts_eRASS,dec_Xcparts_eRASS,frame='fk5',unit='deg')).arcsecond
      if np.size(sep_X)!=1:
         if min(sep_X)<=radius_X:
            flux_Xcpart_eRASS = flux_X_eRASS[np.argmin(sep_X)]
      else:
         if sep_X<=radius_X:
            flux_Xcpart_eRASS = flux_X_eRASS

   if 'flux_Xcpart_Swift' in locals():
      flux_Xcpart = flux_Xcpart_Swift
   else:
      if 'flux_Xcpart_Chandra' in locals():
         flux_Xcpart = flux_Xcpart_Chandra
      else:
         if 'flux_Xcpart_XMM' in locals():
            flux_Xcpart = flux_Xcpart_XMM
         else:
            if 'flux_Xcpart_eRASS' in locals():
               flux_Xcpart = flux_Xcpart_eRASS
            #else:
            #   UL_X = uds.getUpperLimits(position=str(ra)+' '+str(dec), cat='LSXPS')
            #   if 'ULData' in UL_X:
            #      UL_X = uds.getUpperLimits(position=str(ra)+' '+str(dec), cat='LSXPS')['ULData']['Total_UpperLimit'][0]
            #      UL_X = counts_flux*UL_X
            #   else:
            #      del UL_X
               #print("Swift UL server out of service at the moment, try again later!")

   #Check the Kepler eclipsing binaries catalogue, and pick the closest object within a square 6 arcsec x 6 arcsec around the candidate coordinates
   query = """select * from eclipsing_binary_catalog where eb_ra < ("""+str(ra)+"""+0.00083) and eb_ra > ("""+str(ra)+"""-0.00083) and eb_dec < ("""+str(dec)+"""+0.00083) and eb_dec > ("""+str(dec)+"""-0.00083)"""
   Kepler_table = mastcasjobs.MastCasJobs(username="6813", password="Marco6813!", context="kepler").quick(query, task_name="python cone search")
   if np.size(Kepler_table)!=0:
      sep = SkyCoord(ra, dec, frame='fk5', unit='deg').separation(SkyCoord(Kepler_table['eb_ra'].value,Kepler_table['eb_dec'].value,frame='fk5',unit='deg')).arcsecond
      kic_ID = Kepler_table['kic'].value[np.argmin(sep)]
      kic_period = Kepler_table['period'].value[np.argmin(sep)]
      if 'ra_query' not in locals() and 'dec_query' not in locals():
         ra_query, dec_query = [Kepler_table['eb_ra'].value[np.argmin(sep)],Kepler_table['eb_dec'].value[np.argmin(sep)]]
   del Kepler_table

   #Check the ATLAS variable stars catalogue, and pick the closest object within a square 6 arcsec x 6 arcsec around the candidate coordinates
   query = """select * from object where ra < ("""+str(ra)+"""+0.00083) and ra > ("""+str(ra)+"""-0.00083) and dec < ("""+str(dec)+"""+0.00083) and dec > ("""+str(dec)+"""-0.00083)"""
   ATLAS_table = mastcasjobs.MastCasJobs(username="6813", password="Marco6813!", context="HLSP_ATLAS_VAR").quick(query, task_name="python cone search")
   if np.size(ATLAS_table)!=0:
      sep = SkyCoord(ra, dec, frame='fk5', unit='deg').separation(SkyCoord(ATLAS_table['ra'].value,ATLAS_table['dec'].value,frame='fk5',unit='deg')).arcsecond
      ATO_ID = ATLAS_table['ATO_ID'].value[np.argmin(sep)]
      ATO_period = ATLAS_table['fp_LSperiod'].value[np.argmin(sep)]
      if 'ra_query' not in locals() and 'dec_query' not in locals():
         ra_query, dec_query = [ATLAS_table['ra'].value[np.argmin(sep)],ATLAS_table['dec'].value[np.argmin(sep)]]
   del ATLAS_table

   #Check the SIMBAD catalogue, and pick the closest object within a radius of 3 arcsec around the candidate coordinates
   SIMBAD_table = Simbad.query_region(SkyCoord(ra,dec,frame='fk5',unit='deg'),radius=3*u.arcsec)
   #if SIMBAD_table is not None:
   #   sep = [SkyCoord(ra, dec, frame='fk5', unit='deg').separation(SkyCoord(SIMBAD_table[i]['RA']+" "+SIMBAD_table[i]['DEC'],frame='fk5',unit=(u.hourangle, u.deg))).arcsecond for i in range(0,np.size(SIMBAD_table))]
   #   SIMBAD_ID = SIMBAD_table['MAIN_ID'].value[np.argmin(sep)]
   #   if 'ra_query' not in locals() and 'dec_query' not in locals():
   #      ra_query, dec_query = [SkyCoord(SIMBAD_table[np.argmin(sep)]['RA']+" "+SIMBAD_table[np.argmin(sep)]['DEC'],frame='fk5',unit=(u.hourangle, u.deg)).ra.value,SkyCoord(SIMBAD_table[np.argmin(sep)]['RA']+" "+SIMBAD_table[np.argmin(sep)]['DEC'],frame='fk5',unit=(u.hourangle, u.deg)).dec.value]
   #   del SIMBAD_table
      
   
   if os.path.exists(path+'/rfilter/outputcats/doerPhot_V'+str(index)+'.csv') and os.path.exists(path+'/gfilter/outputcats/doerPhot_V'+str(index)+'.csv') and os.path.exists(path+'/ifilter/outputcats/doerPhot_V'+str(index)+'.csv'):
      mjd_g, m_g, sigmam_g = np.genfromtxt(path+'/gfilter/outputcats/doerPhot_V'+str(index)+'.csv', delimiter=',', usecols=(8,-2,-1),unpack=True)
      mjd_r, m_r, sigmam_r = np.genfromtxt(path+'/rfilter/outputcats/doerPhot_V'+str(index)+'.csv', delimiter=',', usecols=(8,-2,-1),unpack=True)
      mjd_i, m_i, sigmam_i = np.genfromtxt(path+'/ifilter/outputcats/doerPhot_V'+str(index)+'.csv', delimiter=',', usecols=(8,-2,-1),unpack=True)
      fwhm = np.genfromtxt(path+'/rfilter/outputcats/doerPhot_V'+str(index)+'.csv', delimiter=',', usecols=(6),unpack=True)

      #Take the light minimum in the r' band as t_0
      t_0 = mjd_r[np.argmax(m_r)]

      if per_search == 'LS':
         periods, periodogram = LombScargle_search(mjd_r, m_r, sigmam_r, P_low, P_up, n_periods_sub, 1)
         best_periods[i] = 2*periods[~np.isnan(periodogram)][np.argmax(periodogram[~np.isnan(periodogram)])]
         #best_periods[i] = periods[~np.isnan(periodogram)][np.argsort(periodogram[~np.isnan(periodogram)])][-2]
         
      elif per_search == 'PDM':
         #m_r = np.random.normal(np.mean(m_r), 0.01, len(mjd_r))
         periods, periodogram = PDM_search(mjd_r, m_r, P_low, P_up, 0., n_periods_sub, n_bin_sub, min_bin)
         best_periods[i] = periods[~np.isnan(periodogram)][np.argmin(periodogram[~np.isnan(periodogram)])]
         #best_periods[i] = periods[~np.isnan(periodogram)][np.argsort(periodogram[~np.isnan(periodogram)])][3]
      else:
         periods, periodogram = PDM_trisearch(mjd_g, m_g, mjd_r, m_r, mjd_i, m_i, P_low, P_up, 0., n_periods_sub, n_bin_sub, min_bin)
         best_periods[i] = periods[~np.isnan(periodogram)][np.argmin(periodogram[~np.isnan(periodogram)])]
      
      fwhmavg = np.mean(fwhm)

      nightdim = np.zeros(N_nights,int)
      for j in range(0,(np.size(mjd_r))):
         if j!=0:
            if (mjd_r[j]-mjd_r[j-1]) > 0.8:
               n_night = np.argmin(np.fabs(mjd_r[j-1]-mjd_ref))
               nightdim[n_night] = j-nightdim[0]
               for k in range(1,n_night):
                  nightdim[n_night] = nightdim[n_night]-nightdim[k]
            if j==np.size(mjd_r)-1:
               n_night = np.argmin(np.fabs(mjd_r[j]-mjd_ref))
               nightdim[n_night] = j+1-nightdim[0]
               for k in range(1,n_night):
                  nightdim[n_night] = nightdim[n_night]-nightdim[k]      
   
      mjd_single = []
      m_single = []
      sigmam_single = []
      n_lc = 0
      for j in range(0,N_nights):
         mjd_single.append(mjd_r[np.sum(nightdim[0:j]):np.sum(nightdim[0:j])+nightdim[j]])
         m_single.append(m_r[np.sum(nightdim[0:j]):np.sum(nightdim[0:j])+nightdim[j]])
         sigmam_single.append(sigmam_r[np.sum(nightdim[0:j]):np.sum(nightdim[0:j])+nightdim[j]])
         if nightdim[j]!=0:
            n_lc = n_lc+1

      figure=plt.figure(figsize=(figxsize,figysize), dpi=150)
      if np.size(sources_3FGL)!=0:
         if 'flux_Xcpart' in locals():
            figure.suptitle(r"Fermi source: "+FoVname+", Instrument: "+telescope+", Source detected in "+str(n_lc)+" nights out of "+str(N_nights)+"\n\nExposure time: %.1f"%deltat+" s, Plate scale: "+str(platescale)+" arcsec/pix, Extraction radius: "+str(fext)+"*FWHM\n\nPeriodic variable #%d"%(i+1)+", $\sigma$: %.2f"%sigma+" mag, $\Delta$m: %.2f"%dm+" mag, Range: [%.2f"%P_low+",%.1f"%min(P_up,mjd_r[-1]-mjd_r[0])+"] d, Steps: %d"%n_periods_sub+", Period folding with $T_{0}=$%.4f"%t_0+" MJD\n\n%.1f"%(sep_4FGL/ell_a_68)+r"$\sigma$ far from the 4FGL source and %.1f"%(sep_3FGL/ell_a_68_3FGL)+r"$\sigma$ far from the 3FGL source, X-ray flux: "+str(np.round(float(str(flux_Xcpart).split('e')[0]),1))+r"$\cdot 10^{"+str(flux_Xcpart).split('e')[1]+"}$ erg s$^{-1}$ cm$^{-2}$", fontsize=tsize)
         else:
            if 'UL_X' in locals():
               figure.suptitle(r"Fermi source: "+FoVname+", Instrument: "+telescope+", Source detected in "+str(n_lc)+" nights out of "+str(N_nights)+"\n\nExposure time: %.1f"%deltat+" s, Plate scale: "+str(platescale)+" arcsec/pix, Extraction radius: "+str(fext)+"*FWHM\n\nPeriodic variable #%d"%(i+1)+", $\sigma$: %.2f"%sigma+" mag, $\Delta$m: %.2f"%dm+" mag, Range: [%.2f"%P_low+",%.1f"%min(P_up,mjd_r[-1]-mjd_r[0])+"] d, Steps: %d"%n_periods_sub+", Period folding with $T_{0}=$%.4f"%t_0+" MJD\n\n%.1f"%(sep_4FGL/ell_a_68)+r"$\sigma$ far from the 4FGL source and %.1f"%(sep_3FGL/ell_a_68_3FGL)+r"$\sigma$ far from the 3FGL source, X-ray flux 3$\sigma$ UL: "+str(np.round(float(str(UL_X).split('e')[0]),1))+r"$\cdot 10^{"+str(UL_X).split('e')[1]+"}$ erg s$^{-1}$ cm$^{-2}$", fontsize=tsize)
            else:
               figure.suptitle(r"Fermi source: "+FoVname+", Instrument: "+telescope+", Source detected in "+str(n_lc)+" nights out of "+str(N_nights)+"\n\nExposure time: %.1f"%deltat+" s, Plate scale: "+str(platescale)+" arcsec/pix, Extraction radius: "+str(fext)+"*FWHM\n\nPeriodic variable #%d"%(i+1)+", $\sigma$: %.2f"%sigma+" mag, $\Delta$m: %.2f"%dm+" mag, Range: [%.2f"%P_low+",%.1f"%min(P_up,mjd_r[-1]-mjd_r[0])+"] d, Steps: %d"%n_periods_sub+", Period folding with $T_{0}=$%.4f"%t_0+" MJD\n\n%.1f"%(sep_4FGL/ell_a_68)+r"$\sigma$ far from the 4FGL source and %.1f"%(sep_3FGL/ell_a_68_3FGL)+r"$\sigma$ far from the 3FGL source", fontsize=tsize)
      else:
         if 'flux_Xcpart' in locals():
            figure.suptitle(r"Fermi source: "+FoVname+", Instrument: "+telescope+", Source detected in "+str(n_lc)+" nights out of "+str(N_nights)+"\n\nExposure time: %.1f"%deltat+" s, Plate scale: "+str(platescale)+" arcsec/pix, Extraction radius: "+str(fext)+"*FWHM\n\nPeriodic variable #%d"%(i+1)+", $\sigma$: %.2f"%sigma+" mag, $\Delta$m: %.2f"%dm+" mag, Range: [%.2f"%P_low+",%.1f"%min(P_up,mjd_r[-1]-mjd_r[0])+"] d, Steps: %d"%n_periods_sub+", Period folding with $T_{0}=$%.4f"%t_0+" MJD\n\n%.1f"%(sep_4FGL/ell_a_68)+r"$\sigma$ far from the 4FGL source, X-ray flux: "+str(np.round(float(str(flux_Xcpart).split('e')[0]),1))+r"$\cdot 10^{"+str(flux_Xcpart).split('e')[1]+"}$ erg s$^{-1}$ cm$^{-2}$", fontsize=tsize)
         else:
            if 'UL_X' in locals():
               figure.suptitle(r"Fermi source: "+FoVname+", Instrument: "+telescope+", Source detected in "+str(n_lc)+" nights out of "+str(N_nights)+"\n\nExposure time: %.1f"%deltat+" s, Plate scale: "+str(platescale)+" arcsec/pix, Extraction radius: "+str(fext)+"*FWHM\n\nPeriodic variable #%d"%(i+1)+", $\sigma$: %.2f"%sigma+" mag, $\Delta$m: %.2f"%dm+" mag, Range: [%.2f"%P_low+",%.1f"%min(P_up,mjd_r[-1]-mjd_r[0])+"] d, Steps: %d"%n_periods_sub+", Period folding with $T_{0}=$%.4f"%t_0+" MJD\n\n%.1f"%(sep_4FGL/ell_a_68)+r"$\sigma$ far from the 4FGL source, X-ray flux 3$\sigma$ UL: "+str(np.round(float(str(UL_X).split('e')[0]),1))+r"$\cdot 10^{"+str(UL_X).split('e')[1]+"}$ erg s$^{-1}$ cm$^{-2}$", fontsize=tsize)
            else:
               figure.suptitle(r"Fermi source: "+FoVname+", Instrument: "+telescope+", Source detected in "+str(n_lc)+" nights out of "+str(N_nights)+"\n\nExposure time: %.1f"%deltat+" s, Plate scale: "+str(platescale)+" arcsec/pix, Extraction radius: "+str(fext)+"*FWHM\n\nPeriodic variable #%d"%(i+1)+", $\sigma$: %.2f"%sigma+" mag, $\Delta$m: %.2f"%dm+" mag, Range: [%.2f"%P_low+",%.1f"%min(P_up,mjd_r[-1]-mjd_r[0])+"] d, Steps: %d"%n_periods_sub+", Period folding with $T_{0}=$%.4f"%t_0+" MJD\n\n%.1f"%(sep_4FGL/ell_a_68)+r"$\sigma$ far from the 4FGL source", fontsize=tsize)

      counter = 0
      for j in range(0,N_nights):
         if nightdim[j]!=0:
            plt.subplot(int(np.trunc(5+n_lc/2)), 2, counter+1)
            plt.plot(mjd_single[j]-ref_ep, m_single[j], 'bo', linestyle='None')
            plt.errorbar(mjd_single[j]-ref_ep, m_single[j], sigmam_single[j], linestyle='None')
            plt.xticks(fontsize=axticksize)
            plt.yticks(fontsize=axticksize)
            plt.gca().invert_yaxis()
            plt.title("Lightcurve in the r band of "+nights[j],fontsize=subsize)
            plt.xlabel(r"Time (MJD-%d)"%ref_ep,fontsize=axlabelsize)
            plt.ylabel(r"Apparent magnitude",fontsize=axlabelsize)
            counter = counter+1

      plt.subplot(int(np.trunc(5+n_lc/2)), 2, n_lc+1)
      plt.xlabel(r'$\Delta$m (mag)',fontsize=axlabelsize)
      plt.ylabel(r'$\sigma$ (mag)',fontsize=axlabelsize)
      plt.title("Dispersion vs differential magnitude in the r band",fontsize=subsize)
      plt.yscale("log")
      plt.plot(sourcesr[:,2],sourcesr[:,3],'ob',markerfacecolor='0.8',markeredgecolor='black',markersize=markersize,alpha=alpha)
      plt.plot(varSources_fin1[:,2],varSources_fin1[:,3],'ob',markerfacecolor='red',markeredgecolor='black',markersize=markersize,alpha=alpha)
      #plt.plot(varSources_fin2[:,2],varSources_fin2[:,3],'ob',markerfacecolor='green',markeredgecolor='black',markersize=markersize,alpha=alpha)
      if nlist==1:
         if nvar!=1:
            plt.plot(perSources_fin[i,2],perSources_fin[i,3],'ob',markerfacecolor='red',markeredgecolor='black',markersize=markersize_main)
         else:
            plt.plot(perSources_fin[2],perSources_fin[3],'ob',markerfacecolor='red',markeredgecolor='black',markersize=markersize_main)
      else:
         if nvar!=1:
            plt.plot(perSources_fin[i,2],perSources_fin[i,3],'ob',markerfacecolor='green',markeredgecolor='black',markersize=markersize_main)
         else:
            plt.plot(perSources_fin[2],perSources_fin[3],'ob',markerfacecolor='green',markeredgecolor='black',markersize=markersize_main)
      plt.xticks(fontsize=axticksize)
      plt.yticks(fontsize=axticksize)

      mjd_gr = []
      g_r = []
      sigma_gr = []
      mjd_ri = []
      r_i = []
      sigma_ri = []
      for j in range(0,np.size(mjd_g)):
         for k in range(0,np.size(mjd_r)):
            if np.size(mjd_g)!=1 and np.size(mjd_r)!=1:
               if (mjd_r[k]-mjd_g[j]) < t_space*deltat/86400. and (mjd_r[k]-mjd_g[j])>0:
                  mjd_gr.append((mjd_g[j]+mjd_r[k])/2.)
                  g_r.append(m_g[j]-m_r[k])
                  sigma_gr.append(np.sqrt(sigmam_g[j]**2+sigmam_r[k]**2))
      for j in range(0,np.size(mjd_r)):
         for k in range(0,np.size(mjd_i)):
            if np.size(mjd_r)!=1 and np.size(mjd_i)!=1:
               if (mjd_i[k]-mjd_r[j]) < t_space*deltat/86400. and (mjd_i[k]-mjd_r[j])>0:
                  mjd_ri.append((mjd_r[j]+mjd_i[k])/2.)
                  r_i.append(m_r[j]-m_i[k])
                  sigma_ri.append(np.sqrt(sigmam_r[j]**2+sigmam_i[k]**2))
      mjd_gr = np.asarray(mjd_gr)
      g_r = np.asarray(g_r)
      sigma_gr = np.asarray(sigma_gr)
      mjd_ri = np.asarray(mjd_ri)
      r_i = np.asarray(r_i)
      sigma_ri = np.asarray(sigma_ri)

      phase_g = ((mjd_g-t_0)/best_periods[i]) % 1
      #sort_index = np.asarray(phase_g.argsort())
      #phase_g = phase_g[sort_index]
      #m_g = m_g[sort_index]
      #sigmam_g = sigmam_g[sort_index]
      phase_r = ((mjd_r-t_0)/best_periods[i]) % 1
      #sort_index = np.asarray(phase_r.argsort())
      #phase_r = phase_r[sort_index]
      #m_r = m_r[sort_index]
      #sigmam_r = sigmam_r[sort_index]
      phase_i = ((mjd_i-t_0)/best_periods[i]) % 1
      #sort_index = np.asarray(phase_i.argsort())
      #phase_i = phase_i[sort_index]
      #m_i = m_i[sort_index]
      #sigmam_i = sigmam_i[sort_index]
      if np.size(mjd_g)!=1 and np.size(mjd_r)!=1:
         phase_gr = ((mjd_gr-t_0)/best_periods[i]) % 1
      #sort_index = np.asarray(phase_gr.argsort())
      #phase_gr = phase_gr[sort_index]
      #g_r = g_r[sort_index]
      #sigma_gr = sigma_gr[sort_index]
      if np.size(mjd_r)!=1 and np.size(mjd_i)!=1:
         phase_ri = ((mjd_ri-t_0)/best_periods[i]) % 1
      #sort_index = np.asarray(phase_ri.argsort())
      #phase_ri = phase_ri[sort_index]
      #r_i = r_i[sort_index]
      #sigma_ri = sigma_ri[sort_index]

      plt.subplot(int(np.trunc(5+n_lc/2)), 2, n_lc+2)
      plt.plot(phase_r, m_r, 'bo', linestyle='None')
      plt.plot(phase_r+1, m_r, 'ro', linestyle='None')
      plt.errorbar(phase_r, m_r, sigmam_r, linestyle='None')
      plt.errorbar(phase_r+1, m_r, sigmam_r, linestyle='None')
      plt.xticks(fontsize=axticksize)
      plt.yticks(fontsize=axticksize)
      plt.gca().invert_yaxis()
      plt.title(r"Folded lightcurve in the r band: $P_{\mathrm{best}}$="+"{0} d".format(best_periods[i]),fontsize=subsize)
      plt.xlabel(r"Phase ($\phi$)",fontsize=axlabelsize)
      plt.ylabel(r"Apparent magnitude",fontsize=axlabelsize)

      plt.subplot(int(np.trunc(5+n_lc/2)), 2, n_lc+3)
      plt.plot(periods,periodogram)
      if not np.isnan(max(periodogram)) and not np.isinf(max(periodogram)) and not np.isnan(min(periodogram)) and not np.isinf(min(periodogram)) and not np.isnan(best_periods[i]) and not np.isinf(best_periods[i]):
         plt.axvline(x=best_periods[i], color='r', ls='--')
      plt.xscale("log")
      plt.xticks(fontsize=axticksize)
      plt.yticks(fontsize=axticksize)
      if per_search == 'LS':
         plt.title("Lomb-Scargle data periodogram in the r band",fontsize=subsize)
         plt.ylabel(r"Power",fontsize=axlabelsize)
         if period_FAPl == 'yes':
            tic_1 = time.perf_counter()
            FAP_periods, FAP_levels = LombScargle_FAP_level(mjd_r, m_r, sigmam_r, P_low, P_up, n_periods, 1, FAPl_factor, conf_level)
            toc_1 = time.perf_counter()
            print(f"FAP levels computed in {toc_1 - tic_1:0.4f} seconds")
            plt.plot(FAP_periods,FAP_levels,color='darkorange',ls='--')
            plt.text(plt.xlim()[0]*1.4, np.mean(FAP_levels)*1.2, r'0.1% FAP', color='darkorange', ha='center',fontsize=axticksize)
      elif per_search == 'PDM':
         plt.title("PDM data periodogram in the r band",fontsize=subsize)
         plt.ylabel(r"Normalized variance",fontsize=axlabelsize)
         if period_FAPl == 'yes':
            tic_1 = time.perf_counter()
            FAP_periods, FAP_levels = PDM_FAP_level(mjd_r, m_r, P_low, P_up, 0., n_periods_sub, n_bin_sub, min_bin, FAPl_factor, conf_level)
            toc_1 = time.perf_counter()
            print(f"FAP levels computed in {toc_1 - tic_1:0.4f} seconds")
            plt.plot(FAP_periods,FAP_levels,color='darkorange',ls='--')
            plt.text(plt.xlim()[0]*1.4, np.mean(FAP_levels)*1.02, r'0.1% FAP', color='darkorange', ha='center',fontsize=axticksize)
      else:
         plt.title("PDM multi-band data periodogram",fontsize=subsize)
         plt.ylabel(r"Sum of normalized variances",fontsize=axlabelsize)
      plt.xlabel(r"Trial period (days)",fontsize=axlabelsize)

      plt.subplot(int(np.trunc(5+n_lc/2)), 2, n_lc+4)
      plt.plot(phase_g, m_g, 'bo', linestyle='None')
      plt.plot(phase_g+1, m_g, 'ro', linestyle='None')
      plt.errorbar(phase_g, m_g, sigmam_g, linestyle='None')
      plt.errorbar(phase_g+1, m_g, sigmam_g, linestyle='None')
      plt.xticks(fontsize=axticksize)
      plt.yticks(fontsize=axticksize)
      plt.gca().invert_yaxis()
      plt.title(r"Folded lightcurve in the g band: $P_{\mathrm{best}}$="+"{0} d".format(best_periods[i]),fontsize=subsize)
      plt.xlabel(r"Phase ($\phi$)",fontsize=axlabelsize)
      plt.ylabel(r"Apparent magnitude",fontsize=axlabelsize)

      plt.subplot(int(np.trunc(5+n_lc/2)), 2, n_lc+5)
      plt.plot(phase_i, m_i, 'bo', linestyle='None')
      plt.plot(phase_i+1, m_i, 'ro', linestyle='None')
      plt.errorbar(phase_i, m_i, sigmam_i, linestyle='None')
      plt.errorbar(phase_i+1, m_i, sigmam_i, linestyle='None')
      plt.xticks(fontsize=axticksize)
      plt.yticks(fontsize=axticksize)
      plt.gca().invert_yaxis()
      plt.title(r"Folded lightcurve in the i band: $P_{\mathrm{best}}$="+"{0} d".format(best_periods[i]),fontsize=subsize)
      plt.xlabel(r"Phase ($\phi$)",fontsize=axlabelsize)
      plt.ylabel(r"Apparent magnitude",fontsize=axlabelsize)

      if np.size(mjd_g)!=1 and np.size(mjd_r)!=1:
         plt.subplot(int(np.trunc(5+n_lc/2)), 2, n_lc+6)
         plt.plot(phase_gr, g_r, 'bo', linestyle='None')
         plt.plot(phase_gr+1, g_r, 'ro', linestyle='None')
         plt.errorbar(phase_gr, g_r, sigma_gr, linestyle='None')
         plt.errorbar(phase_gr+1, g_r, sigma_gr, linestyle='None')
         plt.xticks(fontsize=axticksize)
         plt.yticks(fontsize=axticksize)
         plt.gca().invert_yaxis()
         plt.title("(g-r) folded colour curve",fontsize=subsize)
         plt.xlabel(r"Phase ($\phi$)",fontsize=axlabelsize)
         plt.ylabel(r"Colour",fontsize=axlabelsize)
      
      if np.size(mjd_r)!=1 and np.size(mjd_i)!=1:
         plt.subplot(int(np.trunc(5+n_lc/2)), 2, n_lc+7)
         plt.plot(phase_ri, r_i, 'bo', linestyle='None')
         plt.plot(phase_ri+1, r_i, 'ro', linestyle='None')
         plt.errorbar(phase_ri, r_i, sigma_ri, linestyle='None')
         plt.errorbar(phase_ri+1, r_i, sigma_ri, linestyle='None')
         plt.xticks(fontsize=axticksize)
         plt.yticks(fontsize=axticksize)
         plt.gca().invert_yaxis()
         plt.title("(r-i) folded colour curve",fontsize=subsize)
         plt.xlabel(r"Phase ($\phi$)",fontsize=axlabelsize)
         plt.ylabel(r"Colour",fontsize=axlabelsize)

      ax1 = plt.subplot(int(np.trunc(5+n_lc/2)), 2, n_lc+8)
      plt.title(r'Median r FoV bkg subtracted in the '+str(2*bound)+'$\sigma$ square around the 4FGL coordinates',fontsize=subsize)
      plt.xlim(max(0,xtar-bound*max(ell_a_95,ell_b_95)), min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)))
      plt.ylim(max(0,ytar-bound*max(ell_a_95,ell_b_95)), min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)))
      plt.imshow(data, interpolation='nearest', cmap='gray', vmin=m-colscale*s, vmax=m+colscale*s, origin='lower')
      plt.text(xtar,ytar,FoVname,horizontalalignment='center',verticalalignment='center',color='cyan',fontsize=subsize,alpha=0.8)
      target = Ellipse(xy=(xtar, ytar), width=2*ell_a_95, height=2*ell_b_95, angle=ell_theta_95)
      target.set_facecolor('none')
      target.set_edgecolor('cyan')
      ax1.add_artist(target)
      if np.size(sources_3FGL)!=0:
         for j in range(0,np.shape(sources_3FGL)[0]):
            if x_3FGL[j]>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_3FGL[j]<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_3FGL[j]>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_3FGL[j]<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
               plt.text(x_3FGL[j],y_3FGL[j],sources_3FGL[j].replace("L ","L"),horizontalalignment='center',verticalalignment='center',color='yellow',fontsize=subsize,alpha=0.8)
               a = Ellipse(xy=(x_3FGL[j], y_3FGL[j]), width=2*3600./float(platescale)*coord_3FGL[j,3], height=2*3600./float(platescale)*coord_3FGL[j,4], angle=coord_3FGL[j,5])
               a.set_facecolor('none')
               a.set_edgecolor('yellow')
               ax1.add_artist(a)
      #if np.size(sources_4FGLe)!=0:
      #   for j in range(0,np.shape(sources_4FGLe)[0]):
      #      if x_4FGLe[j]>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_4FGLe[j]<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_4FGLe[j]>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_4FGLe[j]<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
      #         plt.text(x_4FGLe[j],y_4FGLe[j],sources_4FGLe[j].replace("L ","L"),horizontalalignment='center',verticalalignment='center',color='cyan',fontsize=subsize,alpha=0.8)
      #         a = Ellipse(xy=(x_4FGLe[j], y_4FGLe[j]), width=2*3600./float(platescale)*np.asarray(coord_4FGLe)[j,2], height=2*3600./float(platescale)*np.asarray(coord_4FGLe)[j,3], angle=np.asarray(coord_4FGLe)[j,4])
      #         a.set_facecolor('none')
      #         a.set_edgecolor('cyan')
      #         ax1.add_artist(a)
         
      x, y = w.all_world2pix(ra, dec, 1)
      plt.text(x,y+50,str(i+1),horizontalalignment='center',verticalalignment='center',color='r',fontsize=subsize)
      a = Circle(xy=(x, y), radius=fext*fwhmavg)
      a.set_linewidth(0.5)
      a.set_facecolor('none')
      a.set_edgecolor('red')
      ax1.add_artist(a)
      #if hdr['INSTRUME']=='WIFSIP V2.0':
      #   if np.size(obs_Rcparts)!=0:
      #      for j in range(0,np.size(obs_Rcparts)):
      #         if x_Rcparts[j]>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Rcparts[j]<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Rcparts[j]>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Rcparts[j]<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
      #            a = Circle(xy=(x_Rcparts[j], y_Rcparts[j]), radius=8)
      #            a.set_linewidth(0.5)
      #            a.set_facecolor('none')
      #            plt.text(x_Rcparts[j],y_Rcparts[j]+50,obs_Rcparts[j],horizontalalignment='center',verticalalignment='center',color='limegreen',fontsize=1.2*subsize,alpha=0.8)
      #            a.set_edgecolor('limegreen')
      #            ax1.add_artist(a)
      if os.path.exists(path+'/2sxps.cat'):
         if np.size(x_Xcparts_Swift)!=1:
            for j in range(0,np.size(x_Xcparts_Swift)):
               if x_Xcparts_Swift[j]>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Xcparts_Swift[j]<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Xcparts_Swift[j]>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Xcparts_Swift[j]<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
                  a = Circle(xy=(x_Xcparts_Swift[j], y_Xcparts_Swift[j]), radius=err_Xcparts_Swift[j])
                  a.set_linewidth(0.5)
                  a.set_facecolor('none')
                  plt.text(x_Xcparts_Swift[j],y_Xcparts_Swift[j]+50,'swift',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.2*subsize,alpha=1.0)
                  a.set_edgecolor('fuchsia')
                  ax1.add_artist(a)
         else:
            if x_Xcparts_Swift>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Xcparts_Swift<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Xcparts_Swift>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Xcparts_Swift<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
               a = Circle(xy=(x_Xcparts_Swift, y_Xcparts_Swift), radius=err_Xcparts_Swift)
               a.set_linewidth(0.5)
               a.set_facecolor('none')
               plt.text(x_Xcparts_Swift,y_Xcparts_Swift+50,'swift',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.2*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax1.add_artist(a)
      if os.path.exists(path+'/csc2.cat'):
         if np.size(x_Xcparts_Chandra)!=1:
            for j in range(0,np.size(x_Xcparts_Chandra)):
               if x_Xcparts_Chandra[j]>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Xcparts_Chandra[j]<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Xcparts_Chandra[j]>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Xcparts_Chandra[j]<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
                  a = Ellipse(xy=(x_Xcparts_Chandra[j], y_Xcparts_Chandra[j]), width=2*ell_a_Chandra[j], height=2*ell_b_Chandra[j], angle=ell_theta_Chandra[j])
                  a.set_linewidth(0.5)
                  a.set_facecolor('none')
                  plt.text(x_Xcparts_Chandra[j],y_Xcparts_Chandra[j]+50,'chandra',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.2*subsize,alpha=1.0)
                  a.set_edgecolor('fuchsia')
                  ax1.add_artist(a)
         else:
            if x_Xcparts_Chandra>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Xcparts_Chandra<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Xcparts_Chandra>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Xcparts_Chandra<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
               a = Ellipse(xy=(x_Xcparts_Chandra, y_Xcparts_Chandra), width=2*ell_a_Chandra, height=2*ell_b_Chandra, angle=ell_theta_Chandra)
               a.set_linewidth(0.5)
               a.set_facecolor('none')
               plt.text(x_Xcparts_Chandra,y_Xcparts_Chandra+50,'chandra',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.2*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax1.add_artist(a)
      if os.path.exists(path+'/XMM4d13s.cat'):
         if np.size(x_Xcparts_XMM)!=1:
            for j in range(0,np.size(x_Xcparts_XMM)):
               if x_Xcparts_XMM[j]>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Xcparts_XMM[j]<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Xcparts_XMM[j]>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Xcparts_XMM[j]<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
                  a = Circle(xy=(x_Xcparts_XMM[j], y_Xcparts_XMM[j]), radius=err_Xcparts_XMM[j])
                  a.set_linewidth(0.5)
                  a.set_facecolor('none')
                  plt.text(x_Xcparts_XMM[j],y_Xcparts_XMM[j]+50,'XMM',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.2*subsize,alpha=1.0)
                  a.set_edgecolor('fuchsia')
                  ax1.add_artist(a)
         else:
            if x_Xcparts_XMM>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Xcparts_XMM<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Xcparts_XMM>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Xcparts_XMM<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
               a = Circle(xy=(x_Xcparts_XMM, y_Xcparts_XMM), radius=err_Xcparts_XMM)
               a.set_linewidth(0.5)
               a.set_facecolor('none')
               plt.text(x_Xcparts_XMM,y_Xcparts_XMM+50,'XMM',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.2*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax1.add_artist(a)
      if os.path.exists(path+'/eRASS1d1.cat'):
         if np.size(x_Xcparts_eRASS)!=1:
            for j in range(0,np.size(x_Xcparts_eRASS)):
               if x_Xcparts_eRASS[j]>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Xcparts_eRASS[j]<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Xcparts_eRASS[j]>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Xcparts_eRASS[j]<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
                  a = Circle(xy=(x_Xcparts_eRASS[j], y_Xcparts_eRASS[j]), radius=err_Xcparts_eRASS[j])
                  a.set_linewidth(0.5)
                  a.set_facecolor('none')
                  plt.text(x_Xcparts_eRASS[j],y_Xcparts_eRASS[j]+50,'eRASS',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.2*subsize,alpha=1.0)
                  a.set_edgecolor('fuchsia')
                  ax1.add_artist(a)
         else:
            if x_Xcparts_eRASS>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Xcparts_eRASS<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Xcparts_eRASS>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Xcparts_eRASS<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
               a = Circle(xy=(x_Xcparts_eRASS, y_Xcparts_eRASS), radius=err_Xcparts_eRASS)
               a.set_linewidth(0.5)
               a.set_facecolor('none')
               plt.text(x_Xcparts_eRASS,y_Xcparts_eRASS+50,'eRASS',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.2*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax1.add_artist(a)
      if np.size(x_Rcparts) > 0:
         if np.size(x_Rcparts) > 1:
            for j in range(0,np.size(x_Rcparts)):
               if x_Rcparts[j]>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Rcparts[j]<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Rcparts[j]>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Rcparts[j]<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
                  if not np.isnan(errmaj_Rcparts[j]) and not np.isnan(errmin_Rcparts[j]) and not np.isnan(errangle_Rcparts[j]):
                     a = Ellipse(xy=(x_Rcparts[j], y_Rcparts[j]), width=2*errmaj_Rcparts[j], height=2*errmin_Rcparts[j], angle=errangle_Rcparts[j])
                  elif not np.isnan(errmaj_Rcparts[j]) and not np.isnan(errmin_Rcparts[j]) and np.isnan(errangle_Rcparts[j]):
                     a = Circle(xy=(x_Rcparts[j], y_Rcparts[j]),radius=errmaj_Rcparts[j])
                  else:
                     a = Circle(xy=(x_Rcparts[j], y_Rcparts[j]),radius=2.)
                  a.set_linewidth(0.5)
                  a.set_facecolor('none')
                  plt.text(x_Rcparts[j],y_Rcparts[j]+50,cat_R[j],horizontalalignment='center',verticalalignment='center',color='limegreen',fontsize=1.2*subsize,alpha=1.0)
                  a.set_edgecolor('limegreen')
                  ax1.add_artist(a)
         else:
            if x_Rcparts>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Rcparts<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Rcparts>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Rcparts<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
               if not np.isnan(errmaj_Rcparts) and not np.isnan(errmin_Rcparts) and not np.isnan(errangle_Rcparts):
                  a = Ellipse(xy=(x_Rcparts, y_Rcparts), width=2*errmaj_Rcparts, height=2*errmin_Rcparts, angle=errangle_Rcparts)
               elif not np.isnan(errmaj_Rcparts) and not np.isnan(errmin_Rcparts) and np.isnan(errangle_Rcparts):
                  a = Circle(xy=(x_Rcparts, y_Rcparts),radius=errmaj_Rcparts)
               else:
                  a = Circle(xy=(x_Rcparts, y_Rcparts),radius=2.)
               a.set_linewidth(0.5)
               a.set_facecolor('none')
               plt.text(x_Rcparts,y_Rcparts+50,cat_R,horizontalalignment='center',verticalalignment='center',color='limegreen',fontsize=1.2*subsize,alpha=1.0)
               a.set_edgecolor('limegreen')
               ax1.add_artist(a)

      if 'ra_query' in locals() and 'dec_query' in locals():
         x_query, y_query = w.all_world2pix(ra_query, dec_query, 1)
         plt.text(x_query,y_query+50,"query",horizontalalignment='center',verticalalignment='center',color='darkorange',fontsize=1.2*subsize,alpha=1.0)
         a = Circle(xy=(x_query, y_query), radius=8)
         a.set_linewidth(0.5)
         a.set_facecolor('none')
         a.set_edgecolor('darkorange')
         ax1.add_artist(a)
         

      ax2 = plt.subplot(int(np.trunc(5+n_lc/2)), 2, n_lc+9)
      plt.title('Zoom on the variable source: (RA, DEC) = (%.3f'%ra+', %.3f'%dec+') deg',fontsize=subsize)
      plt.xlim(x-zoom,x+zoom)
      plt.ylim(y-zoom,y+zoom)
      plt.imshow(data, interpolation='nearest', cmap='gray', vmin=m-colscale*s, vmax=m+colscale*s, origin='lower')
      rxnew = [rx[j] for j in range(0,len(rx)) if rx[j]>=(x-zoom) and rx[j]<=(x+zoom) and ry[j]>=(y-zoom) and ry[j]<=(y+zoom)]
      rynew = [ry[j] for j in range(0,len(rx)) if rx[j]>=(x-zoom) and rx[j]<=(x+zoom) and ry[j]>=(y-zoom) and ry[j]<=(y+zoom)]
      for j in range(len(rxnew)):
         if not(np.round(x,ndig)==np.round(rxnew[j],ndig) and np.round(y,ndig)==np.round(rynew[j],ndig)):
            a = Circle(xy=(rxnew[j], rynew[j]), radius=fext*fwhmavg)
            a.set_linewidth(1.0)
            a.set_facecolor('none')
            a.set_edgecolor('blue')
            ax2.add_artist(a)
      plt.text(x,y+20,str(i+1),horizontalalignment='center',verticalalignment='center',color='r',fontsize=2*subsize)
      a = Circle(xy=(x, y), radius=fext*fwhmavg)
      a.set_linewidth(1.0)
      a.set_facecolor('none')
      a.set_edgecolor('red')
      ax2.add_artist(a)
      #if hdr['INSTRUME']=='WIFSIP V2.0':
      #   if np.size(obs_Rcparts)!=0:
      #      xnew_Rcparts = [x_Rcparts[j] for j in range(0,len(x_Rcparts)) if x_Rcparts[j]>=(x-zoom) and x_Rcparts[j]<=(x+zoom) and y_Rcparts[j]>=(y-zoom) and y_Rcparts[j]<=(y+zoom)]
      #      ynew_Rcparts = [y_Rcparts[j] for j in range(0,len(x_Rcparts)) if x_Rcparts[j]>=(x-zoom) and x_Rcparts[j]<=(x+zoom) and y_Rcparts[j]>=(y-zoom) and y_Rcparts[j]<=(y+zoom)]
      #      obsnew_Rcparts = [obs_Rcparts[j] for j in range(0,len(x_Rcparts)) if x_Rcparts[j]>=(x-zoom) and x_Rcparts[j]<=(x+zoom) and y_Rcparts[j]>=(y-zoom) and y_Rcparts[j]<=(y+zoom)]
      #      for j in range(0,np.size(xnew_Rcparts)):
      #         a = Circle(xy=(xnew_Rcparts[j], ynew_Rcparts[j]), radius=8)
      #         a.set_linewidth(1.0)
      #         a.set_facecolor('none')
      #         plt.text(xnew_Rcparts[j],ynew_Rcparts[j],obsnew_Rcparts[j],horizontalalignment='center',verticalalignment='center',color='limegreen',fontsize=1.5*subsize,alpha=1.0)
      #         a.set_edgecolor('limegreen')
      #         ax2.add_artist(a)
      if os.path.exists(path+'/2sxps.cat'):
         if np.size(x_Xcparts_Swift)!=1:
            xnew_Xcparts = [x_Xcparts_Swift[j] for j in range(0,len(x_Xcparts_Swift)) if x_Xcparts_Swift[j]>=(x-zoom) and x_Xcparts_Swift[j]<=(x+zoom) and y_Xcparts_Swift[j]>=(y-zoom) and y_Xcparts_Swift[j]<=(y+zoom)]
            ynew_Xcparts = [y_Xcparts_Swift[j] for j in range(0,len(x_Xcparts_Swift)) if x_Xcparts_Swift[j]>=(x-zoom) and x_Xcparts_Swift[j]<=(x+zoom) and y_Xcparts_Swift[j]>=(y-zoom) and y_Xcparts_Swift[j]<=(y+zoom)]
            errnew_Xcparts = [err_Xcparts_Swift[j] for j in range(0,len(x_Xcparts_Swift)) if x_Xcparts_Swift[j]>=(x-zoom) and x_Xcparts_Swift[j]<=(x+zoom) and y_Xcparts_Swift[j]>=(y-zoom) and y_Xcparts_Swift[j]<=(y+zoom)]
            for j in range(0,np.size(xnew_Xcparts)):
               a = Circle(xy=(xnew_Xcparts[j], ynew_Xcparts[j]), radius=errnew_Xcparts[j])
               a.set_linewidth(1.0)
               a.set_facecolor('none')
               plt.text(xnew_Xcparts[j],ynew_Xcparts[j],'swift',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.5*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax2.add_artist(a)
         else:
            if x_Xcparts_Swift>=(x-zoom) and x_Xcparts_Swift<=(x+zoom) and y_Xcparts_Swift>=(y-zoom) and y_Xcparts_Swift<=(y+zoom):
               xnew_Xcparts = x_Xcparts_Swift
               ynew_Xcparts = y_Xcparts_Swift
               errnew_Xcparts = err_Xcparts_Swift
               a = Circle(xy=(xnew_Xcparts, ynew_Xcparts), radius=errnew_Xcparts)
               a.set_linewidth(1.0)
               a.set_facecolor('none')
               plt.text(xnew_Xcparts,ynew_Xcparts,'swift',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.5*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax2.add_artist(a)
      if os.path.exists(path+'/csc2.cat'):
         if np.size(x_Xcparts_Chandra)!=1:
            xnew_Xcparts = [x_Xcparts_Chandra[j] for j in range(0,len(x_Xcparts_Chandra)) if x_Xcparts_Chandra[j]>=(x-zoom) and x_Xcparts_Chandra[j]<=(x+zoom) and y_Xcparts_Chandra[j]>=(y-zoom) and y_Xcparts_Chandra[j]<=(y+zoom)]
            ynew_Xcparts = [y_Xcparts_Chandra[j] for j in range(0,len(x_Xcparts_Chandra)) if x_Xcparts_Chandra[j]>=(x-zoom) and x_Xcparts_Chandra[j]<=(x+zoom) and y_Xcparts_Chandra[j]>=(y-zoom) and y_Xcparts_Chandra[j]<=(y+zoom)]
            ellanew_Xcparts = [ell_a_Chandra[j] for j in range(0,len(x_Xcparts_Chandra)) if x_Xcparts_Chandra[j]>=(x-zoom) and x_Xcparts_Chandra[j]<=(x+zoom) and y_Xcparts_Chandra[j]>=(y-zoom) and y_Xcparts_Chandra[j]<=(y+zoom)]
            ellbnew_Xcparts = [ell_b_Chandra[j] for j in range(0,len(x_Xcparts_Chandra)) if x_Xcparts_Chandra[j]>=(x-zoom) and x_Xcparts_Chandra[j]<=(x+zoom) and y_Xcparts_Chandra[j]>=(y-zoom) and y_Xcparts_Chandra[j]<=(y+zoom)]
            ellthetanew_Xcparts = [ell_theta_Chandra[j] for j in range(0,len(x_Xcparts_Chandra)) if x_Xcparts_Chandra[j]>=(x-zoom) and x_Xcparts_Chandra[j]<=(x+zoom) and y_Xcparts_Chandra[j]>=(y-zoom) and y_Xcparts_Chandra[j]<=(y+zoom)]
            for j in range(0,np.size(xnew_Xcparts)):
               a = Ellipse(xy=(xnew_Xcparts[j], ynew_Xcparts[j]), width=2*ellanew_Xcparts[j], height=2*ellbnew_Xcparts[j], angle=ellthetanew_Xcparts[j])
               a.set_linewidth(1.0)
               a.set_facecolor('none')
               plt.text(xnew_Xcparts[j],ynew_Xcparts[j],'chandra',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.5*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax2.add_artist(a)
         else:
            if x_Xcparts_Chandra>=(x-zoom) and x_Xcparts_Chandra<=(x+zoom) and y_Xcparts_Chandra>=(y-zoom) and y_Xcparts_Chandra<=(y+zoom):
               xnew_Xcparts = x_Xcparts_Chandra
               ynew_Xcparts = y_Xcparts_Chandra
               ellanew_Xcparts = ell_a_Chandra
               ellbnew_Xcparts = ell_b_Chandra
               ellthetanew_Xcparts = ell_theta_Chandra
               a = Ellipse(xy=(xnew_Xcparts, ynew_Xcparts), width=2*ellanew_Xcparts, height=2*ellbnew_Xcparts, angle=ellthetanew_Xcparts)
               a.set_linewidth(1.0)
               a.set_facecolor('none')
               plt.text(xnew_Xcparts,ynew_Xcparts,'chandra',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.5*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax2.add_artist(a)
      if os.path.exists(path+'/XMM4d13s.cat'):
         if np.size(x_Xcparts_XMM)!=1:
            xnew_Xcparts = [x_Xcparts_XMM[j] for j in range(0,len(x_Xcparts_XMM)) if x_Xcparts_XMM[j]>=(x-zoom) and x_Xcparts_XMM[j]<=(x+zoom) and y_Xcparts_XMM[j]>=(y-zoom) and y_Xcparts_XMM[j]<=(y+zoom)]
            ynew_Xcparts = [y_Xcparts_XMM[j] for j in range(0,len(x_Xcparts_XMM)) if x_Xcparts_XMM[j]>=(x-zoom) and x_Xcparts_XMM[j]<=(x+zoom) and y_Xcparts_XMM[j]>=(y-zoom) and y_Xcparts_XMM[j]<=(y+zoom)]
            errnew_Xcparts = [err_Xcparts_XMM[j] for j in range(0,len(x_Xcparts_XMM)) if x_Xcparts_XMM[j]>=(x-zoom) and x_Xcparts_XMM[j]<=(x+zoom) and y_Xcparts_XMM[j]>=(y-zoom) and y_Xcparts_XMM[j]<=(y+zoom)]
            for j in range(0,np.size(xnew_Xcparts)):
               a = Circle(xy=(xnew_Xcparts[j], ynew_Xcparts[j]), radius=errnew_Xcparts[j])
               a.set_linewidth(1.0)
               a.set_facecolor('none')
               plt.text(xnew_Xcparts[j],ynew_Xcparts[j],'XMM',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.5*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax2.add_artist(a)
         else:
            if x_Xcparts_XMM>=(x-zoom) and x_Xcparts_XMM<=(x+zoom) and y_Xcparts_XMM>=(y-zoom) and y_Xcparts_XMM<=(y+zoom):
               xnew_Xcparts = x_Xcparts_XMM
               ynew_Xcparts = y_Xcparts_XMM
               errnew_Xcparts = err_Xcparts_XMM
               a = Circle(xy=(xnew_Xcparts, ynew_Xcparts), radius=errnew_Xcparts)
               a.set_linewidth(1.0)
               a.set_facecolor('none')
               plt.text(xnew_Xcparts,ynew_Xcparts,'XMM',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.5*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax2.add_artist(a)
      if os.path.exists(path+'/eRASS1d1.cat'):
         if np.size(x_Xcparts_eRASS)!=1:
            xnew_Xcparts = [x_Xcparts_eRASS[j] for j in range(0,len(x_Xcparts_eRASS)) if x_Xcparts_eRASS[j]>=(x-zoom) and x_Xcparts_eRASS[j]<=(x+zoom) and y_Xcparts_eRASS[j]>=(y-zoom) and y_Xcparts_eRASS[j]<=(y+zoom)]
            ynew_Xcparts = [y_Xcparts_eRASS[j] for j in range(0,len(x_Xcparts_eRASS)) if x_Xcparts_eRASS[j]>=(x-zoom) and x_Xcparts_eRASS[j]<=(x+zoom) and y_Xcparts_eRASS[j]>=(y-zoom) and y_Xcparts_eRASS[j]<=(y+zoom)]
            errnew_Xcparts = [err_Xcparts_eRASS[j] for j in range(0,len(x_Xcparts_eRASS)) if x_Xcparts_eRASS[j]>=(x-zoom) and x_Xcparts_eRASS[j]<=(x+zoom) and y_Xcparts_eRASS[j]>=(y-zoom) and y_Xcparts_eRASS[j]<=(y+zoom)]
            for j in range(0,np.size(xnew_Xcparts)):
               a = Circle(xy=(xnew_Xcparts[j], ynew_Xcparts[j]), radius=errnew_Xcparts[j])
               a.set_linewidth(1.0)
               a.set_facecolor('none')
               plt.text(xnew_Xcparts[j],ynew_Xcparts[j],'eRASS',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.5*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax2.add_artist(a)
         else:
            if x_Xcparts_eRASS>=(x-zoom) and x_Xcparts_eRASS<=(x+zoom) and y_Xcparts_eRASS>=(y-zoom) and y_Xcparts_eRASS<=(y+zoom):
               xnew_Xcparts = x_Xcparts_eRASS
               ynew_Xcparts = y_Xcparts_eRASS
               errnew_Xcparts = err_Xcparts_eRASS
               a = Circle(xy=(xnew_Xcparts, ynew_Xcparts), radius=errnew_Xcparts)
               a.set_linewidth(1.0)
               a.set_facecolor('none')
               plt.text(xnew_Xcparts,ynew_Xcparts,'eRASS',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.5*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax2.add_artist(a)
      if np.size(x_Rcparts) > 0:
         if np.size(x_Rcparts) > 1:
            xnew_Rcparts = [x_Rcparts[j] for j in range(0,len(x_Rcparts)) if x_Rcparts[j]>=(x-zoom) and x_Rcparts[j]<=(x+zoom) and y_Rcparts[j]>=(y-zoom) and y_Rcparts[j]<=(y+zoom)]
            ynew_Rcparts = [y_Rcparts[j] for j in range(0,len(x_Rcparts)) if x_Rcparts[j]>=(x-zoom) and x_Rcparts[j]<=(x+zoom) and y_Rcparts[j]>=(y-zoom) and y_Rcparts[j]<=(y+zoom)]
            errmajnew_Rcparts = [errmaj_Rcparts[j] for j in range(0,len(x_Rcparts)) if x_Rcparts[j]>=(x-zoom) and x_Rcparts[j]<=(x+zoom) and y_Rcparts[j]>=(y-zoom) and y_Rcparts[j]<=(y+zoom)]
            errminnew_Rcparts = [errmin_Rcparts[j] for j in range(0,len(x_Rcparts)) if x_Rcparts[j]>=(x-zoom) and x_Rcparts[j]<=(x+zoom) and y_Rcparts[j]>=(y-zoom) and y_Rcparts[j]<=(y+zoom)]
            erranglenew_Rcparts = [errangle_Rcparts[j] for j in range(0,len(x_Rcparts)) if x_Rcparts[j]>=(x-zoom) and x_Rcparts[j]<=(x+zoom) and y_Rcparts[j]>=(y-zoom) and y_Rcparts[j]<=(y+zoom)]
            catnew_R = [cat_R[j] for j in range(0,len(x_Rcparts)) if x_Rcparts[j]>=(x-zoom) and x_Rcparts[j]<=(x+zoom) and y_Rcparts[j]>=(y-zoom) and y_Rcparts[j]<=(y+zoom)]
            for j in range(0,np.size(xnew_Rcparts)):
               if not np.isnan(errmajnew_Rcparts[j]) and not np.isnan(errminnew_Rcparts[j]) and not np.isnan(erranglenew_Rcparts[j]):
                  a = Ellipse(xy=(xnew_Rcparts[j], ynew_Rcparts[j]), width=2*errmajnew_Rcparts[j], height=2*errminnew_Rcparts[j], angle=erranglenew_Rcparts[j])
               elif not np.isnan(errmajnew_Rcparts[j]) and not np.isnan(errminnew_Rcparts[j]) and np.isnan(erranglenew_Rcparts[j]):
                  a = Circle(xy=(xnew_Rcparts[j], ynew_Rcparts[j]),radius=errmajnew_Rcparts[j])
               else:
                  a = Circle(xy=(xnew_Rcparts[j], ynew_Rcparts[j]),radius=2.)
               a.set_linewidth(1.0)
               a.set_facecolor('none')
               plt.text(xnew_Rcparts[j],ynew_Rcparts[j],catnew_R[j],horizontalalignment='center',verticalalignment='center',color='limegreen',fontsize=1.5*subsize,alpha=1.0)
               a.set_edgecolor('limegreen')
               ax2.add_artist(a)
         else:
            if x_Rcparts>=(x-zoom) and x_Rcparts<=(x+zoom) and y_Rcparts>=(y-zoom) and y_Rcparts<=(y+zoom):
               xnew_Rcparts = x_Rcparts
               ynew_Rcparts = y_Rcparts
               errmajnew_Rcparts = errmaj_Rcparts
               errminnew_Rcparts = errmin_Rcparts
               erranglenew_Rcparts = errangle_Rcparts
               catnew_R = cat_R
               if not np.isnan(errmajnew_Rcparts) and not np.isnan(errminnew_Rcparts) and not np.isnan(erranglenew_Rcparts):
                  a = Ellipse(xy=(xnew_Rcparts, ynew_Rcparts), width=2*errmajnew_Rcparts, height=2*errminnew_Rcparts, angle=erranglenew_Rcparts)
               elif not np.isnan(errmajnew_Rcparts) and not np.isnan(errminnew_Rcparts) and np.isnan(erranglenew_Rcparts):
                  a = Circle(xy=(xnew_Rcparts, ynew_Rcparts),radius=errmajnew_Rcparts)
               else:
                  a = Circle(xy=(xnew_Rcparts, ynew_Rcparts),radius=2.)
               a.set_linewidth(1.0)
               a.set_facecolor('none')
               plt.text(xnew_Rcparts,ynew_Rcparts,catnew_R,horizontalalignment='center',verticalalignment='center',color='limegreen',fontsize=1.5*subsize,alpha=1.0)
               a.set_edgecolor('limegreen')
               ax2.add_artist(a)

      if 'ra_query' in locals() and 'dec_query' in locals():
         plt.text(x_query,y_query,"query",horizontalalignment='center',verticalalignment='center',color='darkorange',fontsize=1.5*subsize,alpha=0.8)
         a = Circle(xy=(x_query, y_query), radius=8)
         a.set_linewidth(0.5)
         a.set_facecolor('none')
         a.set_edgecolor('darkorange')
         ax2.add_artist(a)

      if 'kic_ID' in locals():
         plt.text(-1.85,5.88,"Identified from Kepler as "+kic_ID+r" with $P_{\mathrm{Kepler}}$="+str(kic_period)+" d",transform=ax2.transAxes,fontsize=tsize)
      elif 'ATO_ID' in locals():
         plt.text(-1.85,5.88,"Identified from ATLAS as "+ATO_ID+r" with $P_{\mathrm{ATLAS}}$="+str(ATO_period)+" d",transform=ax2.transAxes,fontsize=tsize)
      elif 'SIMBAD_ID' in locals():
         plt.text(-1.7,5.88,"Identified from SIMBAD as "+SIMBAD_ID,transform=ax2.transAxes,fontsize=tsize)
      plt.savefig(path+'/finalperiodplots'+str(nlist)+'_'+per_search+'/finalperplot_'+str(i+1)+'_'+per_search+'_'+analysis+'_nper'+str(int(n_periods_sub))+'_double.png',bbox_inches='tight')
      plt.close(figure)

   if os.path.exists(path+'/rfilter/outputcats/doerPhot_V'+str(index)+'.csv') and os.path.exists(path+'/gfilter/outputcats/doerPhot_V'+str(index)+'.csv') and not os.path.exists(path+'/ifilter/outputcats/doerPhot_V'+str(index)+'.csv'):
      mjd_g, m_g, sigmam_g = np.genfromtxt(path+'/gfilter/outputcats/doerPhot_V'+str(index)+'.csv', delimiter=',', usecols=(8,-2,-1),unpack=True)
      mjd_r, m_r, sigmam_r = np.genfromtxt(path+'/rfilter/outputcats/doerPhot_V'+str(index)+'.csv', delimiter=',', usecols=(8,-2,-1),unpack=True)
      fwhm = np.genfromtxt(path+'/rfilter/outputcats/doerPhot_V'+str(index)+'.csv', delimiter=',', usecols=(6),unpack=True)

      #Take the light minimum in the r' band as t_0
      t_0 = mjd_r[np.argmax(m_r)]
      
      if per_search == 'LS':
         periods, periodogram = LombScargle_search(mjd_r, m_r, sigmam_r, P_low, P_up, n_periods_sub, 1)
         best_periods[i] = periods[~np.isnan(periodogram)][np.argmax(periodogram[~np.isnan(periodogram)])]
      elif per_search == 'PDM':
         periods, periodogram = PDM_search(mjd_r, m_r, P_low, P_up, 0., n_periods_sub, n_bin_sub, min_bin)
         best_periods[i] = periods[~np.isnan(periodogram)][np.argmin(periodogram[~np.isnan(periodogram)])]
      else:
         periods, periodogram = PDM_bisearch(mjd_g, m_g, mjd_r, m_r, P_low, P_up, 0., n_periods_sub, n_bin_sub, min_bin)
         best_periods[i] = periods[~np.isnan(periodogram)][np.argmin(periodogram[~np.isnan(periodogram)])]
      
      fwhmavg = np.mean(fwhm)
      
      nightdim = np.zeros(N_nights,int)
      for j in range(0,(np.size(mjd_r))):
         if j!=0:
            if (mjd_r[j]-mjd_r[j-1]) > 0.5:
               n_night = np.argmin(np.fabs(mjd_r[j-1]-mjd_ref))
               nightdim[n_night] = j-nightdim[0]
               for k in range(1,n_night):
                  nightdim[n_night] = nightdim[n_night]-nightdim[k]
            if j==np.size(mjd_r)-1:
               n_night = np.argmin(np.fabs(mjd_r[j]-mjd_ref))
               nightdim[n_night] = j+1-nightdim[0]
               for k in range(1,n_night):
                  nightdim[n_night] = nightdim[n_night]-nightdim[k]      

      mjd_single = []
      m_single = []
      sigmam_single = []
      n_lc = 0
      for j in range(0,N_nights):
         mjd_single.append(mjd_r[np.sum(nightdim[0:j]):np.sum(nightdim[0:j])+nightdim[j]])
         m_single.append(m_r[np.sum(nightdim[0:j]):np.sum(nightdim[0:j])+nightdim[j]])
         sigmam_single.append(sigmam_r[np.sum(nightdim[0:j]):np.sum(nightdim[0:j])+nightdim[j]])
         if nightdim[j]!=0:
            n_lc = n_lc+1

      figure=plt.figure(figsize=(figxsize,figysize), dpi=150)
      if np.size(sources_3FGL)!=0:
         if 'flux_Xcpart' in locals():
            figure.suptitle(r"Fermi source: "+FoVname+", Instrument: "+telescope+", Source detected in "+str(n_lc)+" nights out of "+str(N_nights)+"\n\nExposure time: %.1f"%deltat+" s, Plate scale: "+str(platescale)+" arcsec/pix, Extraction radius: "+str(fext)+"*FWHM\n\nPeriodic variable #%d"%(i+1)+", $\sigma$: %.2f"%sigma+" mag, $\Delta$m: %.2f"%dm+" mag, Range: [%.2f"%P_low+",%.1f"%min(P_up,mjd_r[-1]-mjd_r[0])+"] d, Steps: %d"%n_periods_sub+", Period folding with $T_{0}=$%.4f"%t_0+" MJD\n\n%.1f"%(sep_4FGL/ell_a_68)+r"$\sigma$ far from the 4FGL source and %.1f"%(sep_3FGL/ell_a_68_3FGL)+r"$\sigma$ far from the 3FGL source, X-ray flux: "+str(np.round(float(str(flux_Xcpart).split('e')[0]),1))+"$\cdot 10^{"+str(flux_Xcpart).split('e')[1]+"}$ erg s$^{-1}$ cm$^{-2}$", fontsize=tsize)
         else:
            if 'UL_X' in locals():
               figure.suptitle(r"Fermi source: "+FoVname+", Instrument: "+telescope+", Source detected in "+str(n_lc)+" nights out of "+str(N_nights)+"\n\nExposure time: %.1f"%deltat+" s, Plate scale: "+str(platescale)+" arcsec/pix, Extraction radius: "+str(fext)+"*FWHM\n\nPeriodic variable #%d"%(i+1)+", $\sigma$: %.2f"%sigma+" mag, $\Delta$m: %.2f"%dm+" mag, Range: [%.2f"%P_low+",%.1f"%min(P_up,mjd_r[-1]-mjd_r[0])+"] d, Steps: %d"%n_periods_sub+", Period folding with $T_{0}=$%.4f"%t_0+" MJD\n\n%.1f"%(sep_4FGL/ell_a_68)+r"$\sigma$ far from the 4FGL source and %.1f"%(sep_3FGL/ell_a_68_3FGL)+r"$\sigma$ far from the 3FGL source, X-ray flux 3$\sigma$ UL: "+str(np.round(float(str(UL_X).split('e')[0]),1))+r"$\cdot 10^{"+str(UL_X).split('e')[1]+"}$ erg s$^{-1}$ cm$^{-2}$", fontsize=tsize)
            else:
               figure.suptitle(r"Fermi source: "+FoVname+", Instrument: "+telescope+", Source detected in "+str(n_lc)+" nights out of "+str(N_nights)+"\n\nExposure time: %.1f"%deltat+" s, Plate scale: "+str(platescale)+" arcsec/pix, Extraction radius: "+str(fext)+"*FWHM\n\nPeriodic variable #%d"%(i+1)+", $\sigma$: %.2f"%sigma+" mag, $\Delta$m: %.2f"%dm+" mag, Range: [%.2f"%P_low+",%.1f"%min(P_up,mjd_r[-1]-mjd_r[0])+"] d, Steps: %d"%n_periods_sub+", Period folding with $T_{0}=$%.4f"%t_0+" MJD\n\n%.1f"%(sep_4FGL/ell_a_68)+r"$\sigma$ far from the 4FGL source and %.1f"%(sep_3FGL/ell_a_68_3FGL)+r"$\sigma$ far from the 3FGL source", fontsize=tsize)
      else:
         if 'flux_Xcpart' in locals():
            figure.suptitle(r"Fermi source: "+FoVname+", Instrument: "+telescope+", Source detected in "+str(n_lc)+" nights out of "+str(N_nights)+"\n\nExposure time: %.1f"%deltat+" s, Plate scale: "+str(platescale)+" arcsec/pix, Extraction radius: "+str(fext)+"*FWHM\n\nPeriodic variable #%d"%(i+1)+", $\sigma$: %.2f"%sigma+" mag, $\Delta$m: %.2f"%dm+" mag, Range: [%.2f"%P_low+",%.1f"%min(P_up,mjd_r[-1]-mjd_r[0])+"] d, Steps: %d"%n_periods_sub+", Period folding with $T_{0}=$%.4f"%t_0+" MJD\n\n%.1f"%(sep_4FGL/ell_a_68)+r"$\sigma$ far from the 4FGL source, X-ray flux: "+str(np.round(float(str(flux_Xcpart).split('e')[0]),1))+r"$\cdot 10^{"+str(flux_Xcpart).split('e')[1]+"}$ erg s$^{-1}$ cm$^{-2}$", fontsize=tsize)
         else:
            if 'UL_X' in locals():
               figure.suptitle(r"Fermi source: "+FoVname+", Instrument: "+telescope+", Source detected in "+str(n_lc)+" nights out of "+str(N_nights)+"\n\nExposure time: %.1f"%deltat+" s, Plate scale: "+str(platescale)+" arcsec/pix, Extraction radius: "+str(fext)+"*FWHM\n\nPeriodic variable #%d"%(i+1)+", $\sigma$: %.2f"%sigma+" mag, $\Delta$m: %.2f"%dm+" mag, Range: [%.2f"%P_low+",%.1f"%min(P_up,mjd_r[-1]-mjd_r[0])+"] d, Steps: %d"%n_periods_sub+", Period folding with $T_{0}=$%.4f"%t_0+" MJD\n\n%.1f"%(sep_4FGL/ell_a_68)+r"$\sigma$ far from the 4FGL source, X-ray flux 3$\sigma$ UL: "+str(np.round(float(str(UL_X).split('e')[0]),1))+r"$\cdot 10^{"+str(UL_X).split('e')[1]+"}$ erg s$^{-1}$ cm$^{-2}$", fontsize=tsize)
            else:
               figure.suptitle(r"Fermi source: "+FoVname+", Instrument: "+telescope+", Source detected in "+str(n_lc)+" nights out of "+str(N_nights)+"\n\nExposure time: %.1f"%deltat+" s, Plate scale: "+str(platescale)+" arcsec/pix, Extraction radius: "+str(fext)+"*FWHM\n\nPeriodic variable #%d"%(i+1)+", $\sigma$: %.2f"%sigma+" mag, $\Delta$m: %.2f"%dm+" mag, Range: [%.2f"%P_low+",%.1f"%min(P_up,mjd_r[-1]-mjd_r[0])+"] d, Steps: %d"%n_periods_sub+", Period folding with $T_{0}=$%.4f"%t_0+" MJD\n\n%.1f"%(sep_4FGL/ell_a_68)+r"$\sigma$ far from the 4FGL source", fontsize=tsize)
      
      counter = 0
      for j in range(0,N_nights):
         if nightdim[j]!=0:
            plt.subplot(int(np.trunc(4+n_lc/2)), 2, counter+1)
            plt.plot(mjd_single[j]-ref_ep, m_single[j], 'bo', linestyle='None')
            plt.errorbar(mjd_single[j]-ref_ep, m_single[j], sigmam_single[j], linestyle='None')
            plt.xticks(fontsize=axticksize)
            plt.yticks(fontsize=axticksize)
            plt.gca().invert_yaxis()
            plt.title("Lightcurve in the r band of "+nights[j],fontsize=subsize)
            plt.xlabel(r"Time (MJD-%d)"%ref_ep,fontsize=axlabelsize)
            plt.ylabel(r"Apparent magnitude",fontsize=axlabelsize)
            counter = counter+1

      plt.subplot(int(np.trunc(4+n_lc/2)), 2, n_lc+1)
      plt.xlabel(r'$\Delta$m (mag)',fontsize=axlabelsize)
      plt.ylabel(r'$\sigma$ (mag)',fontsize=axlabelsize)
      plt.title("Dispersion vs differential magnitude in the r band",fontsize=subsize)
      plt.yscale("log")
      plt.plot(sourcesr[:,2],sourcesr[:,3],'ob',markerfacecolor='0.8',markeredgecolor='black',markersize=markersize,alpha=alpha)
      plt.plot(varSources_fin1[:,2],varSources_fin1[:,3],'ob',markerfacecolor='red',markeredgecolor='black',markersize=markersize,alpha=alpha)
      #plt.plot(varSources_fin2[:,2],varSources_fin2[:,3],'ob',markerfacecolor='green',markeredgecolor='black',markersize=markersize,alpha=alpha)
      if nlist==1:
         if nvar!=1:
            plt.plot(perSources_fin[i,2],perSources_fin[i,3],'ob',markerfacecolor='red',markeredgecolor='black',markersize=markersize_main)
         else:
            plt.plot(perSources_fin[2],perSources_fin[3],'ob',markerfacecolor='red',markeredgecolor='black',markersize=markersize_main)
      else:
         if nvar!=1:
            plt.plot(perSources_fin[i,2],perSources_fin[i,3],'ob',markerfacecolor='green',markeredgecolor='black',markersize=markersize_main)
         else:
            plt.plot(perSources_fin[2],perSources_fin[3],'ob',markerfacecolor='green',markeredgecolor='black',markersize=markersize_main)
      plt.xticks(fontsize=axticksize)
      plt.yticks(fontsize=axticksize)

      mjd_gr = []
      g_r = []
      sigma_gr = []
      for j in range(0,np.size(mjd_g)):
         for k in range(0,np.size(mjd_r)):
            if np.size(mjd_g)!=1 and np.size(mjd_r)!=1:
               if (mjd_r[k]-mjd_g[j]) < t_space*deltat/86400. and (mjd_r[k]-mjd_g[j])>0:
                  mjd_gr.append((mjd_g[j]+mjd_r[k])/2.)
                  g_r.append(m_g[j]-m_r[k])
                  sigma_gr.append(np.sqrt(sigmam_g[j]**2+sigmam_r[k]**2))
      mjd_gr = np.asarray(mjd_gr)
      g_r = np.asarray(g_r)
      sigma_gr = np.asarray(sigma_gr)
      
      phase_g = ((mjd_g-t_0)/best_periods[i]) % 1
      #sort_index = np.asarray(phase_g.argsort())
      #phase_g = phase_g[sort_index]
      #m_g = m_g[sort_index]
      #sigmam_g = sigmam_g[sort_index]
      phase_r = ((mjd_r-t_0)/best_periods[i]) % 1
      #sort_index = np.asarray(phase_r.argsort())
      #phase_r = phase_r[sort_index]
      #m_r = m_r[sort_index]
      #sigmam_r = sigmam_r[sort_index]
      if np.size(mjd_g)!=1 and np.size(mjd_r)!=1:
         phase_gr = ((mjd_gr-t_0)/best_periods[i]) % 1
      #sort_index = np.asarray(phase_gr.argsort())
      #phase_gr = phase_gr[sort_index]
      #g_r = g_r[sort_index]
      #sigma_gr = sigma_gr[sort_index]

      plt.subplot(int(np.trunc(4+n_lc/2)), 2, n_lc+2)
      plt.plot(phase_r, m_r, 'bo', linestyle='None')
      plt.plot(phase_r+1, m_r, 'ro', linestyle='None')
      plt.errorbar(phase_r, m_r, sigmam_r, linestyle='None')
      plt.errorbar(phase_r+1, m_r, sigmam_r, linestyle='None')
      plt.xticks(fontsize=axticksize)
      plt.yticks(fontsize=axticksize)
      plt.gca().invert_yaxis()
      plt.title(r"Folded lightcurve in the r band: $P_{\mathrm{best}}$="+"{0} d".format(best_periods[i]),fontsize=subsize)
      plt.xlabel(r"Phase ($\phi$)",fontsize=axlabelsize)
      plt.ylabel(r"Apparent magnitude",fontsize=axlabelsize)

      plt.subplot(int(np.trunc(4+n_lc/2)), 2, n_lc+3)
      plt.plot(periods,periodogram)
      if not np.isnan(max(periodogram)) and not np.isinf(max(periodogram)) and not np.isnan(best_periods[i]) and not np.isinf(best_periods[i]):
         plt.axvline(x=best_periods[i], color='r', ls='--')
      plt.xscale("log")
      plt.xticks(fontsize=axticksize)
      plt.yticks(fontsize=axticksize)
      if per_search == 'LS':
         plt.title("Lomb-Scargle data periodogram in the r band",fontsize=subsize)
         plt.ylabel(r"Power",fontsize=axlabelsize)
         if period_FAPl == 'yes':
            tic_1 = time.perf_counter()
            FAP_periods, FAP_levels = LombScargle_FAP_level(mjd_r, m_r, sigmam_r, P_low, P_up, n_periods, 1, FAPl_factor, conf_level)
            toc_1 = time.perf_counter()
            print(f"FAP levels computed in {toc_1 - tic_1:0.4f} seconds")
            plt.plot(FAP_periods,FAP_levels,color='darkorange',ls='--')
            plt.text(plt.xlim()[0]*1.4, np.mean(FAP_levels)*1.2, r'0.1% FAP', color='darkorange', ha='center',fontsize=axticksize)
      elif per_search == 'PDM':
         plt.title("PDM data periodogram in the r band",fontsize=subsize)
         plt.ylabel(r"Normalized variance",fontsize=axlabelsize)
         if period_FAPl == 'yes':
            tic_1 = time.perf_counter()
            FAP_periods, FAP_levels = PDM_FAP_level(mjd_r, m_r, P_low, P_up, 0., n_periods_sub, n_bin_sub, min_bin, FAPl_factor, conf_level)
            toc_1 = time.perf_counter()
            print(f"FAP levels computed in {toc_1 - tic_1:0.4f} seconds")
            plt.plot(FAP_periods,FAP_levels,color='darkorange',ls='--')
            plt.text(plt.xlim()[0]*1.4, np.mean(FAP_levels)*1.02, r'0.1% FAP', color='darkorange', ha='center',fontsize=axticksize)
      else:
         plt.title("PDM multi-band data periodogram",fontsize=subsize)
         plt.ylabel(r"Sum of normalized variances",fontsize=axlabelsize)
      plt.xlabel(r"Trial period (days)",fontsize=axlabelsize)

      plt.subplot(int(np.trunc(4+n_lc/2)), 2, n_lc+4)
      plt.plot(phase_g, m_g, 'bo', linestyle='None')
      plt.plot(phase_g+1, m_g, 'ro', linestyle='None')
      plt.errorbar(phase_g, m_g, sigmam_g, linestyle='None')
      plt.errorbar(phase_g+1, m_g, sigmam_g, linestyle='None')
      plt.xticks(fontsize=axticksize)
      plt.yticks(fontsize=axticksize)
      plt.gca().invert_yaxis()
      plt.title(r"Folded lightcurve in the g band: $P_{\mathrm{best}}$="+"{0} d".format(best_periods[i]),fontsize=subsize)
      plt.xlabel(r"Phase ($\phi$)",fontsize=axlabelsize)
      plt.ylabel(r"Apparent magnitude",fontsize=axlabelsize)

      if np.size(mjd_g)!=1 and np.size(mjd_r)!=1:
         plt.subplot(int(np.trunc(4+n_lc/2)), 2, n_lc+5)
         plt.plot(phase_gr, g_r, 'bo', linestyle='None')
         plt.plot(phase_gr+1, g_r, 'ro', linestyle='None')
         plt.errorbar(phase_gr, g_r, sigma_gr, linestyle='None')
         plt.errorbar(phase_gr+1, g_r, sigma_gr, linestyle='None')
         plt.xticks(fontsize=axticksize)
         plt.yticks(fontsize=axticksize)
         plt.gca().invert_yaxis()
         plt.title("(g-r) folded colour curve",fontsize=subsize)
         plt.xlabel(r"Phase ($\phi$)",fontsize=axlabelsize)
         plt.ylabel(r"Colour",fontsize=axlabelsize)

      ax1 = plt.subplot(int(np.trunc(4+n_lc/2)), 2, n_lc+6)
      plt.title(r'Median r FoV bkg subtracted in the '+str(2*bound)+'$\sigma$ square around the 4FGL coordinates',fontsize=subsize)
      plt.xlim(max(0,xtar-bound*max(ell_a_95,ell_b_95)), min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)))
      plt.ylim(max(0,ytar-bound*max(ell_a_95,ell_b_95)), min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)))
      plt.imshow(data, interpolation='nearest', cmap='gray', vmin=m-colscale*s, vmax=m+colscale*s, origin='lower')
      plt.text(xtar,ytar,FoVname,horizontalalignment='center',verticalalignment='center',color='cyan',fontsize=subsize,alpha=0.8)
      target = Ellipse(xy=(xtar, ytar), width=2*ell_a_95, height=2*ell_b_95, angle=ell_theta_95)
      target.set_facecolor('none')
      target.set_edgecolor('cyan')
      ax1.add_artist(target)
      if np.size(sources_3FGL)!=0:
         for j in range(0,np.shape(sources_3FGL)[0]):
            if x_3FGL[j]>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_3FGL[j]<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_3FGL[j]>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_3FGL[j]<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
               plt.text(x_3FGL[j],y_3FGL[j],sources_3FGL[j].replace("L ","L"),horizontalalignment='center',verticalalignment='center',color='yellow',fontsize=subsize,alpha=0.8)
               a = Ellipse(xy=(x_3FGL[j], y_3FGL[j]), width=2*3600./float(platescale)*coord_3FGL[j,3], height=2*3600./float(platescale)*coord_3FGL[j,4], angle=coord_3FGL[j,5])
               a.set_facecolor('none')
               a.set_edgecolor('yellow')
               ax1.add_artist(a)
      #if np.size(sources_4FGLe)!=0:
      #   for j in range(0,np.shape(sources_4FGLe)[0]):
      #      if x_4FGLe[j]>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_4FGLe[j]<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_4FGLe[j]>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_4FGLe[j]<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
      #         plt.text(x_4FGLe[j],y_4FGLe[j],sources_4FGLe[j].replace("L ","L"),horizontalalignment='center',verticalalignment='center',color='cyan',fontsize=subsize,alpha=0.8)
      #         a = Ellipse(xy=(x_4FGLe[j], y_4FGLe[j]), width=2*3600./float(platescale)*np.asarray(coord_4FGLe)[j,2], height=2*3600./float(platescale)*np.asarray(coord_4FGLe)[j,3], angle=np.asarray(coord_4FGLe)[j,4])
      #         a.set_facecolor('none')
      #         a.set_edgecolor('cyan')
      #         ax1.add_artist(a)
         
      x, y = w.all_world2pix(ra, dec, 1)
      plt.text(x,y+50,str(i+1),horizontalalignment='center',verticalalignment='center',color='r',fontsize=subsize)
      a = Circle(xy=(x, y), radius=fext*fwhmavg)
      a.set_linewidth(0.5)
      a.set_facecolor('none')
      a.set_edgecolor('red')
      ax1.add_artist(a)
      #if hdr['INSTRUME']=='WIFSIP V2.0':
      #   if np.size(obs_Rcparts)!=0:
      #      for j in range(0,np.size(obs_Rcparts)):
      #         if x_Rcparts[j]>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Rcparts[j]<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Rcparts[j]>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Rcparts[j]<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
      #            a = Circle(xy=(x_Rcparts[j], y_Rcparts[j]), radius=8)
      #            a.set_linewidth(0.5)
      #            a.set_facecolor('none')
      #            plt.text(x_Rcparts[j],y_Rcparts[j]+50,obs_Rcparts[j],horizontalalignment='center',verticalalignment='center',color='limegreen',fontsize=1.2*subsize,alpha=0.8)
      #            a.set_edgecolor('limegreen')
      #            ax1.add_artist(a)
      
      if os.path.exists(path+'/2sxps.cat'):
         if np.size(x_Xcparts_Swift)!=1:
            for j in range(0,np.size(x_Xcparts_Swift)):
               if x_Xcparts_Swift[j]>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Xcparts_Swift[j]<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Xcparts_Swift[j]>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Xcparts_Swift[j]<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
                  a = Circle(xy=(x_Xcparts_Swift[j], y_Xcparts_Swift[j]), radius=err_Xcparts_Swift[j])
                  a.set_linewidth(0.5)
                  a.set_facecolor('none')
                  plt.text(x_Xcparts_Swift[j],y_Xcparts_Swift[j]+50,'swift',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.2*subsize,alpha=1.0)
                  a.set_edgecolor('fuchsia')
                  ax1.add_artist(a)
         else:
            if x_Xcparts_Swift>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Xcparts_Swift<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Xcparts_Swift>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Xcparts_Swift<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
               a = Circle(xy=(x_Xcparts_Swift, y_Xcparts_Swift), radius=err_Xcparts_Swift)
               a.set_linewidth(0.5)
               a.set_facecolor('none')
               plt.text(x_Xcparts_Swift,y_Xcparts_Swift+50,'swift',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.2*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax1.add_artist(a)
      if os.path.exists(path+'/csc2.cat'):
         if np.size(x_Xcparts_Chandra)!=1:
            for j in range(0,np.size(x_Xcparts_Chandra)):
               if x_Xcparts_Chandra[j]>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Xcparts_Chandra[j]<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Xcparts_Chandra[j]>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Xcparts_Chandra[j]<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
                  a = Ellipse(xy=(x_Xcparts_Chandra[j], y_Xcparts_Chandra[j]), width=2*ell_a_Chandra[j], height=2*ell_b_Chandra[j], angle=ell_theta_Chandra[j])
                  a.set_linewidth(0.5)
                  a.set_facecolor('none')
                  plt.text(x_Xcparts_Chandra[j],y_Xcparts_Chandra[j]+50,'chandra',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.2*subsize,alpha=1.0)
                  a.set_edgecolor('fuchsia')
                  ax1.add_artist(a)
         else:
            if x_Xcparts_Chandra>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Xcparts_Chandra<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Xcparts_Chandra>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Xcparts_Chandra<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
               a = Ellipse(xy=(x_Xcparts_Chandra, y_Xcparts_Chandra), width=2*ell_a_Chandra, height=2*ell_b_Chandra, angle=ell_theta_Chandra)
               a.set_linewidth(0.5)
               a.set_facecolor('none')
               plt.text(x_Xcparts_Chandra,y_Xcparts_Chandra+50,'chandra',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.2*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax1.add_artist(a)
      if os.path.exists(path+'/XMM4d13s.cat'):
         if np.size(x_Xcparts_XMM)!=1:
            for j in range(0,np.size(x_Xcparts_XMM)):
               if x_Xcparts_XMM[j]>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Xcparts_XMM[j]<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Xcparts_XMM[j]>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Xcparts_XMM[j]<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
                  a = Circle(xy=(x_Xcparts_XMM[j], y_Xcparts_XMM[j]), radius=err_Xcparts_XMM[j])
                  a.set_linewidth(0.5)
                  a.set_facecolor('none')
                  plt.text(x_Xcparts_XMM[j],y_Xcparts_XMM[j]+50,'XMM',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.2*subsize,alpha=1.0)
                  a.set_edgecolor('fuchsia')
                  ax1.add_artist(a)
         else:
            if x_Xcparts_XMM>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Xcparts_XMM<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Xcparts_XMM>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Xcparts_XMM<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
               a = Circle(xy=(x_Xcparts_XMM, y_Xcparts_XMM), radius=err_Xcparts_XMM)
               a.set_linewidth(0.5)
               a.set_facecolor('none')
               plt.text(x_Xcparts_XMM,y_Xcparts_XMM+50,'XMM',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.2*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax1.add_artist(a)
      if os.path.exists(path+'/eRASS1d1.cat'):
         if np.size(x_Xcparts_eRASS)!=1:
            for j in range(0,np.size(x_Xcparts_eRASS)):
               if x_Xcparts_eRASS[j]>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Xcparts_eRASS[j]<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Xcparts_eRASS[j]>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Xcparts_eRASS[j]<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
                  a = Circle(xy=(x_Xcparts_eRASS[j], y_Xcparts_eRASS[j]), radius=err_Xcparts_eRASS[j])
                  a.set_linewidth(0.5)
                  a.set_facecolor('none')
                  plt.text(x_Xcparts_eRASS[j],y_Xcparts_eRASS[j]+50,'eRASS',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.2*subsize,alpha=1.0)
                  a.set_edgecolor('fuchsia')
                  ax1.add_artist(a)
         else:
            if x_Xcparts_eRASS>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Xcparts_eRASS<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Xcparts_eRASS>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Xcparts_eRASS<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
               a = Circle(xy=(x_Xcparts_eRASS, y_Xcparts_eRASS), radius=err_Xcparts_eRASS)
               a.set_linewidth(0.5)
               a.set_facecolor('none')
               plt.text(x_Xcparts_eRASS,y_Xcparts_eRASS+50,'eRASS',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.2*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax1.add_artist(a)
      if np.size(x_Rcparts) > 0:
         if np.size(x_Rcparts) > 1:
            for j in range(0,np.size(x_Rcparts)):
               if x_Rcparts[j]>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Rcparts[j]<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Rcparts[j]>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Rcparts[j]<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
                  if not np.isnan(errmaj_Rcparts[j]) and not np.isnan(errmin_Rcparts[j]) and not np.isnan(errangle_Rcparts[j]):
                     a = Ellipse(xy=(x_Rcparts[j], y_Rcparts[j]), width=2*errmaj_Rcparts[j], height=2*errmin_Rcparts[j], angle=errangle_Rcparts[j])
                  elif not np.isnan(errmaj_Rcparts[j]) and not np.isnan(errmin_Rcparts[j]) and np.isnan(errangle_Rcparts[j]):
                     a = Circle(xy=(x_Rcparts[j], y_Rcparts[j]),radius=errmaj_Rcparts[j])
                  else:
                     a = Circle(xy=(x_Rcparts[j], y_Rcparts[j]),radius=2.)
                  a.set_linewidth(0.5)
                  a.set_facecolor('none')
                  plt.text(x_Rcparts[j],y_Rcparts[j]+50,cat_R[j],horizontalalignment='center',verticalalignment='center',color='limegreen',fontsize=1.2*subsize,alpha=1.0)
                  a.set_edgecolor('limegreen')
                  ax1.add_artist(a)
         else:
            if x_Rcparts>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Rcparts<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Rcparts>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Rcparts<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
               if not np.isnan(errmaj_Rcparts) and not np.isnan(errmin_Rcparts) and not np.isnan(errangle_Rcparts):
                  a = Ellipse(xy=(x_Rcparts, y_Rcparts), width=2*errmaj_Rcparts, height=2*errmin_Rcparts, angle=errangle_Rcparts)
               elif not np.isnan(errmaj_Rcparts) and not np.isnan(errmin_Rcparts) and np.isnan(errangle_Rcparts):
                  a = Circle(xy=(x_Rcparts, y_Rcparts),radius=errmaj_Rcparts)
               else:
                  a = Circle(xy=(x_Rcparts, y_Rcparts),radius=2.)
               a.set_linewidth(0.5)
               a.set_facecolor('none')
               plt.text(x_Rcparts,y_Rcparts+50,cat_R,horizontalalignment='center',verticalalignment='center',color='limegreen',fontsize=1.2*subsize,alpha=1.0)
               a.set_edgecolor('limegreen')
               ax1.add_artist(a)

      if 'ra_query' in locals() and 'dec_query' in locals():
         x_query, y_query = w.all_world2pix(ra_query, dec_query, 1)
         plt.text(x_query,y_query+50,"query",horizontalalignment='center',verticalalignment='center',color='darkorange',fontsize=1.2*subsize,alpha=1.0)
         a = Circle(xy=(x_query, y_query), radius=8)
         a.set_linewidth(0.5)
         a.set_facecolor('none')
         a.set_edgecolor('darkorange')
         ax1.add_artist(a)
         

      ax2 = plt.subplot(int(np.trunc(5+n_lc/2)), 2, n_lc+9)
      plt.title('Zoom on the variable source: (RA, DEC) = (%.3f'%ra+', %.3f'%dec+') deg',fontsize=subsize)
      plt.xlim(x-zoom,x+zoom)
      plt.ylim(y-zoom,y+zoom)
      plt.imshow(data, interpolation='nearest', cmap='gray', vmin=m-colscale*s, vmax=m+colscale*s, origin='lower')
      rxnew = [rx[j] for j in range(0,len(rx)) if rx[j]>=(x-zoom) and rx[j]<=(x+zoom) and ry[j]>=(y-zoom) and ry[j]<=(y+zoom)]
      rynew = [ry[j] for j in range(0,len(rx)) if rx[j]>=(x-zoom) and rx[j]<=(x+zoom) and ry[j]>=(y-zoom) and ry[j]<=(y+zoom)]
      for j in range(len(rxnew)):
         if not(np.round(x,ndig)==np.round(rxnew[j],ndig) and np.round(y,ndig)==np.round(rynew[j],ndig)):
            a = Circle(xy=(rxnew[j], rynew[j]), radius=fext*fwhmavg)
            a.set_linewidth(1.0)
            a.set_facecolor('none')
            a.set_edgecolor('blue')
            ax2.add_artist(a)
      plt.text(x,y+20,str(i+1),horizontalalignment='center',verticalalignment='center',color='r',fontsize=2*subsize)
      a = Circle(xy=(x, y), radius=fext*fwhmavg)
      a.set_linewidth(1.0)
      a.set_facecolor('none')
      a.set_edgecolor('red')
      ax2.add_artist(a)
      #if hdr['INSTRUME']=='WIFSIP V2.0':
      #   if np.size(obs_Rcparts)!=0:
      #      xnew_Rcparts = [x_Rcparts[j] for j in range(0,len(x_Rcparts)) if x_Rcparts[j]>=(x-zoom) and x_Rcparts[j]<=(x+zoom) and y_Rcparts[j]>=(y-zoom) and y_Rcparts[j]<=(y+zoom)]
      #      ynew_Rcparts = [y_Rcparts[j] for j in range(0,len(x_Rcparts)) if x_Rcparts[j]>=(x-zoom) and x_Rcparts[j]<=(x+zoom) and y_Rcparts[j]>=(y-zoom) and y_Rcparts[j]<=(y+zoom)]
      #      obsnew_Rcparts = [obs_Rcparts[j] for j in range(0,len(x_Rcparts)) if x_Rcparts[j]>=(x-zoom) and x_Rcparts[j]<=(x+zoom) and y_Rcparts[j]>=(y-zoom) and y_Rcparts[j]<=(y+zoom)]
      #      for j in range(0,np.size(xnew_Rcparts)):
      #         a = Circle(xy=(xnew_Rcparts[j], ynew_Rcparts[j]), radius=8)
      #         a.set_linewidth(1.0)
      #         a.set_facecolor('none')
      #         plt.text(xnew_Rcparts[j],ynew_Rcparts[j],obsnew_Rcparts[j],horizontalalignment='center',verticalalignment='center',color='limegreen',fontsize=1.5*subsize,alpha=1.0)
      #         a.set_edgecolor('limegreen')
      #         ax2.add_artist(a)
      if os.path.exists(path+'/2sxps.cat'):
         if np.size(x_Xcparts_Swift)!=1:
            xnew_Xcparts = [x_Xcparts_Swift[j] for j in range(0,len(x_Xcparts_Swift)) if x_Xcparts_Swift[j]>=(x-zoom) and x_Xcparts_Swift[j]<=(x+zoom) and y_Xcparts_Swift[j]>=(y-zoom) and y_Xcparts_Swift[j]<=(y+zoom)]
            ynew_Xcparts = [y_Xcparts_Swift[j] for j in range(0,len(x_Xcparts_Swift)) if x_Xcparts_Swift[j]>=(x-zoom) and x_Xcparts_Swift[j]<=(x+zoom) and y_Xcparts_Swift[j]>=(y-zoom) and y_Xcparts_Swift[j]<=(y+zoom)]
            errnew_Xcparts = [err_Xcparts_Swift[j] for j in range(0,len(x_Xcparts_Swift)) if x_Xcparts_Swift[j]>=(x-zoom) and x_Xcparts_Swift[j]<=(x+zoom) and y_Xcparts_Swift[j]>=(y-zoom) and y_Xcparts_Swift[j]<=(y+zoom)]
            for j in range(0,np.size(xnew_Xcparts)):
               a = Circle(xy=(xnew_Xcparts[j], ynew_Xcparts[j]), radius=errnew_Xcparts[j])
               a.set_linewidth(1.0)
               a.set_facecolor('none')
               plt.text(xnew_Xcparts[j],ynew_Xcparts[j],'swift',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.5*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax2.add_artist(a)
         else:
            if x_Xcparts_Swift>=(x-zoom) and x_Xcparts_Swift<=(x+zoom) and y_Xcparts_Swift>=(y-zoom) and y_Xcparts_Swift<=(y+zoom):
               xnew_Xcparts = x_Xcparts_Swift
               ynew_Xcparts = y_Xcparts_Swift
               errnew_Xcparts = err_Xcparts_Swift
               a = Circle(xy=(xnew_Xcparts, ynew_Xcparts), radius=errnew_Xcparts)
               a.set_linewidth(1.0)
               a.set_facecolor('none')
               plt.text(xnew_Xcparts,ynew_Xcparts,'swift',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.5*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax2.add_artist(a)
      if os.path.exists(path+'/csc2.cat'):
         if np.size(x_Xcparts_Chandra)!=1:
            xnew_Xcparts = [x_Xcparts_Chandra[j] for j in range(0,len(x_Xcparts_Chandra)) if x_Xcparts_Chandra[j]>=(x-zoom) and x_Xcparts_Chandra[j]<=(x+zoom) and y_Xcparts_Chandra[j]>=(y-zoom) and y_Xcparts_Chandra[j]<=(y+zoom)]
            ynew_Xcparts = [y_Xcparts_Chandra[j] for j in range(0,len(x_Xcparts_Chandra)) if x_Xcparts_Chandra[j]>=(x-zoom) and x_Xcparts_Chandra[j]<=(x+zoom) and y_Xcparts_Chandra[j]>=(y-zoom) and y_Xcparts_Chandra[j]<=(y+zoom)]
            ellanew_Xcparts = [ell_a_Chandra[j] for j in range(0,len(x_Xcparts_Chandra)) if x_Xcparts_Chandra[j]>=(x-zoom) and x_Xcparts_Chandra[j]<=(x+zoom) and y_Xcparts_Chandra[j]>=(y-zoom) and y_Xcparts_Chandra[j]<=(y+zoom)]
            ellbnew_Xcparts = [ell_b_Chandra[j] for j in range(0,len(x_Xcparts_Chandra)) if x_Xcparts_Chandra[j]>=(x-zoom) and x_Xcparts_Chandra[j]<=(x+zoom) and y_Xcparts_Chandra[j]>=(y-zoom) and y_Xcparts_Chandra[j]<=(y+zoom)]
            ellthetanew_Xcparts = [ell_theta_Chandra[j] for j in range(0,len(x_Xcparts_Chandra)) if x_Xcparts_Chandra[j]>=(x-zoom) and x_Xcparts_Chandra[j]<=(x+zoom) and y_Xcparts_Chandra[j]>=(y-zoom) and y_Xcparts_Chandra[j]<=(y+zoom)]
            for j in range(0,np.size(xnew_Xcparts)):
               a = Ellipse(xy=(xnew_Xcparts[j], ynew_Xcparts[j]), width=2*ellanew_Xcparts[j], height=2*ellbnew_Xcparts[j], angle=ellthetanew_Xcparts[j])
               a.set_linewidth(1.0)
               a.set_facecolor('none')
               plt.text(xnew_Xcparts[j],ynew_Xcparts[j],'chandra',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.5*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax2.add_artist(a)
         else:
            if x_Xcparts_Chandra>=(x-zoom) and x_Xcparts_Chandra<=(x+zoom) and y_Xcparts_Chandra>=(y-zoom) and y_Xcparts_Chandra<=(y+zoom):
               xnew_Xcparts = x_Xcparts_Chandra
               ynew_Xcparts = y_Xcparts_Chandra
               ellanew_Xcparts = ell_a_Chandra
               ellbnew_Xcparts = ell_b_Chandra
               ellthetanew_Xcparts = ell_theta_Chandra
               a = Ellipse(xy=(xnew_Xcparts, ynew_Xcparts), width=2*ellanew_Xcparts, height=2*ellbnew_Xcparts, angle=ellthetanew_Xcparts)
               a.set_linewidth(1.0)
               a.set_facecolor('none')
               plt.text(xnew_Xcparts,ynew_Xcparts,'chandra',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.5*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax2.add_artist(a)
      if os.path.exists(path+'/XMM4d13s.cat'):
         if np.size(x_Xcparts_XMM)!=1:
            xnew_Xcparts = [x_Xcparts_XMM[j] for j in range(0,len(x_Xcparts_XMM)) if x_Xcparts_XMM[j]>=(x-zoom) and x_Xcparts_XMM[j]<=(x+zoom) and y_Xcparts_XMM[j]>=(y-zoom) and y_Xcparts_XMM[j]<=(y+zoom)]
            ynew_Xcparts = [y_Xcparts_XMM[j] for j in range(0,len(x_Xcparts_XMM)) if x_Xcparts_XMM[j]>=(x-zoom) and x_Xcparts_XMM[j]<=(x+zoom) and y_Xcparts_XMM[j]>=(y-zoom) and y_Xcparts_XMM[j]<=(y+zoom)]
            errnew_Xcparts = [err_Xcparts_XMM[j] for j in range(0,len(x_Xcparts_XMM)) if x_Xcparts_XMM[j]>=(x-zoom) and x_Xcparts_XMM[j]<=(x+zoom) and y_Xcparts_XMM[j]>=(y-zoom) and y_Xcparts_XMM[j]<=(y+zoom)]
            for j in range(0,np.size(xnew_Xcparts)):
               a = Circle(xy=(xnew_Xcparts[j], ynew_Xcparts[j]), radius=errnew_Xcparts[j])
               a.set_linewidth(1.0)
               a.set_facecolor('none')
               plt.text(xnew_Xcparts[j],ynew_Xcparts[j],'XMM',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.5*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax2.add_artist(a)
         else:
            if x_Xcparts_XMM>=(x-zoom) and x_Xcparts_XMM<=(x+zoom) and y_Xcparts_XMM>=(y-zoom) and y_Xcparts_XMM<=(y+zoom):
               xnew_Xcparts = x_Xcparts_XMM
               ynew_Xcparts = y_Xcparts_XMM
               errnew_Xcparts = err_Xcparts_XMM
               a = Circle(xy=(xnew_Xcparts, ynew_Xcparts), radius=errnew_Xcparts)
               a.set_linewidth(1.0)
               a.set_facecolor('none')
               plt.text(xnew_Xcparts,ynew_Xcparts,'XMM',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.5*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax2.add_artist(a)
      if os.path.exists(path+'/eRASS1d1.cat'):
         if np.size(x_Xcparts_eRASS)!=1:
            xnew_Xcparts = [x_Xcparts_eRASS[j] for j in range(0,len(x_Xcparts_eRASS)) if x_Xcparts_eRASS[j]>=(x-zoom) and x_Xcparts_eRASS[j]<=(x+zoom) and y_Xcparts_eRASS[j]>=(y-zoom) and y_Xcparts_eRASS[j]<=(y+zoom)]
            ynew_Xcparts = [y_Xcparts_eRASS[j] for j in range(0,len(x_Xcparts_eRASS)) if x_Xcparts_eRASS[j]>=(x-zoom) and x_Xcparts_eRASS[j]<=(x+zoom) and y_Xcparts_eRASS[j]>=(y-zoom) and y_Xcparts_eRASS[j]<=(y+zoom)]
            errnew_Xcparts = [err_Xcparts_eRASS[j] for j in range(0,len(x_Xcparts_eRASS)) if x_Xcparts_eRASS[j]>=(x-zoom) and x_Xcparts_eRASS[j]<=(x+zoom) and y_Xcparts_eRASS[j]>=(y-zoom) and y_Xcparts_eRASS[j]<=(y+zoom)]
            for j in range(0,np.size(xnew_Xcparts)):
               a = Circle(xy=(xnew_Xcparts[j], ynew_Xcparts[j]), radius=errnew_Xcparts[j])
               a.set_linewidth(1.0)
               a.set_facecolor('none')
               plt.text(xnew_Xcparts[j],ynew_Xcparts[j],'eRASS',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.5*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax2.add_artist(a)
         else:
            if x_Xcparts_eRASS>=(x-zoom) and x_Xcparts_eRASS<=(x+zoom) and y_Xcparts_eRASS>=(y-zoom) and y_Xcparts_eRASS<=(y+zoom):
               xnew_Xcparts = x_Xcparts_eRASS
               ynew_Xcparts = y_Xcparts_eRASS
               errnew_Xcparts = err_Xcparts_eRASS
               a = Circle(xy=(xnew_Xcparts, ynew_Xcparts), radius=errnew_Xcparts)
               a.set_linewidth(1.0)
               a.set_facecolor('none')
               plt.text(xnew_Xcparts,ynew_Xcparts,'eRASS',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.5*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax2.add_artist(a)
      if np.size(x_Rcparts) > 0:
         if np.size(x_Rcparts) > 1:
            xnew_Rcparts = [x_Rcparts[j] for j in range(0,len(x_Rcparts)) if x_Rcparts[j]>=(x-zoom) and x_Rcparts[j]<=(x+zoom) and y_Rcparts[j]>=(y-zoom) and y_Rcparts[j]<=(y+zoom)]
            ynew_Rcparts = [y_Rcparts[j] for j in range(0,len(x_Rcparts)) if x_Rcparts[j]>=(x-zoom) and x_Rcparts[j]<=(x+zoom) and y_Rcparts[j]>=(y-zoom) and y_Rcparts[j]<=(y+zoom)]
            errmajnew_Rcparts = [errmaj_Rcparts[j] for j in range(0,len(x_Rcparts)) if x_Rcparts[j]>=(x-zoom) and x_Rcparts[j]<=(x+zoom) and y_Rcparts[j]>=(y-zoom) and y_Rcparts[j]<=(y+zoom)]
            errminnew_Rcparts = [errmin_Rcparts[j] for j in range(0,len(x_Rcparts)) if x_Rcparts[j]>=(x-zoom) and x_Rcparts[j]<=(x+zoom) and y_Rcparts[j]>=(y-zoom) and y_Rcparts[j]<=(y+zoom)]
            erranglenew_Rcparts = [errangle_Rcparts[j] for j in range(0,len(x_Rcparts)) if x_Rcparts[j]>=(x-zoom) and x_Rcparts[j]<=(x+zoom) and y_Rcparts[j]>=(y-zoom) and y_Rcparts[j]<=(y+zoom)]
            catnew_R = [cat_R[j] for j in range(0,len(x_Rcparts)) if x_Rcparts[j]>=(x-zoom) and x_Rcparts[j]<=(x+zoom) and y_Rcparts[j]>=(y-zoom) and y_Rcparts[j]<=(y+zoom)]
            for j in range(0,np.size(xnew_Rcparts)):
               if not np.isnan(errmajnew_Rcparts[j]) and not np.isnan(errminnew_Rcparts[j]) and not np.isnan(erranglenew_Rcparts[j]):
                  a = Ellipse(xy=(xnew_Rcparts[j], ynew_Rcparts[j]), width=2*errmajnew_Rcparts[j], height=2*errminnew_Rcparts[j], angle=erranglenew_Rcparts[j])
               elif not np.isnan(errmajnew_Rcparts[j]) and not np.isnan(errminnew_Rcparts[j]) and np.isnan(erranglenew_Rcparts[j]):
                  a = Circle(xy=(xnew_Rcparts[j], ynew_Rcparts[j]),radius=errmajnew_Rcparts[j])
               else:
                  a = Circle(xy=(xnew_Rcparts[j], ynew_Rcparts[j]),radius=2.)
               a.set_linewidth(1.0)
               a.set_facecolor('none')
               plt.text(xnew_Rcparts[j],ynew_Rcparts[j],catnew_R[j],horizontalalignment='center',verticalalignment='center',color='limegreen',fontsize=1.5*subsize,alpha=1.0)
               a.set_edgecolor('limegreen')
               ax2.add_artist(a)
         else:
            if x_Rcparts>=(x-zoom) and x_Rcparts<=(x+zoom) and y_Rcparts>=(y-zoom) and y_Rcparts<=(y+zoom):
               xnew_Rcparts = x_Rcparts
               ynew_Rcparts = y_Rcparts
               errmajnew_Rcparts = errmaj_Rcparts
               errminnew_Rcparts = errmin_Rcparts
               erranglenew_Rcparts = errangle_Rcparts
               catnew_R = cat_R
               if not np.isnan(errmajnew_Rcparts) and not np.isnan(errminnew_Rcparts) and not np.isnan(erranglenew_Rcparts):
                  a = Ellipse(xy=(xnew_Rcparts, ynew_Rcparts), width=2*errmajnew_Rcparts, height=2*errminnew_Rcparts, angle=erranglenew_Rcparts)
               elif not np.isnan(errmajnew_Rcparts) and not np.isnan(errminnew_Rcparts) and np.isnan(erranglenew_Rcparts):
                  a = Circle(xy=(xnew_Rcparts, ynew_Rcparts),radius=errmajnew_Rcparts)
               else:
                  a = Circle(xy=(xnew_Rcparts, ynew_Rcparts),radius=2.)
               a.set_linewidth(1.0)
               a.set_facecolor('none')
               plt.text(xnew_Rcparts,ynew_Rcparts,catnew_R,horizontalalignment='center',verticalalignment='center',color='limegreen',fontsize=1.5*subsize,alpha=1.0)
               a.set_edgecolor('limegreen')
               ax2.add_artist(a)
      
      if 'ra_query' in locals() and 'dec_query' in locals():
         plt.text(x_query,y_query,"query",horizontalalignment='center',verticalalignment='center',color='darkorange',fontsize=1.5*subsize,alpha=0.8)
         a = Circle(xy=(x_query, y_query), radius=8)
         a.set_linewidth(0.5)
         a.set_facecolor('none')
         a.set_edgecolor('darkorange')
         ax2.add_artist(a)

      if 'kic_ID' in locals():
         plt.text(-1.85,5.88,"Identified from Kepler as "+kic_ID+r" with $P_{\mathrm{Kepler}}$="+str(kic_period)+" d",transform=ax2.transAxes,fontsize=tsize)
      elif 'ATO_ID' in locals():
         plt.text(-1.85,5.88,"Identified from ATLAS as "+ATO_ID+r" with $P_{\mathrm{ATLAS}}$="+str(ATO_period)+" d",transform=ax2.transAxes,fontsize=tsize)
      elif 'SIMBAD_ID' in locals():
         plt.text(-1.7,5.88,"Identified from SIMBAD as "+SIMBAD_ID,transform=ax2.transAxes,fontsize=tsize)
      plt.savefig(path+'/finalperiodplots'+str(nlist)+'_'+per_search+'/finalperplot_'+str(i+1)+'_'+per_search+'_'+analysis+'_nper'+str(int(n_periods_sub))+'.png',bbox_inches='tight')
      plt.close(figure)

   elif os.path.exists(path+'/rfilter/outputcats/doerPhot_V'+str(index)+'.csv') and not os.path.exists(path+'/gfilter/outputcats/doerPhot_V'+str(index)+'.csv') and os.path.exists(path+'/ifilter/outputcats/doerPhot_V'+str(index)+'.csv'):
      mjd_r, m_r, sigmam_r = np.genfromtxt(path+'/rfilter/outputcats/doerPhot_V'+str(index)+'.csv', delimiter=',', usecols=(8,-2,-1),unpack=True)
      mjd_i, m_i, sigmam_i = np.genfromtxt(path+'/ifilter/outputcats/doerPhot_V'+str(index)+'.csv', delimiter=',', usecols=(8,-2,-1),unpack=True)
      fwhm = np.genfromtxt(path+'/rfilter/outputcats/doerPhot_V'+str(index)+'.csv', delimiter=',', usecols=(6),unpack=True)

      #Take the light minimum in the r' band as t_0
      t_0 = mjd_r[np.argmax(m_r)]
      
      if per_search == 'LS':
         periods, periodogram = LombScargle_search(mjd_r, m_r, sigmam_r, P_low, P_up, n_periods_sub, 1)
         best_periods[i] = periods[~np.isnan(periodogram)][np.argmax(periodogram[~np.isnan(periodogram)])]
      elif per_search == 'PDM':
         periods, periodogram = PDM_search(mjd_r, m_r, P_low, P_up, 0., n_periods_sub, n_bin_sub, min_bin)
         best_periods[i] = periods[~np.isnan(periodogram)][np.argmin(periodogram[~np.isnan(periodogram)])]
      else:
         periods, periodogram = PDM_bisearch(mjd_r, m_r, mjd_i, m_i, P_low, P_up, 0., n_periods_sub, n_bin_sub, min_bin)
         best_periods[i] = periods[~np.isnan(periodogram)][np.argmin(periodogram[~np.isnan(periodogram)])]
      
      fwhmavg = np.mean(fwhm)
      
      nightdim = np.zeros(N_nights,int)
      for j in range(0,(np.size(mjd_r))):
         if j!=0:
            if (mjd_r[j]-mjd_r[j-1]) > 0.5:
               n_night = np.argmin(np.fabs(mjd_r[j-1]-mjd_ref))
               nightdim[n_night] = j-nightdim[0]
               for k in range(1,n_night):
                  nightdim[n_night] = nightdim[n_night]-nightdim[k]
            if j==np.size(mjd_r)-1:
               n_night = np.argmin(np.fabs(mjd_r[j]-mjd_ref))
               nightdim[n_night] = j+1-nightdim[0]
               for k in range(1,n_night):
                  nightdim[n_night] = nightdim[n_night]-nightdim[k]      

      mjd_single = []
      m_single = []
      sigmam_single = []
      n_lc = 0
      for j in range(0,N_nights):
         mjd_single.append(mjd_r[np.sum(nightdim[0:j]):np.sum(nightdim[0:j])+nightdim[j]])
         m_single.append(m_r[np.sum(nightdim[0:j]):np.sum(nightdim[0:j])+nightdim[j]])
         sigmam_single.append(sigmam_r[np.sum(nightdim[0:j]):np.sum(nightdim[0:j])+nightdim[j]])
         if nightdim[j]!=0:
            n_lc = n_lc+1

      figure=plt.figure(figsize=(figxsize,figysize), dpi=150)
      if np.size(sources_3FGL)!=0:
         if 'flux_Xcpart' in locals():
            figure.suptitle(r"Fermi source: "+FoVname+", Instrument: "+telescope+", Source detected in "+str(n_lc)+" nights out of "+str(N_nights)+"\n\nExposure time: %.1f"%deltat+" s, Plate scale: "+str(platescale)+" arcsec/pix, Extraction radius: "+str(fext)+"*FWHM\n\nPeriodic variable #%d"%(i+1)+", $\sigma$: %.2f"%sigma+" mag, $\Delta$m: %.2f"%dm+" mag, Range: [%.2f"%P_low+",%.1f"%min(P_up,mjd_r[-1]-mjd_r[0])+"] d, Steps: %d"%n_periods_sub+", Period folding with $T_{0}=$%.4f"%t_0+" MJD\n\n%.1f"%(sep_4FGL/ell_a_68)+r"$\sigma$ far from the 4FGL source and %.1f"%(sep_3FGL/ell_a_68_3FGL)+r"$\sigma$ far from the 3FGL source, X-ray flux: "+str(np.round(float(str(flux_Xcpart).split('e')[0]),1))+r"$\cdot 10^{"+str(flux_Xcpart).split('e')[1]+"}$ erg s$^{-1}$ cm$^{-2}$", fontsize=tsize)
         else:
            if 'UL_X' in locals():
               figure.suptitle(r"Fermi source: "+FoVname+", Instrument: "+telescope+", Source detected in "+str(n_lc)+" nights out of "+str(N_nights)+"\n\nExposure time: %.1f"%deltat+" s, Plate scale: "+str(platescale)+" arcsec/pix, Extraction radius: "+str(fext)+"*FWHM\n\nPeriodic variable #%d"%(i+1)+", $\sigma$: %.2f"%sigma+" mag, $\Delta$m: %.2f"%dm+" mag, Range: [%.2f"%P_low+",%.1f"%min(P_up,mjd_r[-1]-mjd_r[0])+"] d, Steps: %d"%n_periods_sub+", Period folding with $T_{0}=$%.4f"%t_0+" MJD\n\n%.1f"%(sep_4FGL/ell_a_68)+r"$\sigma$ far from the 4FGL source and %.1f"%(sep_3FGL/ell_a_68_3FGL)+r"$\sigma$ far from the 3FGL source, X-ray flux 3$\sigma$ UL: "+str(np.round(float(str(UL_X).split('e')[0]),1))+r"$\cdot 10^{"+str(UL_X).split('e')[1]+"}$ erg s$^{-1}$ cm$^{-2}$", fontsize=tsize)
            else:
               figure.suptitle(r"Fermi source: "+FoVname+", Instrument: "+telescope+", Source detected in "+str(n_lc)+" nights out of "+str(N_nights)+"\n\nExposure time: %.1f"%deltat+" s, Plate scale: "+str(platescale)+" arcsec/pix, Extraction radius: "+str(fext)+"*FWHM\n\nPeriodic variable #%d"%(i+1)+", $\sigma$: %.2f"%sigma+" mag, $\Delta$m: %.2f"%dm+" mag, Range: [%.2f"%P_low+",%.1f"%min(P_up,mjd_r[-1]-mjd_r[0])+"] d, Steps: %d"%n_periods_sub+", Period folding with $T_{0}=$%.4f"%t_0+" MJD\n\n%.1f"%(sep_4FGL/ell_a_68)+r"$\sigma$ far from the 4FGL source and %.1f"%(sep_3FGL/ell_a_68_3FGL)+r"$\sigma$ far from the 3FGL source", fontsize=tsize)
      else:
         if 'flux_Xcpart' in locals():
            figure.suptitle(r"Fermi source: "+FoVname+", Instrument: "+telescope+", Source detected in "+str(n_lc)+" nights out of "+str(N_nights)+"\n\nExposure time: %.1f"%deltat+" s, Plate scale: "+str(platescale)+" arcsec/pix, Extraction radius: "+str(fext)+"*FWHM\\n\neriodic variable #%d"%(i+1)+", $\sigma$: %.2f"%sigma+" mag, $\Delta$m: %.2f"%dm+" mag, Range: [%.2f"%P_low+",%.1f"%min(P_up,mjd_r[-1]-mjd_r[0])+"] d, Steps: %d"%n_periods_sub+", Period folding with $T_{0}=$%.4f"%t_0+" MJD\n\n%.1f"%(sep_4FGL/ell_a_68)+r"$\sigma$ far from the 4FGL source, X-ray flux: "+str(np.round(float(str(flux_Xcpart).split('e')[0]),1))+"$\cdot 10^{"+str(flux_Xcpart).split('e')[1]+"}$ erg s$^{-1}$ cm$^{-2}$", fontsize=tsize)
         else:
            if 'UL_X' in locals():
               figure.suptitle(r"Fermi source: "+FoVname+", Instrument: "+telescope+", Source detected in "+str(n_lc)+" nights out of "+str(N_nights)+"\n\nExposure time: %.1f"%deltat+" s, Plate scale: "+str(platescale)+" arcsec/pix, Extraction radius: "+str(fext)+"*FWHM\n\nPeriodic variable #%d"%(i+1)+", $\sigma$: %.2f"%sigma+" mag, $\Delta$m: %.2f"%dm+" mag, Range: [%.2f"%P_low+",%.1f"%min(P_up,mjd_r[-1]-mjd_r[0])+"] d, Steps: %d"%n_periods_sub+", Period folding with $T_{0}=$%.4f"%t_0+" MJD\n\n%.1f"%(sep_4FGL/ell_a_68)+r"$\sigma$ far from the 4FGL source, X-ray flux 3$\sigma$ UL: "+str(np.round(float(str(UL_X).split('e')[0]),1))+"$\cdot 10^{"+str(UL_X).split('e')[1]+"}$ erg s$^{-1}$ cm$^{-2}$", fontsize=tsize)
            else:
               figure.suptitle(r"Fermi source: "+FoVname+", Instrument: "+telescope+", Source detected in "+str(n_lc)+" nights out of "+str(N_nights)+"\n\nExposure time: %.1f"%deltat+" s, Plate scale: "+str(platescale)+" arcsec/pix, Extraction radius: "+str(fext)+"*FWHM\n\nPeriodic variable #%d"%(i+1)+", $\sigma$: %.2f"%sigma+" mag, $\Delta$m: %.2f"%dm+" mag, Range: [%.2f"%P_low+",%.1f"%min(P_up,mjd_r[-1]-mjd_r[0])+"] d, Steps: %d"%n_periods_sub+", Period folding with $T_{0}=$%.4f"%t_0+" MJD\n\n%.1f"%(sep_4FGL/ell_a_68)+r"$\sigma$ far from the 4FGL source", fontsize=tsize)
      
      counter = 0
      for j in range(0,N_nights):
         if nightdim[j]!=0:
            plt.subplot(int(np.trunc(4+n_lc/2)), 2, counter+1)
            plt.plot(mjd_single[j]-ref_ep, m_single[j], 'bo', linestyle='None')
            plt.errorbar(mjd_single[j]-ref_ep, m_single[j], sigmam_single[j], linestyle='None')
            plt.xticks(fontsize=axticksize)
            plt.yticks(fontsize=axticksize)
            plt.gca().invert_yaxis()
            plt.title("Lightcurve in the r band of "+nights[j],fontsize=subsize)
            plt.xlabel(r"Time (MJD-%d)"%ref_ep,fontsize=axlabelsize)
            plt.ylabel(r"Apparent magnitude",fontsize=axlabelsize)
            counter = counter+1

      plt.subplot(int(np.trunc(4+n_lc/2)), 2, n_lc+1)
      plt.xlabel(r'$\Delta$m (mag)',fontsize=axlabelsize)
      plt.ylabel(r'$\sigma$ (mag)',fontsize=axlabelsize)
      plt.title("Dispersion vs differential magnitude in the r band",fontsize=subsize)
      plt.yscale("log")
      plt.plot(sourcesr[:,2],sourcesr[:,3],'ob',markerfacecolor='0.8',markeredgecolor='black',markersize=markersize,alpha=alpha)
      plt.plot(varSources_fin1[:,2],varSources_fin1[:,3],'ob',markerfacecolor='red',markeredgecolor='black',markersize=markersize,alpha=alpha)
      #plt.plot(varSources_fin2[:,2],varSources_fin2[:,3],'ob',markerfacecolor='green',markeredgecolor='black',markersize=markersize,alpha=alpha)
      if nlist==1:
         if nvar!=1:
            plt.plot(perSources_fin[i,2],perSources_fin[i,3],'ob',markerfacecolor='red',markeredgecolor='black',markersize=markersize_main)
         else:
            plt.plot(perSources_fin[2],perSources_fin[3],'ob',markerfacecolor='red',markeredgecolor='black',markersize=markersize_main)
      else:
         if nvar!=1:
            plt.plot(perSources_fin[i,2],perSources_fin[i,3],'ob',markerfacecolor='green',markeredgecolor='black',markersize=markersize_main)
         else:
            plt.plot(perSources_fin[2],perSources_fin[3],'ob',markerfacecolor='green',markeredgecolor='black',markersize=markersize_main)
      plt.xticks(fontsize=axticksize)
      plt.yticks(fontsize=axticksize)

      mjd_ri = []
      r_i = []
      sigma_ri = []
      for j in range(0,np.size(mjd_r)):
         for k in range(0,np.size(mjd_i)):
            if np.size(mjd_r)!=1 and np.size(mjd_i)!=1:
               if (mjd_i[k]-mjd_r[j]) < t_space*deltat/86400. and (mjd_i[k]-mjd_r[j])>0:
                  mjd_ri.append((mjd_r[j]+mjd_i[k])/2.)
                  r_i.append(m_r[j]-m_i[k])
                  sigma_ri.append(np.sqrt(sigmam_r[j]**2+sigmam_i[k]**2))
      mjd_ri = np.asarray(mjd_ri)
      r_i = np.asarray(r_i)
      sigma_ri = np.asarray(sigma_ri)
      
      phase_i = ((mjd_i-t_0)/best_periods[i]) % 1
      #sort_index = np.asarray(phase_i.argsort())
      #phase_i = phase_i[sort_index]
      #m_i = m_i[sort_index]
      #sigmam_i = sigmam_i[sort_index]
      phase_r = ((mjd_r-t_0)/best_periods[i]) % 1
      #sort_index = np.asarray(phase_r.argsort())
      #phase_r = phase_r[sort_index]
      #m_r = m_r[sort_index]
      #sigmam_r = sigmam_r[sort_index]
      if np.size(mjd_r)!=1 and np.size(mjd_i)!=1:
         phase_ri = ((mjd_ri-t_0)/best_periods[i]) % 1
      #sort_index = np.asarray(phase_ri.argsort())
      #phase_ri = phase_ri[sort_index]
      #r_i = r_i[sort_index]
      #sigma_ri = sigma_ri[sort_index]

      plt.subplot(int(np.trunc(4+n_lc/2)), 2, n_lc+2)
      plt.plot(phase_r, m_r, 'bo', linestyle='None')
      plt.plot(phase_r+1, m_r, 'ro', linestyle='None')
      plt.errorbar(phase_r, m_r, sigmam_r, linestyle='None')
      plt.errorbar(phase_r+1, m_r, sigmam_r, linestyle='None')
      plt.xticks(fontsize=axticksize)
      plt.yticks(fontsize=axticksize)
      plt.gca().invert_yaxis()
      plt.title(r"Folded lightcurve in the r band: $P_{\mathrm{best}}$="+"{0} d".format(best_periods[i]),fontsize=subsize)
      plt.xlabel(r"Phase ($\phi$)",fontsize=axlabelsize)
      plt.ylabel(r"Apparent magnitude",fontsize=axlabelsize)

      plt.subplot(int(np.trunc(4+n_lc/2)), 2, n_lc+3)
      plt.plot(periods,periodogram)
      if not np.isnan(max(periodogram)) and not np.isinf(max(periodogram)) and not np.isnan(best_periods[i]) and not np.isinf(best_periods[i]):
         plt.axvline(x=best_periods[i], color='r', ls='--')
      plt.xscale("log")
      plt.xticks(fontsize=axticksize)
      plt.yticks(fontsize=axticksize)
      if per_search == 'LS':
         plt.title("Lomb-Scargle data periodogram in the r band",fontsize=subsize)
         plt.ylabel(r"Power",fontsize=axlabelsize)
         if period_FAPl == 'yes':
            tic_1 = time.perf_counter()
            FAP_periods, FAP_levels = LombScargle_FAP_level(mjd_r, m_r, sigmam_r, P_low, P_up, n_periods, 1, FAPl_factor, conf_level)
            toc_1 = time.perf_counter()
            print(f"FAP levels computed in {toc_1 - tic_1:0.4f} seconds")
            plt.plot(FAP_periods,FAP_levels,color='darkorange',ls='--')
            plt.text(plt.xlim()[0]*1.4, np.mean(FAP_levels)*1.2, r'0.1% FAP', color='darkorange', ha='center',fontsize=axticksize)
      elif per_search == 'PDM':
         plt.title("PDM data periodogram in the r band",fontsize=subsize)
         plt.ylabel(r"Normalized variance",fontsize=axlabelsize)
         if period_FAPl == 'yes':
            tic_1 = time.perf_counter()
            FAP_periods, FAP_levels = PDM_FAP_level(mjd_r, m_r, P_low, P_up, 0., n_periods_sub, n_bin_sub, min_bin, FAPl_factor, conf_level)
            toc_1 = time.perf_counter()
            print(f"FAP levels computed in {toc_1 - tic_1:0.4f} seconds")
            plt.plot(FAP_periods,FAP_levels,color='darkorange',ls='--')
            plt.text(plt.xlim()[0]*1.4, np.mean(FAP_levels)*1.02, r'0.1% FAP', color='darkorange', ha='center',fontsize=axticksize)
      else:
         plt.title("PDM multi-band data periodogram",fontsize=subsize)
         plt.ylabel(r"Sum of normalized variances",fontsize=axlabelsize)
      plt.xlabel(r"Trial period (days)",fontsize=axlabelsize)

      plt.subplot(int(np.trunc(4+n_lc/2)), 2, n_lc+4)
      plt.plot(phase_i, m_i, 'bo', linestyle='None')
      plt.plot(phase_i+1, m_i, 'ro', linestyle='None')
      plt.errorbar(phase_i, m_i, sigmam_i, linestyle='None')
      plt.errorbar(phase_i+1, m_i, sigmam_i, linestyle='None')
      plt.xticks(fontsize=axticksize)
      plt.yticks(fontsize=axticksize)
      plt.gca().invert_yaxis()
      plt.title(r"Folded lightcurve in the i band: $P_{\mathrm{best}}$="+"{0} d".format(best_periods[i]),fontsize=subsize)
      plt.xlabel(r"Phase ($\phi$)",fontsize=axlabelsize)
      plt.ylabel(r"Apparent magnitude",fontsize=axlabelsize)

      if np.size(mjd_r)!=1 and np.size(mjd_i)!=1:
         plt.subplot(int(np.trunc(4+n_lc/2)), 2, n_lc+5)
         plt.plot(phase_ri, r_i, 'bo', linestyle='None')
         plt.plot(phase_ri+1, r_i, 'ro', linestyle='None')
         plt.errorbar(phase_ri, r_i, sigma_ri, linestyle='None')
         plt.errorbar(phase_ri+1, r_i, sigma_ri, linestyle='None')
         plt.xticks(fontsize=axticksize)
         plt.yticks(fontsize=axticksize)
         plt.gca().invert_yaxis()
         plt.title("(r-i) folded colour curve",fontsize=subsize)
         plt.xlabel(r"Phase ($\phi$)",fontsize=axlabelsize)
         plt.ylabel(r"Colour",fontsize=axlabelsize)

      ax1 = plt.subplot(int(np.trunc(4+n_lc/2)), 2, n_lc+6)
      plt.title(r'Median r FoV bkg subtracted in the '+str(2*bound)+'$\sigma$ square around the 4FGL coordinates',fontsize=subsize)
      plt.xlim(max(0,xtar-bound*max(ell_a_95,ell_b_95)), min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)))
      plt.ylim(max(0,ytar-bound*max(ell_a_95,ell_b_95)), min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)))
      plt.imshow(data, interpolation='nearest', cmap='gray', vmin=m-colscale*s, vmax=m+colscale*s, origin='lower')
      plt.text(xtar,ytar,FoVname,horizontalalignment='center',verticalalignment='center',color='cyan',fontsize=subsize,alpha=0.8)
      target = Ellipse(xy=(xtar, ytar), width=2*ell_a_95, height=2*ell_b_95, angle=ell_theta_95)
      target.set_facecolor('none')
      target.set_edgecolor('cyan')
      ax1.add_artist(target)
      if np.size(sources_3FGL)!=0:
         for j in range(0,np.shape(sources_3FGL)[0]):
            if x_3FGL[j]>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_3FGL[j]<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_3FGL[j]>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_3FGL[j]<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
               plt.text(x_3FGL[j],y_3FGL[j],sources_3FGL[j].replace("L ","L"),horizontalalignment='center',verticalalignment='center',color='yellow',fontsize=subsize,alpha=0.8)
               a = Ellipse(xy=(x_3FGL[j], y_3FGL[j]), width=2*3600./float(platescale)*coord_3FGL[j,3], height=2*3600./float(platescale)*coord_3FGL[j,4], angle=coord_3FGL[j,5])
               a.set_facecolor('none')
               a.set_edgecolor('yellow')
               ax1.add_artist(a)
      #if np.size(sources_4FGLe)!=0:
      #   for j in range(0,np.shape(sources_4FGLe)[0]):
      #      if x_4FGLe[j]>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_4FGLe[j]<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_4FGLe[j]>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_4FGLe[j]<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
      #         plt.text(x_4FGLe[j],y_4FGLe[j],sources_4FGLe[j].replace("L ","L"),horizontalalignment='center',verticalalignment='center',color='cyan',fontsize=subsize,alpha=0.8)
      #         a = Ellipse(xy=(x_4FGLe[j], y_4FGLe[j]), width=2*3600./float(platescale)*np.asarray(coord_4FGLe)[j,2], height=2*3600./float(platescale)*np.asarray(coord_4FGLe)[j,3], angle=np.asarray(coord_4FGLe)[j,4])
      #         a.set_facecolor('none')
      #         a.set_edgecolor('cyan')
      #         ax1.add_artist(a)
         
      x, y = w.all_world2pix(ra, dec, 1)
      plt.text(x,y+50,str(i+1),horizontalalignment='center',verticalalignment='center',color='r',fontsize=subsize)
      a = Circle(xy=(x, y), radius=fext*fwhmavg)
      a.set_linewidth(0.5)
      a.set_facecolor('none')
      a.set_edgecolor('red')
      ax1.add_artist(a)
      #if hdr['INSTRUME']=='WIFSIP V2.0':
      #   if np.size(obs_Rcparts)!=0:
      #      for j in range(0,np.size(obs_Rcparts)):
      #         if x_Rcparts[j]>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Rcparts[j]<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Rcparts[j]>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Rcparts[j]<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
      #            a = Circle(xy=(x_Rcparts[j], y_Rcparts[j]), radius=8)
      #            a.set_linewidth(0.5)
      #            a.set_facecolor('none')
      #            plt.text(x_Rcparts[j],y_Rcparts[j]+50,obs_Rcparts[j],horizontalalignment='center',verticalalignment='center',color='limegreen',fontsize=1.2*subsize,alpha=0.8)
      #            a.set_edgecolor('limegreen')
      #            ax1.add_artist(a)
      
      if os.path.exists(path+'/2sxps.cat'):
         if np.size(x_Xcparts_Swift)!=1:
            for j in range(0,np.size(x_Xcparts_Swift)):
               if x_Xcparts_Swift[j]>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Xcparts_Swift[j]<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Xcparts_Swift[j]>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Xcparts_Swift[j]<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
                  a = Circle(xy=(x_Xcparts_Swift[j], y_Xcparts_Swift[j]), radius=err_Xcparts_Swift[j])
                  a.set_linewidth(0.5)
                  a.set_facecolor('none')
                  plt.text(x_Xcparts_Swift[j],y_Xcparts_Swift[j]+50,'swift',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.2*subsize,alpha=1.0)
                  a.set_edgecolor('fuchsia')
                  ax1.add_artist(a)
         else:
            if x_Xcparts_Swift>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Xcparts_Swift<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Xcparts_Swift>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Xcparts_Swift<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
               a = Circle(xy=(x_Xcparts_Swift, y_Xcparts_Swift), radius=err_Xcparts_Swift)
               a.set_linewidth(0.5)
               a.set_facecolor('none')
               plt.text(x_Xcparts_Swift,y_Xcparts_Swift+50,'swift',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.2*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax1.add_artist(a)
      if os.path.exists(path+'/csc2.cat'):
         if np.size(x_Xcparts_Chandra)!=1:
            for j in range(0,np.size(x_Xcparts_Chandra)):
               if x_Xcparts_Chandra[j]>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Xcparts_Chandra[j]<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Xcparts_Chandra[j]>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Xcparts_Chandra[j]<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
                  a = Ellipse(xy=(x_Xcparts_Chandra[j], y_Xcparts_Chandra[j]), width=2*ell_a_Chandra[j], height=2*ell_b_Chandra[j], angle=ell_theta_Chandra[j])
                  a.set_linewidth(0.5)
                  a.set_facecolor('none')
                  plt.text(x_Xcparts_Chandra[j],y_Xcparts_Chandra[j]+50,'chandra',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.2*subsize,alpha=1.0)
                  a.set_edgecolor('fuchsia')
                  ax1.add_artist(a)
         else:
            if x_Xcparts_Chandra>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Xcparts_Chandra<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Xcparts_Chandra>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Xcparts_Chandra<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
               a = Ellipse(xy=(x_Xcparts_Chandra, y_Xcparts_Chandra), width=2*ell_a_Chandra, height=2*ell_b_Chandra, angle=ell_theta_Chandra)
               a.set_linewidth(0.5)
               a.set_facecolor('none')
               plt.text(x_Xcparts_Chandra,y_Xcparts_Chandra+50,'chandra',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.2*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax1.add_artist(a)
      if os.path.exists(path+'/XMM4d13s.cat'):
         if np.size(x_Xcparts_XMM)!=1:
            for j in range(0,np.size(x_Xcparts_XMM)):
               if x_Xcparts_XMM[j]>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Xcparts_XMM[j]<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Xcparts_XMM[j]>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Xcparts_XMM[j]<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
                  a = Circle(xy=(x_Xcparts_XMM[j], y_Xcparts_XMM[j]), radius=err_Xcparts_XMM[j])
                  a.set_linewidth(0.5)
                  a.set_facecolor('none')
                  plt.text(x_Xcparts_XMM[j],y_Xcparts_XMM[j]+50,'XMM',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.2*subsize,alpha=1.0)
                  a.set_edgecolor('fuchsia')
                  ax1.add_artist(a)
         else:
            if x_Xcparts_XMM>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Xcparts_XMM<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Xcparts_XMM>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Xcparts_XMM<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
               a = Circle(xy=(x_Xcparts_XMM, y_Xcparts_XMM), radius=err_Xcparts_XMM)
               a.set_linewidth(0.5)
               a.set_facecolor('none')
               plt.text(x_Xcparts_XMM,y_Xcparts_XMM+50,'XMM',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.2*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax1.add_artist(a)
      if os.path.exists(path+'/eRASS1d1.cat'):
         if np.size(x_Xcparts_eRASS)!=1:
            for j in range(0,np.size(x_Xcparts_eRASS)):
               if x_Xcparts_eRASS[j]>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Xcparts_eRASS[j]<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Xcparts_eRASS[j]>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Xcparts_eRASS[j]<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
                  a = Circle(xy=(x_Xcparts_eRASS[j], y_Xcparts_eRASS[j]), radius=err_Xcparts_eRASS[j])
                  a.set_linewidth(0.5)
                  a.set_facecolor('none')
                  plt.text(x_Xcparts_eRASS[j],y_Xcparts_eRASS[j]+50,'eRASS',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.2*subsize,alpha=1.0)
                  a.set_edgecolor('fuchsia')
                  ax1.add_artist(a)
         else:
            if x_Xcparts_eRASS>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Xcparts_eRASS<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Xcparts_eRASS>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Xcparts_eRASS<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
               a = Circle(xy=(x_Xcparts_eRASS, y_Xcparts_eRASS), radius=err_Xcparts_eRASS)
               a.set_linewidth(0.5)
               a.set_facecolor('none')
               plt.text(x_Xcparts_eRASS,y_Xcparts_eRASS+50,'eRASS',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.2*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax1.add_artist(a)
      if np.size(x_Rcparts) > 0:
         if np.size(x_Rcparts) > 1:
            for j in range(0,np.size(x_Rcparts)):
               if x_Rcparts[j]>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Rcparts[j]<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Rcparts[j]>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Rcparts[j]<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
                  if not np.isnan(errmaj_Rcparts[j]) and not np.isnan(errmin_Rcparts[j]) and not np.isnan(errangle_Rcparts[j]):
                     a = Ellipse(xy=(x_Rcparts[j], y_Rcparts[j]), width=2*errmaj_Rcparts[j], height=2*errmin_Rcparts[j], angle=errangle_Rcparts[j])
                  elif not np.isnan(errmaj_Rcparts[j]) and not np.isnan(errmin_Rcparts[j]) and np.isnan(errangle_Rcparts[j]):
                     a = Circle(xy=(x_Rcparts[j], y_Rcparts[j]),radius=errmaj_Rcparts[j])
                  else:
                     a = Circle(xy=(x_Rcparts[j], y_Rcparts[j]),radius=2.)
                  a.set_linewidth(0.5)
                  a.set_facecolor('none')
                  plt.text(x_Rcparts[j],y_Rcparts[j]+50,cat_R[j],horizontalalignment='center',verticalalignment='center',color='limegreen',fontsize=1.2*subsize,alpha=1.0)
                  a.set_edgecolor('limegreen')
                  ax1.add_artist(a)
         else:
            if x_Rcparts>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Rcparts<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Rcparts>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Rcparts<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
               if not np.isnan(errmaj_Rcparts) and not np.isnan(errmin_Rcparts) and not np.isnan(errangle_Rcparts):
                  a = Ellipse(xy=(x_Rcparts, y_Rcparts), width=2*errmaj_Rcparts, height=2*errmin_Rcparts, angle=errangle_Rcparts)
               elif not np.isnan(errmaj_Rcparts) and not np.isnan(errmin_Rcparts) and np.isnan(errangle_Rcparts):
                  a = Circle(xy=(x_Rcparts, y_Rcparts),radius=errmaj_Rcparts)
               else:
                  a = Circle(xy=(x_Rcparts, y_Rcparts),radius=2.)
               a.set_linewidth(0.5)
               a.set_facecolor('none')
               plt.text(x_Rcparts,y_Rcparts+50,cat_R,horizontalalignment='center',verticalalignment='center',color='limegreen',fontsize=1.2*subsize,alpha=1.0)
               a.set_edgecolor('limegreen')
               ax1.add_artist(a)

      if 'ra_query' in locals() and 'dec_query' in locals():
         x_query, y_query = w.all_world2pix(ra_query, dec_query, 1)
         plt.text(x_query,y_query+50,"query",horizontalalignment='center',verticalalignment='center',color='darkorange',fontsize=1.2*subsize,alpha=1.0)
         a = Circle(xy=(x_query, y_query), radius=8)
         a.set_linewidth(0.5)
         a.set_facecolor('none')
         a.set_edgecolor('darkorange')
         ax1.add_artist(a)
         

      ax2 = plt.subplot(int(np.trunc(4+n_lc/2)), 2, n_lc+6)
      plt.title('Zoom on the variable source: (RA, DEC) = (%.3f'%ra+', %.3f'%dec+') deg',fontsize=subsize)
      plt.xlim(x-zoom,x+zoom)
      plt.ylim(y-zoom,y+zoom)
      plt.imshow(data, interpolation='nearest', cmap='gray', vmin=m-colscale*s, vmax=m+colscale*s, origin='lower')
      rxnew = [rx[j] for j in range(0,len(rx)) if rx[j]>=(x-zoom) and rx[j]<=(x+zoom) and ry[j]>=(y-zoom) and ry[j]<=(y+zoom)]
      rynew = [ry[j] for j in range(0,len(rx)) if rx[j]>=(x-zoom) and rx[j]<=(x+zoom) and ry[j]>=(y-zoom) and ry[j]<=(y+zoom)]
      for j in range(len(rxnew)):
         if not(np.round(x,ndig)==np.round(rxnew[j],ndig) and np.round(y,ndig)==np.round(rynew[j],ndig)):
            a = Circle(xy=(rxnew[j], rynew[j]), radius=fext*fwhmavg)
            a.set_linewidth(1.0)
            a.set_facecolor('none')
            a.set_edgecolor('blue')
            ax2.add_artist(a)
      plt.text(x,y+20,str(i+1),horizontalalignment='center',verticalalignment='center',color='r',fontsize=2*subsize)
      a = Circle(xy=(x, y), radius=fext*fwhmavg)
      a.set_linewidth(1.0)
      a.set_facecolor('none')
      a.set_edgecolor('red')
      ax2.add_artist(a)
      #if hdr['INSTRUME']=='WIFSIP V2.0':
      #   if np.size(obs_Rcparts)!=0:
      #      xnew_Rcparts = [x_Rcparts[j] for j in range(0,len(x_Rcparts)) if x_Rcparts[j]>=(x-zoom) and x_Rcparts[j]<=(x+zoom) and y_Rcparts[j]>=(y-zoom) and y_Rcparts[j]<=(y+zoom)]
      #      ynew_Rcparts = [y_Rcparts[j] for j in range(0,len(x_Rcparts)) if x_Rcparts[j]>=(x-zoom) and x_Rcparts[j]<=(x+zoom) and y_Rcparts[j]>=(y-zoom) and y_Rcparts[j]<=(y+zoom)]
      #      obsnew_Rcparts = [obs_Rcparts[j] for j in range(0,len(x_Rcparts)) if x_Rcparts[j]>=(x-zoom) and x_Rcparts[j]<=(x+zoom) and y_Rcparts[j]>=(y-zoom) and y_Rcparts[j]<=(y+zoom)]
      #      for j in range(0,np.size(xnew_Rcparts)):
      #         a = Circle(xy=(xnew_Rcparts[j], ynew_Rcparts[j]), radius=8)
      #         a.set_linewidth(1.0)
      #         a.set_facecolor('none')
      #         plt.text(xnew_Rcparts[j],ynew_Rcparts[j],obsnew_Rcparts[j],horizontalalignment='center',verticalalignment='center',color='limegreen',fontsize=1.5*subsize,alpha=1.0)
      #         a.set_edgecolor('limegreen')
      #         ax2.add_artist(a)
      if os.path.exists(path+'/2sxps.cat'):
         if np.size(x_Xcparts_Swift)!=1:
            xnew_Xcparts = [x_Xcparts_Swift[j] for j in range(0,len(x_Xcparts_Swift)) if x_Xcparts_Swift[j]>=(x-zoom) and x_Xcparts_Swift[j]<=(x+zoom) and y_Xcparts_Swift[j]>=(y-zoom) and y_Xcparts_Swift[j]<=(y+zoom)]
            ynew_Xcparts = [y_Xcparts_Swift[j] for j in range(0,len(x_Xcparts_Swift)) if x_Xcparts_Swift[j]>=(x-zoom) and x_Xcparts_Swift[j]<=(x+zoom) and y_Xcparts_Swift[j]>=(y-zoom) and y_Xcparts_Swift[j]<=(y+zoom)]
            errnew_Xcparts = [err_Xcparts_Swift[j] for j in range(0,len(x_Xcparts_Swift)) if x_Xcparts_Swift[j]>=(x-zoom) and x_Xcparts_Swift[j]<=(x+zoom) and y_Xcparts_Swift[j]>=(y-zoom) and y_Xcparts_Swift[j]<=(y+zoom)]
            for j in range(0,np.size(xnew_Xcparts)):
               a = Circle(xy=(xnew_Xcparts[j], ynew_Xcparts[j]), radius=errnew_Xcparts[j])
               a.set_linewidth(1.0)
               a.set_facecolor('none')
               plt.text(xnew_Xcparts[j],ynew_Xcparts[j],'swift',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.5*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax2.add_artist(a)
         else:
            if x_Xcparts_Swift>=(x-zoom) and x_Xcparts_Swift<=(x+zoom) and y_Xcparts_Swift>=(y-zoom) and y_Xcparts_Swift<=(y+zoom):
               xnew_Xcparts = x_Xcparts_Swift
               ynew_Xcparts = y_Xcparts_Swift
               errnew_Xcparts = err_Xcparts_Swift
               a = Circle(xy=(xnew_Xcparts, ynew_Xcparts), radius=errnew_Xcparts)
               a.set_linewidth(1.0)
               a.set_facecolor('none')
               plt.text(xnew_Xcparts,ynew_Xcparts,'swift',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.5*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax2.add_artist(a)
      if os.path.exists(path+'/csc2.cat'):
         if np.size(x_Xcparts_Chandra)!=1:
            xnew_Xcparts = [x_Xcparts_Chandra[j] for j in range(0,len(x_Xcparts_Chandra)) if x_Xcparts_Chandra[j]>=(x-zoom) and x_Xcparts_Chandra[j]<=(x+zoom) and y_Xcparts_Chandra[j]>=(y-zoom) and y_Xcparts_Chandra[j]<=(y+zoom)]
            ynew_Xcparts = [y_Xcparts_Chandra[j] for j in range(0,len(x_Xcparts_Chandra)) if x_Xcparts_Chandra[j]>=(x-zoom) and x_Xcparts_Chandra[j]<=(x+zoom) and y_Xcparts_Chandra[j]>=(y-zoom) and y_Xcparts_Chandra[j]<=(y+zoom)]
            ellanew_Xcparts = [ell_a_Chandra[j] for j in range(0,len(x_Xcparts_Chandra)) if x_Xcparts_Chandra[j]>=(x-zoom) and x_Xcparts_Chandra[j]<=(x+zoom) and y_Xcparts_Chandra[j]>=(y-zoom) and y_Xcparts_Chandra[j]<=(y+zoom)]
            ellbnew_Xcparts = [ell_b_Chandra[j] for j in range(0,len(x_Xcparts_Chandra)) if x_Xcparts_Chandra[j]>=(x-zoom) and x_Xcparts_Chandra[j]<=(x+zoom) and y_Xcparts_Chandra[j]>=(y-zoom) and y_Xcparts_Chandra[j]<=(y+zoom)]
            ellthetanew_Xcparts = [ell_theta_Chandra[j] for j in range(0,len(x_Xcparts_Chandra)) if x_Xcparts_Chandra[j]>=(x-zoom) and x_Xcparts_Chandra[j]<=(x+zoom) and y_Xcparts_Chandra[j]>=(y-zoom) and y_Xcparts_Chandra[j]<=(y+zoom)]
            for j in range(0,np.size(xnew_Xcparts)):
               a = Ellipse(xy=(xnew_Xcparts[j], ynew_Xcparts[j]), width=2*ellanew_Xcparts[j], height=2*ellbnew_Xcparts[j], angle=ellthetanew_Xcparts[j])
               a.set_linewidth(1.0)
               a.set_facecolor('none')
               plt.text(xnew_Xcparts[j],ynew_Xcparts[j],'chandra',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.5*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax2.add_artist(a)
         else:
            if x_Xcparts_Chandra>=(x-zoom) and x_Xcparts_Chandra<=(x+zoom) and y_Xcparts_Chandra>=(y-zoom) and y_Xcparts_Chandra<=(y+zoom):
               xnew_Xcparts = x_Xcparts_Chandra
               ynew_Xcparts = y_Xcparts_Chandra
               ellanew_Xcparts = ell_a_Chandra
               ellbnew_Xcparts = ell_b_Chandra
               ellthetanew_Xcparts = ell_theta_Chandra
               a = Ellipse(xy=(xnew_Xcparts, ynew_Xcparts), width=2*ellanew_Xcparts, height=2*ellbnew_Xcparts, angle=ellthetanew_Xcparts)
               a.set_linewidth(1.0)
               a.set_facecolor('none')
               plt.text(xnew_Xcparts,ynew_Xcparts,'chandra',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.5*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax2.add_artist(a)
      if os.path.exists(path+'/XMM4d13s.cat'):
         if np.size(x_Xcparts_XMM)!=1:
            xnew_Xcparts = [x_Xcparts_XMM[j] for j in range(0,len(x_Xcparts_XMM)) if x_Xcparts_XMM[j]>=(x-zoom) and x_Xcparts_XMM[j]<=(x+zoom) and y_Xcparts_XMM[j]>=(y-zoom) and y_Xcparts_XMM[j]<=(y+zoom)]
            ynew_Xcparts = [y_Xcparts_XMM[j] for j in range(0,len(x_Xcparts_XMM)) if x_Xcparts_XMM[j]>=(x-zoom) and x_Xcparts_XMM[j]<=(x+zoom) and y_Xcparts_XMM[j]>=(y-zoom) and y_Xcparts_XMM[j]<=(y+zoom)]
            errnew_Xcparts = [err_Xcparts_XMM[j] for j in range(0,len(x_Xcparts_XMM)) if x_Xcparts_XMM[j]>=(x-zoom) and x_Xcparts_XMM[j]<=(x+zoom) and y_Xcparts_XMM[j]>=(y-zoom) and y_Xcparts_XMM[j]<=(y+zoom)]
            for j in range(0,np.size(xnew_Xcparts)):
               a = Circle(xy=(xnew_Xcparts[j], ynew_Xcparts[j]), radius=errnew_Xcparts[j])
               a.set_linewidth(1.0)
               a.set_facecolor('none')
               plt.text(xnew_Xcparts[j],ynew_Xcparts[j],'XMM',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.5*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax2.add_artist(a)
         else:
            if x_Xcparts_XMM>=(x-zoom) and x_Xcparts_XMM<=(x+zoom) and y_Xcparts_XMM>=(y-zoom) and y_Xcparts_XMM<=(y+zoom):
               xnew_Xcparts = x_Xcparts_XMM
               ynew_Xcparts = y_Xcparts_XMM
               errnew_Xcparts = err_Xcparts_XMM
               a = Circle(xy=(xnew_Xcparts, ynew_Xcparts), radius=errnew_Xcparts)
               a.set_linewidth(1.0)
               a.set_facecolor('none')
               plt.text(xnew_Xcparts,ynew_Xcparts,'XMM',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.5*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax2.add_artist(a)
      if os.path.exists(path+'/eRASS1d1.cat'):
         if np.size(x_Xcparts_eRASS)!=1:
            xnew_Xcparts = [x_Xcparts_eRASS[j] for j in range(0,len(x_Xcparts_eRASS)) if x_Xcparts_eRASS[j]>=(x-zoom) and x_Xcparts_eRASS[j]<=(x+zoom) and y_Xcparts_eRASS[j]>=(y-zoom) and y_Xcparts_eRASS[j]<=(y+zoom)]
            ynew_Xcparts = [y_Xcparts_eRASS[j] for j in range(0,len(x_Xcparts_eRASS)) if x_Xcparts_eRASS[j]>=(x-zoom) and x_Xcparts_eRASS[j]<=(x+zoom) and y_Xcparts_eRASS[j]>=(y-zoom) and y_Xcparts_eRASS[j]<=(y+zoom)]
            errnew_Xcparts = [err_Xcparts_eRASS[j] for j in range(0,len(x_Xcparts_eRASS)) if x_Xcparts_eRASS[j]>=(x-zoom) and x_Xcparts_eRASS[j]<=(x+zoom) and y_Xcparts_eRASS[j]>=(y-zoom) and y_Xcparts_eRASS[j]<=(y+zoom)]
            for j in range(0,np.size(xnew_Xcparts)):
               a = Circle(xy=(xnew_Xcparts[j], ynew_Xcparts[j]), radius=errnew_Xcparts[j])
               a.set_linewidth(1.0)
               a.set_facecolor('none')
               plt.text(xnew_Xcparts[j],ynew_Xcparts[j],'eRASS',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.5*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax2.add_artist(a)
         else:
            if x_Xcparts_eRASS>=(x-zoom) and x_Xcparts_eRASS<=(x+zoom) and y_Xcparts_eRASS>=(y-zoom) and y_Xcparts_eRASS<=(y+zoom):
               xnew_Xcparts = x_Xcparts_eRASS
               ynew_Xcparts = y_Xcparts_eRASS
               errnew_Xcparts = err_Xcparts_eRASS
               a = Circle(xy=(xnew_Xcparts, ynew_Xcparts), radius=errnew_Xcparts)
               a.set_linewidth(1.0)
               a.set_facecolor('none')
               plt.text(xnew_Xcparts,ynew_Xcparts,'eRASS',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.5*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax2.add_artist(a)
      if np.size(x_Rcparts) > 0:
         if np.size(x_Rcparts) > 1:
            xnew_Rcparts = [x_Rcparts[j] for j in range(0,len(x_Rcparts)) if x_Rcparts[j]>=(x-zoom) and x_Rcparts[j]<=(x+zoom) and y_Rcparts[j]>=(y-zoom) and y_Rcparts[j]<=(y+zoom)]
            ynew_Rcparts = [y_Rcparts[j] for j in range(0,len(x_Rcparts)) if x_Rcparts[j]>=(x-zoom) and x_Rcparts[j]<=(x+zoom) and y_Rcparts[j]>=(y-zoom) and y_Rcparts[j]<=(y+zoom)]
            errmajnew_Rcparts = [errmaj_Rcparts[j] for j in range(0,len(x_Rcparts)) if x_Rcparts[j]>=(x-zoom) and x_Rcparts[j]<=(x+zoom) and y_Rcparts[j]>=(y-zoom) and y_Rcparts[j]<=(y+zoom)]
            errminnew_Rcparts = [errmin_Rcparts[j] for j in range(0,len(x_Rcparts)) if x_Rcparts[j]>=(x-zoom) and x_Rcparts[j]<=(x+zoom) and y_Rcparts[j]>=(y-zoom) and y_Rcparts[j]<=(y+zoom)]
            erranglenew_Rcparts = [errangle_Rcparts[j] for j in range(0,len(x_Rcparts)) if x_Rcparts[j]>=(x-zoom) and x_Rcparts[j]<=(x+zoom) and y_Rcparts[j]>=(y-zoom) and y_Rcparts[j]<=(y+zoom)]
            catnew_R = [cat_R[j] for j in range(0,len(x_Rcparts)) if x_Rcparts[j]>=(x-zoom) and x_Rcparts[j]<=(x+zoom) and y_Rcparts[j]>=(y-zoom) and y_Rcparts[j]<=(y+zoom)]
            for j in range(0,np.size(xnew_Rcparts)):
               if not np.isnan(errmajnew_Rcparts[j]) and not np.isnan(errminnew_Rcparts[j]) and not np.isnan(erranglenew_Rcparts[j]):
                  a = Ellipse(xy=(xnew_Rcparts[j], ynew_Rcparts[j]), width=2*errmajnew_Rcparts[j], height=2*errminnew_Rcparts[j], angle=erranglenew_Rcparts[j])
               elif not np.isnan(errmajnew_Rcparts[j]) and not np.isnan(errminnew_Rcparts[j]) and np.isnan(erranglenew_Rcparts[j]):
                  a = Circle(xy=(xnew_Rcparts[j], ynew_Rcparts[j]),radius=errmajnew_Rcparts[j])
               else:
                  a = Circle(xy=(xnew_Rcparts[j], ynew_Rcparts[j]),radius=2.)
               a.set_linewidth(1.0)
               a.set_facecolor('none')
               plt.text(xnew_Rcparts[j],ynew_Rcparts[j],catnew_R[j],horizontalalignment='center',verticalalignment='center',color='limegreen',fontsize=1.5*subsize,alpha=1.0)
               a.set_edgecolor('limegreen')
               ax2.add_artist(a)
         else:
            if x_Rcparts>=(x-zoom) and x_Rcparts<=(x+zoom) and y_Rcparts>=(y-zoom) and y_Rcparts<=(y+zoom):
               xnew_Rcparts = x_Rcparts
               ynew_Rcparts = y_Rcparts
               errmajnew_Rcparts = errmaj_Rcparts
               errminnew_Rcparts = errmin_Rcparts
               erranglenew_Rcparts = errangle_Rcparts
               catnew_R = cat_R
               if not np.isnan(errmajnew_Rcparts) and not np.isnan(errminnew_Rcparts) and not np.isnan(erranglenew_Rcparts):
                  a = Ellipse(xy=(xnew_Rcparts, ynew_Rcparts), width=2*errmajnew_Rcparts, height=2*errminnew_Rcparts, angle=erranglenew_Rcparts)
               elif not np.isnan(errmajnew_Rcparts) and not np.isnan(errminnew_Rcparts) and np.isnan(erranglenew_Rcparts):
                  a = Circle(xy=(xnew_Rcparts, ynew_Rcparts),radius=errmajnew_Rcparts)
               else:
                  a = Circle(xy=(xnew_Rcparts, ynew_Rcparts),radius=2.)
               a.set_linewidth(1.0)
               a.set_facecolor('none')
               plt.text(xnew_Rcparts,ynew_Rcparts,catnew_R,horizontalalignment='center',verticalalignment='center',color='limegreen',fontsize=1.5*subsize,alpha=1.0)
               a.set_edgecolor('limegreen')
               ax2.add_artist(a)

      if 'ra_query' in locals() and 'dec_query' in locals():
         plt.text(x_query,y_query,"query",horizontalalignment='center',verticalalignment='center',color='darkorange',fontsize=1.5*subsize,alpha=0.8)
         a = Circle(xy=(x_query, y_query), radius=8)
         a.set_linewidth(0.5)
         a.set_facecolor('none')
         a.set_edgecolor('darkorange')
         ax2.add_artist(a)

      if 'kic_ID' in locals():
         plt.text(-1.85,5.88,"Identified from Kepler as "+kic_ID+r" with $P_{\mathrm{Kepler}}$="+str(kic_period)+" d",transform=ax2.transAxes,fontsize=tsize)
      elif 'ATO_ID' in locals():
         plt.text(-1.85,5.88,"Identified from ATLAS as "+ATO_ID+r" with $P_{\mathrm{ATLAS}}$="+str(ATO_period)+" d",transform=ax2.transAxes,fontsize=tsize)
      elif 'SIMBAD_ID' in locals():
         plt.text(-1.7,5.88,"Identified from SIMBAD as "+SIMBAD_ID,transform=ax2.transAxes,fontsize=tsize)
         
      plt.savefig(path+'/finalperiodplots'+str(nlist)+'_'+per_search+'/finalperplot_'+str(i+1)+'_'+per_search+'_'+analysis+'_nper'+str(int(n_periods_sub))+'.png',bbox_inches='tight')
      plt.close(figure) 

   if os.path.exists(path+'/rfilter/outputcats/doerPhot_V'+str(index)+'.csv') and not os.path.exists(path+'/gfilter/outputcats/doerPhot_V'+str(index)+'.csv') and not os.path.exists(path+'/ifilter/outputcats/doerPhot_V'+str(index)+'.csv'):
      mjd_r, m_r, sigmam_r = np.genfromtxt(path+'/rfilter/outputcats/doerPhot_V'+str(index)+'.csv', delimiter=',', usecols=(8,-2,-1),unpack=True)
      fwhm = np.genfromtxt(path+'/rfilter/outputcats/doerPhot_V'+str(index)+'.csv', delimiter=',', usecols=(6),unpack=True)

      #Take the light minimum in the r' band as t_0
      t_0 = mjd_r[np.argmax(m_r)]
      
      if per_search == 'LS':
         periods, periodogram = LombScargle_search(mjd_r, m_r, sigmam_r, P_low, P_up, n_periods_sub, 1)
         best_periods[i] = periods[~np.isnan(periodogram)][np.argmax(periodogram[~np.isnan(periodogram)])]
      else:
         periods, periodogram = PDM_search(mjd_r, m_r, P_low, P_up, 0., n_periods_sub, n_bin_sub, min_bin)
         best_periods[i] = periods[~np.isnan(periodogram)][np.argmin(periodogram[~np.isnan(periodogram)])]
      
      fwhmavg = np.mean(fwhm)
      
      nightdim = np.zeros(N_nights,int)
      for j in range(0,(np.size(mjd_r))):
         if j!=0:
            if (mjd_r[j]-mjd_r[j-1]) > 0.5:
               n_night = np.argmin(np.fabs(mjd_r[j-1]-mjd_ref))
               nightdim[n_night] = j-nightdim[0]
               for k in range(1,n_night):
                  nightdim[n_night] = nightdim[n_night]-nightdim[k]
            if j==np.size(mjd_r)-1:
               n_night = np.argmin(np.fabs(mjd_r[j]-mjd_ref))
               nightdim[n_night] = j+1-nightdim[0]
               for k in range(1,n_night):
                  nightdim[n_night] = nightdim[n_night]-nightdim[k]      

      mjd_single = []
      m_single = []
      sigmam_single = []
      n_lc = 0
      for j in range(0,N_nights):
         mjd_single.append(mjd_r[np.sum(nightdim[0:j]):np.sum(nightdim[0:j])+nightdim[j]])
         m_single.append(m_r[np.sum(nightdim[0:j]):np.sum(nightdim[0:j])+nightdim[j]])
         sigmam_single.append(sigmam_r[np.sum(nightdim[0:j]):np.sum(nightdim[0:j])+nightdim[j]])
         if nightdim[j]!=0:
            n_lc = n_lc+1

      figure=plt.figure(figsize=(figxsize,figysize), dpi=150)
      if np.size(sources_3FGL)!=0:
         if 'flux_Xcpart' in locals():
            figure.suptitle(r"Fermi source: "+FoVname+", Instrument: "+telescope+", Source detected in "+str(n_lc)+" nights out of "+str(N_nights)+"\n\nExposure time: %.1f"%deltat+" s, Plate scale: "+str(platescale)+" arcsec/pix, Extraction radius: "+str(fext)+"*FWHM\n\nPeriodic variable #%d"%(i+1)+", $\sigma$: %.2f"%sigma+" mag, $\Delta$m: %.2f"%dm+" mag, Range: [%.2f"%P_low+",%.1f"%min(P_up,mjd_r[-1]-mjd_r[0])+"] d, Steps: %d"%n_periods_sub+", Period folding with $T_{0}=$%.4f"%t_0+" MJD\n\n%.1f"%(sep_4FGL/ell_a_68)+r"$\sigma$ far from the 4FGL source and %.1f"%(sep_3FGL/ell_a_68_3FGL)+r"$\sigma$ far from the 3FGL source, X-ray flux: "+str(np.round(float(str(flux_Xcpart).split('e')[0]),1))+"$\cdot 10^{"+str(flux_Xcpart).split('e')[1]+"}$ erg s$^{-1}$ cm$^{-2}$", fontsize=tsize)
         else:
            if 'UL_X' in locals():
               figure.suptitle(r"Fermi source: "+FoVname+", Instrument: "+telescope+", Source detected in "+str(n_lc)+" nights out of "+str(N_nights)+"\n\nExposure time: %.1f"%deltat+" s, Plate scale: "+str(platescale)+" arcsec/pix, Extraction radius: "+str(fext)+"*FWHM\n\nPeriodic variable #%d"%(i+1)+", $\sigma$: %.2f"%sigma+" mag, $\Delta$m: %.2f"%dm+" mag, Range: [%.2f"%P_low+",%.1f"%min(P_up,mjd_r[-1]-mjd_r[0])+"] d, Steps: %d"%n_periods_sub+", Period folding with $T_{0}=$%.4f"%t_0+" MJD\n\n%.1f"%(sep_4FGL/ell_a_68)+r"$\sigma$ far from the 4FGL source and %.1f"%(sep_3FGL/ell_a_68_3FGL)+r"$\sigma$ far from the 3FGL source, X-ray flux 3$\sigma$ UL: "+str(np.round(float(str(UL_X).split('e')[0]),1))+"$\cdot 10^{"+str(UL_X).split('e')[1]+"}$ erg s$^{-1}$ cm$^{-2}$", fontsize=tsize)
            else:
               figure.suptitle(r"Fermi source: "+FoVname+", Instrument: "+telescope+", Source detected in "+str(n_lc)+" nights out of "+str(N_nights)+"\n\nExposure time: %.1f"%deltat+" s, Plate scale: "+str(platescale)+" arcsec/pix, Extraction radius: "+str(fext)+"*FWHM\n\nPeriodic variable #%d"%(i+1)+", $\sigma$: %.2f"%sigma+" mag, $\Delta$m: %.2f"%dm+" mag, Range: [%.2f"%P_low+",%.1f"%min(P_up,mjd_r[-1]-mjd_r[0])+"] d, Steps: %d"%n_periods_sub+", Period folding with $T_{0}=$%.4f"%t_0+" MJD\n\n%.1f"%(sep_4FGL/ell_a_68)+r"$\sigma$ far from the 4FGL source and %.1f"%(sep_3FGL/ell_a_68_3FGL)+r"$\sigma$ far from the 3FGL source", fontsize=tsize)
      else:
         if 'flux_Xcpart' in locals():
            figure.suptitle(r"Fermi source: "+FoVname+", Instrument: "+telescope+", Source detected in "+str(n_lc)+" nights out of "+str(N_nights)+"\n\nExposure time: %.1f"%deltat+" s, Plate scale: "+str(platescale)+" arcsec/pix, Extraction radius: "+str(fext)+"*FWHM\n\nPeriodic variable #%d"%(i+1)+", $\sigma$: %.2f"%sigma+" mag, $\Delta$m: %.2f"%dm+" mag, Range: [%.2f"%P_low+",%.1f"%min(P_up,mjd_r[-1]-mjd_r[0])+"] d, Steps: %d"%n_periods_sub+", Period folding with $T_{0}=$%.4f"%t_0+" MJD\n\n%.1f"%(sep_4FGL/ell_a_68)+r"$\sigma$ far from the 4FGL source, X-ray flux: "+str(np.round(float(str(flux_Xcpart).split('e')[0]),1))+"$\cdot 10^{"+str(flux_Xcpart).split('e')[1]+"}$ erg s$^{-1}$ cm$^{-2}$", fontsize=tsize)
         else:
            if 'UL_X' in locals():
               figure.suptitle(r"Fermi source: "+FoVname+", Instrument: "+telescope+", Source detected in "+str(n_lc)+" nights out of "+str(N_nights)+"\n\nExposure time: %.1f"%deltat+" s, Plate scale: "+str(platescale)+" arcsec/pix, Extraction radius: "+str(fext)+"*FWHM\n\nPeriodic variable #%d"%(i+1)+", $\sigma$: %.2f"%sigma+" mag, $\Delta$m: %.2f"%dm+" mag, Range: [%.2f"%P_low+",%.1f"%min(P_up,mjd_r[-1]-mjd_r[0])+"] d, Steps: %d"%n_periods_sub+", Period folding with $T_{0}=$%.4f"%t_0+" MJD\n\n%.1f"%(sep_4FGL/ell_a_68)+r"$\sigma$ far from the 4FGL source, X-ray flux 3$\sigma$ UL: "+str(np.round(float(str(UL_X).split('e')[0]),1))+"$\cdot 10^{"+str(UL_X).split('e')[1]+"}$ erg s$^{-1}$ cm$^{-2}$", fontsize=tsize)
            else:
               figure.suptitle(r"Fermi source: "+FoVname+", Instrument: "+telescope+", Source detected in "+str(n_lc)+" nights out of "+str(N_nights)+"\n\nExposure time: %.1f"%deltat+" s, Plate scale: "+str(platescale)+" arcsec/pix, Extraction radius: "+str(fext)+"*FWHM\n\nPeriodic variable #%d"%(i+1)+", $\sigma$: %.2f"%sigma+" mag, $\Delta$m: %.2f"%dm+" mag, Range: [%.2f"%P_low+",%.1f"%min(P_up,mjd_r[-1]-mjd_r[0])+"] d, Steps: %d"%n_periods_sub+", Period folding with $T_{0}=$%.4f"%t_0+" MJD\n\n%.1f"%(sep_4FGL/ell_a_68)+r"$\sigma$ far from the 4FGL source", fontsize=tsize)
      
      counter = 0
      for j in range(0,N_nights):
         if nightdim[j]!=0:
            plt.subplot(int(np.trunc(3+n_lc/2)), 2, counter+1)
            plt.plot(mjd_single[j]-ref_ep, m_single[j], 'bo', linestyle='None')
            plt.errorbar(mjd_single[j]-ref_ep, m_single[j], sigmam_single[j], linestyle='None')
            plt.xticks(fontsize=axticksize)
            plt.yticks(fontsize=axticksize)
            plt.gca().invert_yaxis()
            plt.title("Lightcurve in the r band of "+nights[j],fontsize=subsize)
            plt.xlabel(r"Time (MJD-%d)"%ref_ep,fontsize=axlabelsize)
            plt.ylabel(r"Apparent magnitude",fontsize=axlabelsize)
            counter = counter+1

      plt.subplot(int(np.trunc(3+n_lc/2)), 2, n_lc+1)
      plt.xlabel(r'$\Delta$m (mag)',fontsize=axlabelsize)
      plt.ylabel(r'$\sigma$ (mag)',fontsize=axlabelsize)
      plt.title("Dispersion vs differential magnitude in the r band",fontsize=subsize)
      plt.yscale("log")
      plt.plot(sourcesr[:,2],sourcesr[:,3],'ob',markerfacecolor='0.8',markeredgecolor='black',markersize=markersize,alpha=alpha)
      plt.plot(varSources_fin1[:,2],varSources_fin1[:,3],'ob',markerfacecolor='red',markeredgecolor='black',markersize=markersize,alpha=alpha)
      #plt.plot(varSources_fin2[:,2],varSources_fin2[:,3],'ob',markerfacecolor='green',markeredgecolor='black',markersize=markersize,alpha=alpha)
      if nlist==1:
         if nvar!=1:
            plt.plot(perSources_fin[i,2],perSources_fin[i,3],'ob',markerfacecolor='red',markeredgecolor='black',markersize=markersize_main)
         else:
            plt.plot(perSources_fin[2],perSources_fin[3],'ob',markerfacecolor='red',markeredgecolor='black',markersize=markersize_main)
      else:
         if nvar!=1:
            plt.plot(perSources_fin[i,2],perSources_fin[i,3],'ob',markerfacecolor='green',markeredgecolor='black',markersize=markersize_main)
         else:
            plt.plot(perSources_fin[2],perSources_fin[3],'ob',markerfacecolor='green',markeredgecolor='black',markersize=markersize_main)
      plt.xticks(fontsize=axticksize)
      plt.yticks(fontsize=axticksize)

      phase_r = ((mjd_r-t_0)/best_periods[i]) % 1
      #sort_index = np.asarray(phase_r.argsort())
      #phase_r = phase_r[sort_index]
      #m_r = m_r[sort_index]
      #sigmam_r = sigmam_r[sort_index]

      plt.subplot(int(np.trunc(3+n_lc/2)), 2, n_lc+2)
      plt.plot(phase_r, m_r, 'bo', linestyle='None')
      plt.plot(phase_r+1, m_r, 'ro', linestyle='None')
      plt.errorbar(phase_r, m_r, sigmam_r, linestyle='None')
      plt.errorbar(phase_r+1, m_r, sigmam_r, linestyle='None')
      plt.xticks(fontsize=axticksize)
      plt.yticks(fontsize=axticksize)
      plt.gca().invert_yaxis()
      plt.title(r"Folded lightcurve in the r band: $P_{\mathrm{best}}$="+"{0} d".format(best_periods[i]),fontsize=subsize)
      plt.xlabel(r"Phase ($\phi$)",fontsize=axlabelsize)
      plt.ylabel(r"Apparent magnitude",fontsize=axlabelsize)

      plt.subplot(int(np.trunc(3+n_lc/2)), 2, n_lc+3)
      plt.plot(periods,periodogram)
      if not np.isnan(max(periodogram)) and not np.isinf(max(periodogram)) and not np.isnan(best_periods[i]) and not np.isinf(best_periods[i]):
         plt.axvline(x=best_periods[i], color='r', ls='--')
      plt.xscale("log")
      plt.xticks(fontsize=axticksize)
      plt.yticks(fontsize=axticksize)
      if per_search == 'LS':
         plt.title("Lomb-Scargle data periodogram in the r band",fontsize=subsize)
         plt.ylabel(r"Power",fontsize=axlabelsize)
         if period_FAPl == 'yes':
            tic_1 = time.perf_counter()
            FAP_periods, FAP_levels = LombScargle_FAP_level(mjd_r, m_r, sigmam_r, P_low, P_up, n_periods, 1, FAPl_factor, conf_level)
            toc_1 = time.perf_counter()
            print(f"FAP levels computed in {toc_1 - tic_1:0.4f} seconds")
            plt.plot(FAP_periods,FAP_levels,color='darkorange',ls='--')
            plt.text(plt.xlim()[0]*1.4, np.mean(FAP_levels)*1.2, r'0.1% FAP', color='darkorange', ha='center',fontsize=axticksize)
      else:
         plt.title("PDM data periodogram in the r band",fontsize=subsize)
         plt.ylabel(r"Normalized variance",fontsize=axlabelsize)
         if period_FAPl == 'yes':
            tic_1 = time.perf_counter()
            FAP_periods, FAP_levels = PDM_FAP_level(mjd_r, m_r, P_low, P_up, 0., n_periods_sub, n_bin_sub, min_bin, FAPl_factor, conf_level)
            toc_1 = time.perf_counter()
            print(f"FAP levels computed in {toc_1 - tic_1:0.4f} seconds")
            plt.plot(FAP_periods,FAP_levels,color='darkorange',ls='--')
            plt.text(plt.xlim()[0]*1.4, np.mean(FAP_levels)*1.02, r'0.1% FAP', color='darkorange', ha='center',fontsize=axticksize)
      plt.xlabel(r"Trial period (days)",fontsize=axlabelsize)

      ax1 = plt.subplot(int(np.trunc(3+n_lc/2)), 2, n_lc+4)
      plt.title(r'Median r FoV bkg subtracted in the '+str(2*bound)+'$\sigma$ square around the 4FGL coordinates',fontsize=subsize)
      plt.xlim(max(0,xtar-bound*max(ell_a_95,ell_b_95)), min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)))
      plt.ylim(max(0,ytar-bound*max(ell_a_95,ell_b_95)), min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)))
      plt.imshow(data, interpolation='nearest', cmap='gray', vmin=m-colscale*s, vmax=m+colscale*s, origin='lower')
      plt.text(xtar,ytar,FoVname,horizontalalignment='center',verticalalignment='center',color='cyan',fontsize=subsize,alpha=0.8)
      target = Ellipse(xy=(xtar, ytar), width=2*ell_a_95, height=2*ell_b_95, angle=ell_theta_95)
      target.set_facecolor('none')
      target.set_edgecolor('cyan')
      ax1.add_artist(target)
      if np.size(sources_3FGL)!=0:
         for j in range(0,np.shape(sources_3FGL)[0]):
            if x_3FGL[j]>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_3FGL[j]<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_3FGL[j]>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_3FGL[j]<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
               plt.text(x_3FGL[j],y_3FGL[j],sources_3FGL[j].replace("L ","L"),horizontalalignment='center',verticalalignment='center',color='yellow',fontsize=subsize,alpha=0.8)
               a = Ellipse(xy=(x_3FGL[j], y_3FGL[j]), width=2*3600./float(platescale)*coord_3FGL[j,3], height=2*3600./float(platescale)*coord_3FGL[j,4], angle=coord_3FGL[j,5])
               a.set_facecolor('none')
               a.set_edgecolor('yellow')
               ax1.add_artist(a)
      #if np.size(sources_4FGLe)!=0:
      #   for j in range(0,np.shape(sources_4FGLe)[0]):
      #      if x_4FGLe[j]>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_4FGLe[j]<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_4FGLe[j]>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_4FGLe[j]<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
      #         plt.text(x_4FGLe[j],y_4FGLe[j],sources_4FGLe[j].replace("L ","L"),horizontalalignment='center',verticalalignment='center',color='cyan',fontsize=subsize,alpha=0.8)
      #         a = Ellipse(xy=(x_4FGLe[j], y_4FGLe[j]), width=2*3600./float(platescale)*np.asarray(coord_4FGLe)[j,2], height=2*3600./float(platescale)*np.asarray(coord_4FGLe)[j,3], angle=np.asarray(coord_4FGLe)[j,4])
      #         a.set_facecolor('none')
      #         a.set_edgecolor('cyan')
      #         ax1.add_artist(a)
         
      x, y = w.all_world2pix(ra, dec, 1)
      plt.text(x,y+50,str(i+1),horizontalalignment='center',verticalalignment='center',color='r',fontsize=subsize)
      a = Circle(xy=(x, y), radius=fext*fwhmavg)
      a.set_linewidth(0.5)
      a.set_facecolor('none')
      a.set_edgecolor('red')
      ax1.add_artist(a)
      #if hdr['INSTRUME']=='WIFSIP V2.0':
      #   if np.size(obs_Rcparts)!=0:
      #      for j in range(0,np.size(obs_Rcparts)):
      #         if x_Rcparts[j]>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Rcparts[j]<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Rcparts[j]>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Rcparts[j]<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
      #            a = Circle(xy=(x_Rcparts[j], y_Rcparts[j]), radius=8)
      #            a.set_linewidth(0.5)
      #            a.set_facecolor('none')
      #            plt.text(x_Rcparts[j],y_Rcparts[j]+50,obs_Rcparts[j],horizontalalignment='center',verticalalignment='center',color='limegreen',fontsize=1.2*subsize,alpha=0.8)
      #            a.set_edgecolor('limegreen')
      #            ax1.add_artist(a)
      
      if os.path.exists(path+'/2sxps.cat'):
         if np.size(x_Xcparts_Swift)!=1:
            for j in range(0,np.size(x_Xcparts_Swift)):
               if x_Xcparts_Swift[j]>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Xcparts_Swift[j]<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Xcparts_Swift[j]>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Xcparts_Swift[j]<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
                  a = Circle(xy=(x_Xcparts_Swift[j], y_Xcparts_Swift[j]), radius=err_Xcparts_Swift[j])
                  a.set_linewidth(0.5)
                  a.set_facecolor('none')
                  plt.text(x_Xcparts_Swift[j],y_Xcparts_Swift[j]+50,'swift',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.2*subsize,alpha=1.0)
                  a.set_edgecolor('fuchsia')
                  ax1.add_artist(a)
         else:
            if x_Xcparts_Swift>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Xcparts_Swift<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Xcparts_Swift>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Xcparts_Swift<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
               a = Circle(xy=(x_Xcparts_Swift, y_Xcparts_Swift), radius=err_Xcparts_Swift)
               a.set_linewidth(0.5)
               a.set_facecolor('none')
               plt.text(x_Xcparts_Swift,y_Xcparts_Swift+50,'swift',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.2*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax1.add_artist(a)
      if os.path.exists(path+'/csc2.cat'):
         if np.size(x_Xcparts_Chandra)!=1:
            for j in range(0,np.size(x_Xcparts_Chandra)):
               if x_Xcparts_Chandra[j]>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Xcparts_Chandra[j]<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Xcparts_Chandra[j]>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Xcparts_Chandra[j]<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
                  a = Ellipse(xy=(x_Xcparts_Chandra[j], y_Xcparts_Chandra[j]), width=2*ell_a_Chandra[j], height=2*ell_b_Chandra[j], angle=ell_theta_Chandra[j])
                  a.set_linewidth(0.5)
                  a.set_facecolor('none')
                  plt.text(x_Xcparts_Chandra[j],y_Xcparts_Chandra[j]+50,'chandra',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.2*subsize,alpha=1.0)
                  a.set_edgecolor('fuchsia')
                  ax1.add_artist(a)
         else:
            if x_Xcparts_Chandra>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Xcparts_Chandra<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Xcparts_Chandra>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Xcparts_Chandra<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
               a = Ellipse(xy=(x_Xcparts_Chandra, y_Xcparts_Chandra), width=2*ell_a_Chandra, height=2*ell_b_Chandra, angle=ell_theta_Chandra)
               a.set_linewidth(0.5)
               a.set_facecolor('none')
               plt.text(x_Xcparts_Chandra,y_Xcparts_Chandra+50,'chandra',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.2*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax1.add_artist(a)
      if os.path.exists(path+'/XMM4d13s.cat'):
         if np.size(x_Xcparts_XMM)!=1:
            for j in range(0,np.size(x_Xcparts_XMM)):
               if x_Xcparts_XMM[j]>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Xcparts_XMM[j]<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Xcparts_XMM[j]>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Xcparts_XMM[j]<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
                  a = Circle(xy=(x_Xcparts_XMM[j], y_Xcparts_XMM[j]), radius=err_Xcparts_XMM[j])
                  a.set_linewidth(0.5)
                  a.set_facecolor('none')
                  plt.text(x_Xcparts_XMM[j],y_Xcparts_XMM[j]+50,'XMM',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.2*subsize,alpha=1.0)
                  a.set_edgecolor('fuchsia')
                  ax1.add_artist(a)
         else:
            if x_Xcparts_XMM>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Xcparts_XMM<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Xcparts_XMM>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Xcparts_XMM<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
               a = Circle(xy=(x_Xcparts_XMM, y_Xcparts_XMM), radius=err_Xcparts_XMM)
               a.set_linewidth(0.5)
               a.set_facecolor('none')
               plt.text(x_Xcparts_XMM,y_Xcparts_XMM+50,'XMM',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.2*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax1.add_artist(a)
      if os.path.exists(path+'/eRASS1d1.cat'):
         if np.size(x_Xcparts_eRASS)!=1:
            for j in range(0,np.size(x_Xcparts_eRASS)):
               if x_Xcparts_eRASS[j]>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Xcparts_eRASS[j]<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Xcparts_eRASS[j]>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Xcparts_eRASS[j]<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
                  a = Circle(xy=(x_Xcparts_eRASS[j], y_Xcparts_eRASS[j]), radius=err_Xcparts_eRASS[j])
                  a.set_linewidth(0.5)
                  a.set_facecolor('none')
                  plt.text(x_Xcparts_eRASS[j],y_Xcparts_eRASS[j]+50,'eRASS',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.2*subsize,alpha=1.0)
                  a.set_edgecolor('fuchsia')
                  ax1.add_artist(a)
         else:
            if x_Xcparts_eRASS>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Xcparts_eRASS<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Xcparts_eRASS>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Xcparts_eRASS<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
               a = Circle(xy=(x_Xcparts_eRASS, y_Xcparts_eRASS), radius=err_Xcparts_eRASS)
               a.set_linewidth(0.5)
               a.set_facecolor('none')
               plt.text(x_Xcparts_eRASS,y_Xcparts_eRASS+50,'eRASS',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.2*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax1.add_artist(a)
      if np.size(x_Rcparts) > 0:
         if np.size(x_Rcparts) > 1:
            for j in range(0,np.size(x_Rcparts)):
               if x_Rcparts[j]>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Rcparts[j]<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Rcparts[j]>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Rcparts[j]<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
                  if not np.isnan(errmaj_Rcparts[j]) and not np.isnan(errmin_Rcparts[j]) and not np.isnan(errangle_Rcparts[j]):
                     a = Ellipse(xy=(x_Rcparts[j], y_Rcparts[j]), width=2*errmaj_Rcparts[j], height=2*errmin_Rcparts[j], angle=errangle_Rcparts[j])
                  elif not np.isnan(errmaj_Rcparts[j]) and not np.isnan(errmin_Rcparts[j]) and np.isnan(errangle_Rcparts[j]):
                     a = Circle(xy=(x_Rcparts[j], y_Rcparts[j]),radius=errmaj_Rcparts[j])
                  else:
                     a = Circle(xy=(x_Rcparts[j], y_Rcparts[j]),radius=2.)
                  a.set_linewidth(0.5)
                  a.set_facecolor('none')
                  plt.text(x_Rcparts[j],y_Rcparts[j]+50,cat_R[j],horizontalalignment='center',verticalalignment='center',color='limegreen',fontsize=1.2*subsize,alpha=1.0)
                  a.set_edgecolor('limegreen')
                  ax1.add_artist(a)
         else:
            if x_Rcparts>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Rcparts<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Rcparts>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Rcparts<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
               if not np.isnan(errmaj_Rcparts) and not np.isnan(errmin_Rcparts) and not np.isnan(errangle_Rcparts):
                  a = Ellipse(xy=(x_Rcparts, y_Rcparts), width=2*errmaj_Rcparts, height=2*errmin_Rcparts, angle=errangle_Rcparts)
               elif not np.isnan(errmaj_Rcparts) and not np.isnan(errmin_Rcparts) and np.isnan(errangle_Rcparts):
                  a = Circle(xy=(x_Rcparts, y_Rcparts),radius=errmaj_Rcparts)
               else:
                  a = Circle(xy=(x_Rcparts, y_Rcparts),radius=2.)
               a.set_linewidth(0.5)
               a.set_facecolor('none')
               plt.text(x_Rcparts,y_Rcparts+50,cat_R,horizontalalignment='center',verticalalignment='center',color='limegreen',fontsize=1.2*subsize,alpha=1.0)
               a.set_edgecolor('limegreen')
               ax1.add_artist(a)

      if 'ra_query' in locals() and 'dec_query' in locals():
         x_query, y_query = w.all_world2pix(ra_query, dec_query, 1)
         plt.text(x_query,y_query+50,"query",horizontalalignment='center',verticalalignment='center',color='darkorange',fontsize=1.2*subsize,alpha=1.0)
         a = Circle(xy=(x_query, y_query), radius=8)
         a.set_linewidth(0.5)
         a.set_facecolor('none')
         a.set_edgecolor('darkorange')
         ax1.add_artist(a)
         

      ax2 = plt.subplot(int(np.trunc(3+n_lc/2)), 2, n_lc+4)
      plt.title('Zoom on the variable source: (RA, DEC) = (%.3f'%ra+', %.3f'%dec+') deg',fontsize=subsize)
      plt.xlim(x-zoom,x+zoom)
      plt.ylim(y-zoom,y+zoom)
      plt.imshow(data, interpolation='nearest', cmap='gray', vmin=m-colscale*s, vmax=m+colscale*s, origin='lower')
      rxnew = [rx[j] for j in range(0,len(rx)) if rx[j]>=(x-zoom) and rx[j]<=(x+zoom) and ry[j]>=(y-zoom) and ry[j]<=(y+zoom)]
      rynew = [ry[j] for j in range(0,len(rx)) if rx[j]>=(x-zoom) and rx[j]<=(x+zoom) and ry[j]>=(y-zoom) and ry[j]<=(y+zoom)]
      for j in range(len(rxnew)):
         if not(np.round(x,ndig)==np.round(rxnew[j],ndig) and np.round(y,ndig)==np.round(rynew[j],ndig)):
            a = Circle(xy=(rxnew[j], rynew[j]), radius=fext*fwhmavg)
            a.set_linewidth(1.0)
            a.set_facecolor('none')
            a.set_edgecolor('blue')
            ax2.add_artist(a)
      plt.text(x,y+20,str(i+1),horizontalalignment='center',verticalalignment='center',color='r',fontsize=2*subsize)
      a = Circle(xy=(x, y), radius=fext*fwhmavg)
      a.set_linewidth(1.0)
      a.set_facecolor('none')
      a.set_edgecolor('red')
      ax2.add_artist(a)
      #if hdr['INSTRUME']=='WIFSIP V2.0':
      #   if np.size(obs_Rcparts)!=0:
      #      xnew_Rcparts = [x_Rcparts[j] for j in range(0,len(x_Rcparts)) if x_Rcparts[j]>=(x-zoom) and x_Rcparts[j]<=(x+zoom) and y_Rcparts[j]>=(y-zoom) and y_Rcparts[j]<=(y+zoom)]
      #      ynew_Rcparts = [y_Rcparts[j] for j in range(0,len(x_Rcparts)) if x_Rcparts[j]>=(x-zoom) and x_Rcparts[j]<=(x+zoom) and y_Rcparts[j]>=(y-zoom) and y_Rcparts[j]<=(y+zoom)]
      #      obsnew_Rcparts = [obs_Rcparts[j] for j in range(0,len(x_Rcparts)) if x_Rcparts[j]>=(x-zoom) and x_Rcparts[j]<=(x+zoom) and y_Rcparts[j]>=(y-zoom) and y_Rcparts[j]<=(y+zoom)]
      #      for j in range(0,np.size(xnew_Rcparts)):
      #         a = Circle(xy=(xnew_Rcparts[j], ynew_Rcparts[j]), radius=8)
      #         a.set_linewidth(1.0)
      #         a.set_facecolor('none')
      #         plt.text(xnew_Rcparts[j],ynew_Rcparts[j],obsnew_Rcparts[j],horizontalalignment='center',verticalalignment='center',color='limegreen',fontsize=1.5*subsize,alpha=1.0)
      #         a.set_edgecolor('limegreen')
      #         ax2.add_artist(a)
      if os.path.exists(path+'/2sxps.cat'):
         if np.size(x_Xcparts_Swift)!=1:
            xnew_Xcparts = [x_Xcparts_Swift[j] for j in range(0,len(x_Xcparts_Swift)) if x_Xcparts_Swift[j]>=(x-zoom) and x_Xcparts_Swift[j]<=(x+zoom) and y_Xcparts_Swift[j]>=(y-zoom) and y_Xcparts_Swift[j]<=(y+zoom)]
            ynew_Xcparts = [y_Xcparts_Swift[j] for j in range(0,len(x_Xcparts_Swift)) if x_Xcparts_Swift[j]>=(x-zoom) and x_Xcparts_Swift[j]<=(x+zoom) and y_Xcparts_Swift[j]>=(y-zoom) and y_Xcparts_Swift[j]<=(y+zoom)]
            errnew_Xcparts = [err_Xcparts_Swift[j] for j in range(0,len(x_Xcparts_Swift)) if x_Xcparts_Swift[j]>=(x-zoom) and x_Xcparts_Swift[j]<=(x+zoom) and y_Xcparts_Swift[j]>=(y-zoom) and y_Xcparts_Swift[j]<=(y+zoom)]
            for j in range(0,np.size(xnew_Xcparts)):
               a = Circle(xy=(xnew_Xcparts[j], ynew_Xcparts[j]), radius=errnew_Xcparts[j])
               a.set_linewidth(1.0)
               a.set_facecolor('none')
               plt.text(xnew_Xcparts[j],ynew_Xcparts[j],'swift',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.5*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax2.add_artist(a)
         else:
            if x_Xcparts_Swift>=(x-zoom) and x_Xcparts_Swift<=(x+zoom) and y_Xcparts_Swift>=(y-zoom) and y_Xcparts_Swift<=(y+zoom):
               xnew_Xcparts = x_Xcparts_Swift
               ynew_Xcparts = y_Xcparts_Swift
               errnew_Xcparts = err_Xcparts_Swift
               a = Circle(xy=(xnew_Xcparts, ynew_Xcparts), radius=errnew_Xcparts)
               a.set_linewidth(1.0)
               a.set_facecolor('none')
               plt.text(xnew_Xcparts,ynew_Xcparts,'swift',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.5*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax2.add_artist(a)
      if os.path.exists(path+'/csc2.cat'):
         if np.size(x_Xcparts_Chandra)!=1:
            xnew_Xcparts = [x_Xcparts_Chandra[j] for j in range(0,len(x_Xcparts_Chandra)) if x_Xcparts_Chandra[j]>=(x-zoom) and x_Xcparts_Chandra[j]<=(x+zoom) and y_Xcparts_Chandra[j]>=(y-zoom) and y_Xcparts_Chandra[j]<=(y+zoom)]
            ynew_Xcparts = [y_Xcparts_Chandra[j] for j in range(0,len(x_Xcparts_Chandra)) if x_Xcparts_Chandra[j]>=(x-zoom) and x_Xcparts_Chandra[j]<=(x+zoom) and y_Xcparts_Chandra[j]>=(y-zoom) and y_Xcparts_Chandra[j]<=(y+zoom)]
            ellanew_Xcparts = [ell_a_Chandra[j] for j in range(0,len(x_Xcparts_Chandra)) if x_Xcparts_Chandra[j]>=(x-zoom) and x_Xcparts_Chandra[j]<=(x+zoom) and y_Xcparts_Chandra[j]>=(y-zoom) and y_Xcparts_Chandra[j]<=(y+zoom)]
            ellbnew_Xcparts = [ell_b_Chandra[j] for j in range(0,len(x_Xcparts_Chandra)) if x_Xcparts_Chandra[j]>=(x-zoom) and x_Xcparts_Chandra[j]<=(x+zoom) and y_Xcparts_Chandra[j]>=(y-zoom) and y_Xcparts_Chandra[j]<=(y+zoom)]
            ellthetanew_Xcparts = [ell_theta_Chandra[j] for j in range(0,len(x_Xcparts_Chandra)) if x_Xcparts_Chandra[j]>=(x-zoom) and x_Xcparts_Chandra[j]<=(x+zoom) and y_Xcparts_Chandra[j]>=(y-zoom) and y_Xcparts_Chandra[j]<=(y+zoom)]
            for j in range(0,np.size(xnew_Xcparts)):
               a = Ellipse(xy=(xnew_Xcparts[j], ynew_Xcparts[j]), width=2*ellanew_Xcparts[j], height=2*ellbnew_Xcparts[j], angle=ellthetanew_Xcparts[j])
               a.set_linewidth(1.0)
               a.set_facecolor('none')
               plt.text(xnew_Xcparts[j],ynew_Xcparts[j],'chandra',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.5*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax2.add_artist(a)
         else:
            if x_Xcparts_Chandra>=(x-zoom) and x_Xcparts_Chandra<=(x+zoom) and y_Xcparts_Chandra>=(y-zoom) and y_Xcparts_Chandra<=(y+zoom):
               xnew_Xcparts = x_Xcparts_Chandra
               ynew_Xcparts = y_Xcparts_Chandra
               ellanew_Xcparts = ell_a_Chandra
               ellbnew_Xcparts = ell_b_Chandra
               ellthetanew_Xcparts = ell_theta_Chandra
               a = Ellipse(xy=(xnew_Xcparts, ynew_Xcparts), width=2*ellanew_Xcparts, height=2*ellbnew_Xcparts, angle=ellthetanew_Xcparts)
               a.set_linewidth(1.0)
               a.set_facecolor('none')
               plt.text(xnew_Xcparts,ynew_Xcparts,'chandra',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.5*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax2.add_artist(a)
      if os.path.exists(path+'/XMM4d13s.cat'):
         if np.size(x_Xcparts_XMM)!=1:
            xnew_Xcparts = [x_Xcparts_XMM[j] for j in range(0,len(x_Xcparts_XMM)) if x_Xcparts_XMM[j]>=(x-zoom) and x_Xcparts_XMM[j]<=(x+zoom) and y_Xcparts_XMM[j]>=(y-zoom) and y_Xcparts_XMM[j]<=(y+zoom)]
            ynew_Xcparts = [y_Xcparts_XMM[j] for j in range(0,len(x_Xcparts_XMM)) if x_Xcparts_XMM[j]>=(x-zoom) and x_Xcparts_XMM[j]<=(x+zoom) and y_Xcparts_XMM[j]>=(y-zoom) and y_Xcparts_XMM[j]<=(y+zoom)]
            errnew_Xcparts = [err_Xcparts_XMM[j] for j in range(0,len(x_Xcparts_XMM)) if x_Xcparts_XMM[j]>=(x-zoom) and x_Xcparts_XMM[j]<=(x+zoom) and y_Xcparts_XMM[j]>=(y-zoom) and y_Xcparts_XMM[j]<=(y+zoom)]
            for j in range(0,np.size(xnew_Xcparts)):
               a = Circle(xy=(xnew_Xcparts[j], ynew_Xcparts[j]), radius=errnew_Xcparts[j])
               a.set_linewidth(1.0)
               a.set_facecolor('none')
               plt.text(xnew_Xcparts[j],ynew_Xcparts[j],'XMM',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.5*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax2.add_artist(a)
         else:
            if x_Xcparts_XMM>=(x-zoom) and x_Xcparts_XMM<=(x+zoom) and y_Xcparts_XMM>=(y-zoom) and y_Xcparts_XMM<=(y+zoom):
               xnew_Xcparts = x_Xcparts_XMM
               ynew_Xcparts = y_Xcparts_XMM
               errnew_Xcparts = err_Xcparts_XMM
               a = Circle(xy=(xnew_Xcparts, ynew_Xcparts), radius=errnew_Xcparts)
               a.set_linewidth(1.0)
               a.set_facecolor('none')
               plt.text(xnew_Xcparts,ynew_Xcparts,'XMM',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.5*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax2.add_artist(a)
      if os.path.exists(path+'/eRASS1d1.cat'):
         if np.size(x_Xcparts_eRASS)!=1:
            xnew_Xcparts = [x_Xcparts_eRASS[j] for j in range(0,len(x_Xcparts_eRASS)) if x_Xcparts_eRASS[j]>=(x-zoom) and x_Xcparts_eRASS[j]<=(x+zoom) and y_Xcparts_eRASS[j]>=(y-zoom) and y_Xcparts_eRASS[j]<=(y+zoom)]
            ynew_Xcparts = [y_Xcparts_eRASS[j] for j in range(0,len(x_Xcparts_eRASS)) if x_Xcparts_eRASS[j]>=(x-zoom) and x_Xcparts_eRASS[j]<=(x+zoom) and y_Xcparts_eRASS[j]>=(y-zoom) and y_Xcparts_eRASS[j]<=(y+zoom)]
            errnew_Xcparts = [err_Xcparts_eRASS[j] for j in range(0,len(x_Xcparts_eRASS)) if x_Xcparts_eRASS[j]>=(x-zoom) and x_Xcparts_eRASS[j]<=(x+zoom) and y_Xcparts_eRASS[j]>=(y-zoom) and y_Xcparts_eRASS[j]<=(y+zoom)]
            for j in range(0,np.size(xnew_Xcparts)):
               a = Circle(xy=(xnew_Xcparts[j], ynew_Xcparts[j]), radius=errnew_Xcparts[j])
               a.set_linewidth(1.0)
               a.set_facecolor('none')
               plt.text(xnew_Xcparts[j],ynew_Xcparts[j],'eRASS',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.5*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax2.add_artist(a)
         else:
            if x_Xcparts_eRASS>=(x-zoom) and x_Xcparts_eRASS<=(x+zoom) and y_Xcparts_eRASS>=(y-zoom) and y_Xcparts_eRASS<=(y+zoom):
               xnew_Xcparts = x_Xcparts_eRASS
               ynew_Xcparts = y_Xcparts_eRASS
               errnew_Xcparts = err_Xcparts_eRASS
               a = Circle(xy=(xnew_Xcparts, ynew_Xcparts), radius=errnew_Xcparts)
               a.set_linewidth(1.0)
               a.set_facecolor('none')
               plt.text(xnew_Xcparts,ynew_Xcparts,'eRASS',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.5*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax2.add_artist(a)
      if np.size(x_Rcparts) > 0:
         if np.size(x_Rcparts) > 1:
            xnew_Rcparts = [x_Rcparts[j] for j in range(0,len(x_Rcparts)) if x_Rcparts[j]>=(x-zoom) and x_Rcparts[j]<=(x+zoom) and y_Rcparts[j]>=(y-zoom) and y_Rcparts[j]<=(y+zoom)]
            ynew_Rcparts = [y_Rcparts[j] for j in range(0,len(x_Rcparts)) if x_Rcparts[j]>=(x-zoom) and x_Rcparts[j]<=(x+zoom) and y_Rcparts[j]>=(y-zoom) and y_Rcparts[j]<=(y+zoom)]
            errmajnew_Rcparts = [errmaj_Rcparts[j] for j in range(0,len(x_Rcparts)) if x_Rcparts[j]>=(x-zoom) and x_Rcparts[j]<=(x+zoom) and y_Rcparts[j]>=(y-zoom) and y_Rcparts[j]<=(y+zoom)]
            errminnew_Rcparts = [errmin_Rcparts[j] for j in range(0,len(x_Rcparts)) if x_Rcparts[j]>=(x-zoom) and x_Rcparts[j]<=(x+zoom) and y_Rcparts[j]>=(y-zoom) and y_Rcparts[j]<=(y+zoom)]
            erranglenew_Rcparts = [errangle_Rcparts[j] for j in range(0,len(x_Rcparts)) if x_Rcparts[j]>=(x-zoom) and x_Rcparts[j]<=(x+zoom) and y_Rcparts[j]>=(y-zoom) and y_Rcparts[j]<=(y+zoom)]
            catnew_R = [cat_R[j] for j in range(0,len(x_Rcparts)) if x_Rcparts[j]>=(x-zoom) and x_Rcparts[j]<=(x+zoom) and y_Rcparts[j]>=(y-zoom) and y_Rcparts[j]<=(y+zoom)]
            for j in range(0,np.size(xnew_Rcparts)):
               if not np.isnan(errmajnew_Rcparts[j]) and not np.isnan(errminnew_Rcparts[j]) and not np.isnan(erranglenew_Rcparts[j]):
                  a = Ellipse(xy=(xnew_Rcparts[j], ynew_Rcparts[j]), width=2*errmajnew_Rcparts[j], height=2*errminnew_Rcparts[j], angle=erranglenew_Rcparts[j])
               elif not np.isnan(errmajnew_Rcparts[j]) and not np.isnan(errminnew_Rcparts[j]) and np.isnan(erranglenew_Rcparts[j]):
                  a = Circle(xy=(xnew_Rcparts[j], ynew_Rcparts[j]),radius=errmajnew_Rcparts[j])
               else:
                  a = Circle(xy=(xnew_Rcparts[j], ynew_Rcparts[j]),radius=2.)
               a.set_linewidth(1.0)
               a.set_facecolor('none')
               plt.text(xnew_Rcparts[j],ynew_Rcparts[j],catnew_R[j],horizontalalignment='center',verticalalignment='center',color='limegreen',fontsize=1.5*subsize,alpha=1.0)
               a.set_edgecolor('limegreen')
               ax2.add_artist(a)
         else:
            if x_Rcparts>=(x-zoom) and x_Rcparts<=(x+zoom) and y_Rcparts>=(y-zoom) and y_Rcparts<=(y+zoom):
               xnew_Rcparts = x_Rcparts
               ynew_Rcparts = y_Rcparts
               errmajnew_Rcparts = errmaj_Rcparts
               errminnew_Rcparts = errmin_Rcparts
               erranglenew_Rcparts = errangle_Rcparts
               catnew_R = cat_R
               if not np.isnan(errmajnew_Rcparts) and not np.isnan(errminnew_Rcparts) and not np.isnan(erranglenew_Rcparts):
                  a = Ellipse(xy=(xnew_Rcparts, ynew_Rcparts), width=2*errmajnew_Rcparts, height=2*errminnew_Rcparts, angle=erranglenew_Rcparts)
               elif not np.isnan(errmajnew_Rcparts) and not np.isnan(errminnew_Rcparts) and np.isnan(erranglenew_Rcparts):
                  a = Circle(xy=(xnew_Rcparts, ynew_Rcparts),radius=errmajnew_Rcparts)
               else:
                  a = Circle(xy=(xnew_Rcparts, ynew_Rcparts),radius=2.)
               a.set_linewidth(1.0)
               a.set_facecolor('none')
               plt.text(xnew_Rcparts,ynew_Rcparts,catnew_R,horizontalalignment='center',verticalalignment='center',color='limegreen',fontsize=1.5*subsize,alpha=1.0)
               a.set_edgecolor('limegreen')
               ax2.add_artist(a)

      if 'ra_query' in locals() and 'dec_query' in locals():
         plt.text(x_query,y_query,"query",horizontalalignment='center',verticalalignment='center',color='darkorange',fontsize=1.5*subsize,alpha=0.8)
         a = Circle(xy=(x_query, y_query), radius=8)
         a.set_linewidth(0.5)
         a.set_facecolor('none')
         a.set_edgecolor('darkorange')
         ax2.add_artist(a)
         
      if 'kic_ID' in locals():
         plt.text(-0.07,4.7,"Identified from Kepler as "+kic_ID+r" with $P_{\mathrm{Kepler}}$="+str(kic_period)+" d",transform=ax2.transAxes,fontsize=tsize)
      elif 'ATO_ID' in locals():
         plt.text(-0.07,4.7,"Identified from ATLAS as "+ATO_ID+r" with $P_{\mathrm{ATLAS}}$="+str(ATO_period)+" d",transform=ax2.transAxes,fontsize=tsize)
      elif 'SIMBAD_ID' in locals():
         plt.text(-0.07,4.7,"Identified from SIMBAD as "+SIMBAD_ID,transform=ax2.transAxes,fontsize=tsize)
      if per_search == 'LS':
         plt.savefig(path+'/finalperiodplots'+str(nlist)+'_'+per_search+'/finalperplot_'+str(i+1)+'_'+per_search+'_'+analysis+'_nper'+str(int(n_periods_sub))+'.png',bbox_inches='tight')
      else:
         plt.savefig(path+'/finalperiodplots'+str(nlist)+'_'+per_search+'/finalperplot_'+str(i+1)+'_PDM_'+analysis+'_nper'+str(int(n_periods_sub))+'.png',bbox_inches='tight')
      plt.close(figure)

   #Download ZTF (Zwicky transient catalog) light curves in g and r bands of the closest object within 0.001 degrees
   #df = lightcurve.LCQuery.download_data(circle=[ra, dec,0.001], bandname="g")
   if "df['oid']" in locals():
      g_check = 1
      ztf_oid = df['oid'].to_numpy()
      ztf_mjd = df['mjd'].to_numpy()
      ztf_mag = df['mag'].to_numpy()
      ztf_magerr = df['magerr'].to_numpy()
      ztf_ra = df['ra'].to_numpy()
      ztf_dec = df['dec'].to_numpy()
      ztf_airmass = df['airmass'].to_numpy()
      list_id = []
      list_eqcoord = []
      for j in range(0,len(ztf_oid)):
         if j==0:
            list_id.append(ztf_oid[j])
            list_eqcoord.append([ztf_ra[j],ztf_dec[j]])
         else:
            if ztf_oid[j]!=ztf_oid[j-1]:
               list_id.append(ztf_oid[j])
               list_eqcoord.append([ztf_ra[j],ztf_dec[j]])
      list_eqcoord = np.asarray(list_eqcoord)
      sep = SkyCoord(ra, dec, frame='fk5', unit='deg').separation(SkyCoord(list_eqcoord[:,0],list_eqcoord[:,1],frame='fk5',unit='deg')).arcsecond
      ztf_mjd_g = np.array([ztf_mjd[j] for j in range(0,len(ztf_oid)) if ztf_oid[j] == ztf_oid[np.argmin(sep)] and ztf_airmass[j]<=2.])
      ztf_mag_g = np.array([ztf_mag[j] for j in range(0,len(ztf_oid)) if ztf_oid[j] == ztf_oid[np.argmin(sep)] and ztf_airmass[j]<=2.])
      ztf_magerr_g = np.array([ztf_magerr[j] for j in range(0,len(ztf_oid)) if ztf_oid[j] == ztf_oid[np.argmin(sep)] and ztf_airmass[j]<=2.])
      #ztf_mjd_g = np.array([ztf_mjd_g[j] for j in range(0,len(ztf_mag_g)) if ztf_mag_g[j] <= 17.45])
      #ztf_magerr_g = np.array([ztf_magerr_g[j] for j in range(0,len(ztf_mag_g)) if ztf_mag_g[j] <= 17.45])
      #ztf_mag_g = np.array([ztf_mag_g[j] for j in range(0,len(ztf_mag_g)) if ztf_mag_g[j] <= 17.45])
   else:
      g_check = 0

   #df = lightcurve.LCQuery.download_data(circle=[ra, dec,0.001], bandname="r")
   if "df['oid']" in locals():
      r_check = 1
      ztf_oid = df['oid'].to_numpy()
      ztf_mjd = df['mjd'].to_numpy()
      ztf_mag = df['mag'].to_numpy()
      ztf_magerr = df['magerr'].to_numpy()
      ztf_ra = df['ra'].to_numpy()
      ztf_dec = df['dec'].to_numpy()
      ztf_airmass = df['airmass'].to_numpy()
      list_id = []
      list_eqcoord = []
      for j in range(0,len(ztf_oid)):
         if j==0:
            list_id.append(ztf_oid[j])
            list_eqcoord.append([ztf_ra[j],ztf_dec[j]])
         else:
            if ztf_oid[j]!=ztf_oid[j-1]:
               list_id.append(ztf_oid[j])
               list_eqcoord.append([ztf_ra[j],ztf_dec[j]])
      list_eqcoord = np.asarray(list_eqcoord)
      sep = SkyCoord(ra, dec, frame='fk5', unit='deg').separation(SkyCoord(list_eqcoord[:,0],list_eqcoord[:,1],frame='fk5',unit='deg')).arcsecond
      ztf_x, ztf_y = w.all_world2pix(list_eqcoord[np.argmin(sep),0], list_eqcoord[np.argmin(sep),1], 1)
      ztf_mjd_r = np.array([ztf_mjd[j] for j in range(0,len(ztf_oid)) if ztf_oid[j] == ztf_oid[np.argmin(sep)] and ztf_airmass[j]<=2.])
      ztf_mag_r = np.array([ztf_mag[j] for j in range(0,len(ztf_oid)) if ztf_oid[j] == ztf_oid[np.argmin(sep)] and ztf_airmass[j]<=2.])
      ztf_magerr_r = np.array([ztf_magerr[j] for j in range(0,len(ztf_oid)) if ztf_oid[j] == ztf_oid[np.argmin(sep)] and ztf_airmass[j]<=2.])
      #ztf_mjd_r = np.array([ztf_mjd_r[j] for j in range(0,len(ztf_mag_r)) if ztf_mag_r[j] <= 16.2])
      #ztf_magerr_r = np.array([ztf_magerr_r[j] for j in range(0,len(ztf_mag_r)) if ztf_mag_r[j] <= 16.2])
      #ztf_mag_r = np.array([ztf_mag_r[j] for j in range(0,len(ztf_mag_r)) if ztf_mag_r[j] <= 16.2])
      sep = np.min(sep)
   else:
      r_check = 0
   if g_check == 1 or r_check==1:
      del df, ztf_oid, ztf_mjd, ztf_mag, ztf_magerr, ztf_ra, ztf_dec, ztf_airmass

   if g_check == 1:
      print("ZTF query in the g band succedeed")
   else:
      print("No ZTF data in the g band")
   if r_check == 1:
      print("ZTF query in the r band succedeed")
   else:
      print("No ZTF data in the r band")
   
   #Plot ZTF lightcurves and color curves folded at best_periods[i], best_periods[i]/2, 2*best_periods[0] and t_0 in g and r bands together with our light curves, with also a link for the SIMBAD query (in pdf format) and the period from the eventual corresponding ATLAS variable source
   if r_check == 1:
      ZTF_std = np.std(ztf_mag_r)
      ZTF_mjd_start = ztf_mjd_r[0]
      ZTF_mjd_end = ztf_mjd_r[-1]
   else:
      if g_check == 1:
         ZTF_std = np.std(ztf_mag_g)
         ZTF_mjd_start = ztf_mjd_g[0]
         ZTF_mjd_end = ztf_mjd_g[-1]
      else:
         ZTF_std = 0.
   
   for j in range(0,2):
      figure=plt.figure(figsize=(figxsize,figysize), dpi=150)
      if np.size(sources_3FGL)!=0:
         if 'flux_Xcpart' in locals():
            figure.suptitle(r"Fermi source: "+FoVname+", Periodic variable #%d"%(i+1)+", $\sigma_{\mathrm{ZTF}}$: %.2f"%ZTF_std+r" mag, $\sigma_{\mathrm{"+telescope.split("/")[0]+"}}$: %.2f"%sigma+" mag\n\nPeriod folding with "+r"$P_{\mathrm{best}}=$"+"{0} d, ".format(best_periods[i])+"$T_{0}=$%.4f"%t_0+" MJD\n\n%.1f"%(sep_4FGL/ell_a_68)+r"$\sigma$ far from the 4FGL source and %.1f"%(sep_3FGL/ell_a_68_3FGL)+r"$\sigma$ far from the 3FGL source, X-ray flux: "+str(np.round(float(str(flux_Xcpart).split('e')[0]),1))+"$\cdot 10^{"+str(flux_Xcpart).split('e')[1]+"}$ erg s$^{-1}$ cm$^{-2}$", fontsize=tsize)
         else:
            if 'UL_X' in locals():
               figure.suptitle(r"Fermi source: "+FoVname+", Periodic variable #%d"%(i+1)+", $\sigma_{\mathrm{ZTF}}$: %.2f"%ZTF_std+r" mag, $\sigma_{\mathrm{"+telescope.split("/")[0]+"}}$: %.2f"%sigma+" mag\n\nPeriod folding with "+r"$P_{\mathrm{best}}=$"+"{0} d, ".format(best_periods[i])+"$T_{0}=$%.4f"%t_0+" MJD\n\n%.1f"%(sep_4FGL/ell_a_68)+r"$\sigma$ far from the 4FGL source and %.1f"%(sep_3FGL/ell_a_68_3FGL)+r"$\sigma$ far from the 3FGL source, X-ray flux 3$\sigma$ UL: "+str(np.round(float(str(UL_X).split('e')[0]),1))+"$\cdot 10^{"+str(UL_X).split('e')[1]+"}$ erg s$^{-1}$ cm$^{-2}$", fontsize=tsize)
            else:
               figure.suptitle(r"Fermi source: "+FoVname+", Periodic variable #%d"%(i+1)+", $\sigma_{\mathrm{ZTF}}$: %.2f"%ZTF_std+r" mag, $\sigma_{\mathrm{"+telescope.split("/")[0]+"}}$: %.2f"%sigma+" mag\n\nPeriod folding with "+r"$P_{\mathrm{best}}=$"+"{0} d, ".format(best_periods[i])+"$T_{0}=$%.4f"%t_0+" MJD\n\n%.1f"%(sep_4FGL/ell_a_68)+r"$\sigma$ far from the 4FGL source and %.1f"%(sep_3FGL/ell_a_68_3FGL)+r"$\sigma$ far from the 3FGL source", fontsize=tsize)
      else:
         if 'flux_Xcpart' in locals():
            figure.suptitle(r"Fermi source: "+FoVname+", Periodic variable #%d"%(i+1)+", $\sigma_{\mathrm{ZTF}}$: %.2f"%ZTF_std+r" mag, $\sigma_{\mathrm{"+telescope.split("/")[0]+"}}$: %.2f"%sigma+" mag\n\nPeriod folding with "+r"$P_{\mathrm{best}}=$"+"{0} d, ".format(best_periods[i])+"$T_{0}=$%.4f"%t_0+" MJD\n\n%.1f"%(sep_4FGL/ell_a_68)+r"$\sigma$ far from the 4FGL source, X-ray flux: "+str(np.round(float(str(flux_Xcpart).split('e')[0]),1))+"$\cdot 10^{"+str(flux_Xcpart).split('e')[1]+"}$ erg s$^{-1}$ cm$^{-2}$", fontsize=tsize)
         else:
            if 'UL_X' in locals():
               figure.suptitle(r"Fermi source: "+FoVname+", Periodic variable #%d"%(i+1)+", $\sigma_{\mathrm{ZTF}}$: %.2f"%ZTF_std+r" mag, $\sigma_{\mathrm{"+telescope.split("/")[0]+"}}$: %.2f"%sigma+" mag\n\nPeriod folding with "+r"$P_{\mathrm{best}}=$"+"{0} d, ".format(best_periods[i])+"$T_{0}=$%.4f"%t_0+" MJD\n\n%.1f"%(sep_4FGL/ell_a_68)+r"$\sigma$ far from the 4FGL source, X-ray flux 3$\sigma$ UL: "+str(np.round(float(str(UL_X).split('e')[0]),1))+"$\cdot 10^{"+str(UL_X).split('e')[1]+"}$ erg s$^{-1}$ cm$^{-2}$", fontsize=tsize)
            else:
               figure.suptitle(r"Fermi source: "+FoVname+", Periodic variable #%d"%(i+1)+", $\sigma_{\mathrm{ZTF}}$: %.2f"%ZTF_std+r" mag, $\sigma_{\mathrm{"+telescope.split("/")[0]+"}}$: %.2f"%sigma+" mag\n\nPeriod folding with "+r"$P_{\mathrm{best}}=$"+"{0} d, ".format(best_periods[i])+"$T_{0}=$%.4f"%t_0+" MJD\n\n%.1f"%(sep_4FGL/ell_a_68)+r"$\sigma$ far from the 4FGL source", fontsize=tsize)

      if g_check == 1:
         ztf_phase_g = ((ztf_mjd_g-t_0)/best_periods[i]) % 1
         plt.subplot(5, 2, 1)
         plt.plot(ztf_phase_g, ztf_mag_g, 'bo', linestyle='None')
         plt.plot(ztf_phase_g+1, ztf_mag_g, 'ro', linestyle='None')
         plt.errorbar(ztf_phase_g, ztf_mag_g, ztf_magerr_g, linestyle='None')
         plt.errorbar(ztf_phase_g+1, ztf_mag_g, ztf_magerr_g, linestyle='None')
         plt.xticks(fontsize=axticksize)
         plt.yticks(fontsize=axticksize)
         plt.gca().invert_yaxis()
         plt.title(r"Lightcurve in the g band folded at $P_{\mathrm{best}}$: ZTF data",fontsize=subsize)
         plt.xlabel(r"Phase ($\phi$)",fontsize=axlabelsize)
         plt.ylabel(r"Apparent magnitude",fontsize=axlabelsize)

      if os.path.exists(path+'/gfilter/periods/V'+str(index)+'_LombScargle_N1_Lightcurve.csv'):
         phase_g = ((mjd_g-t_0)/best_periods[i]) % 1
         plt.subplot(5, 2, 2)
         plt.plot(phase_g, m_g, 'bo', linestyle='None')
         plt.plot(phase_g+1, m_g, 'ro', linestyle='None')
         plt.errorbar(phase_g, m_g, sigmam_g, linestyle='None')
         plt.errorbar(phase_g+1, m_g, sigmam_g, linestyle='None')
         plt.xticks(fontsize=axticksize)
         plt.yticks(fontsize=axticksize)
         plt.gca().invert_yaxis()
         plt.title(r"Lightcurve in the g band folded at $P_{\mathrm{best}}$: "+telescope.split("/")[0]+" data",fontsize=subsize)
         plt.xlabel(r"Phase ($\phi$)",fontsize=axlabelsize)
         plt.ylabel(r"Apparent magnitude",fontsize=axlabelsize)

      if r_check == 1:
         ztf_phase_r = ((ztf_mjd_r-t_0)/best_periods[i]) % 1
         plt.subplot(5, 2, 3)
         plt.plot(ztf_phase_r, ztf_mag_r, 'bo', linestyle='None')
         plt.plot(ztf_phase_r+1, ztf_mag_r, 'ro', linestyle='None')
         plt.errorbar(ztf_phase_r, ztf_mag_r, ztf_magerr_r, linestyle='None')
         plt.errorbar(ztf_phase_r+1, ztf_mag_r, ztf_magerr_r, linestyle='None')
         plt.xticks(fontsize=axticksize)
         plt.yticks(fontsize=axticksize)
         plt.gca().invert_yaxis()
         plt.title(r"Lightcurve in the r band folded at $P_{\mathrm{best}}$: ZTF data",fontsize=subsize)
         plt.xlabel(r"Phase ($\phi$)",fontsize=axlabelsize)
         plt.ylabel(r"Apparent magnitude",fontsize=axlabelsize)

      if os.path.exists(path+'/rfilter/periods/V'+str(index)+'_LombScargle_N1_Lightcurve.csv'):
         phase_r = ((mjd_r-t_0)/best_periods[i]) % 1
         plt.subplot(5, 2, 4)
         plt.plot(phase_r, m_r, 'bo', linestyle='None')
         plt.plot(phase_r+1, m_r, 'ro', linestyle='None')
         plt.errorbar(phase_r, m_r, sigmam_r, linestyle='None')
         plt.errorbar(phase_r+1, m_r, sigmam_r, linestyle='None')
         plt.xticks(fontsize=axticksize)
         plt.yticks(fontsize=axticksize)
         plt.gca().invert_yaxis()
         plt.title(r"Lightcurve in the r band folded at $P_{\mathrm{best}}$: "+telescope.split("/")[0]+" data",fontsize=subsize)
         plt.xlabel(r"Phase ($\phi$)",fontsize=axlabelsize)
         plt.ylabel(r"Apparent magnitude",fontsize=axlabelsize)

      if r_check == 1:
         ztf_phase_r = ((ztf_mjd_r-t_0)/(0.5*best_periods[i])) % 1
         plt.subplot(5, 2, 5)
         plt.plot(ztf_phase_r, ztf_mag_r, 'bo', linestyle='None')
         plt.plot(ztf_phase_r+1, ztf_mag_r, 'ro', linestyle='None')
         plt.errorbar(ztf_phase_r, ztf_mag_r, ztf_magerr_r, linestyle='None')
         plt.errorbar(ztf_phase_r+1, ztf_mag_r, ztf_magerr_r, linestyle='None')
         plt.xticks(fontsize=axticksize)
         plt.yticks(fontsize=axticksize)
         plt.gca().invert_yaxis()
         plt.title(r"Lightcurve in the r band folded at $P_{\mathrm{best}}/2$: ZTF data",fontsize=subsize)
         plt.xlabel(r"Phase ($\phi$)",fontsize=axlabelsize)
         plt.ylabel(r"Apparent magnitude",fontsize=axlabelsize)

      if os.path.exists(path+'/rfilter/periods/V'+str(index)+'_LombScargle_N1_Lightcurve.csv'):
         phase_r = ((mjd_r-t_0)/(0.5*best_periods[i])) % 1
         plt.subplot(5, 2, 6)
         plt.plot(phase_r, m_r, 'bo', linestyle='None')
         plt.plot(phase_r+1, m_r, 'ro', linestyle='None')
         plt.errorbar(phase_r, m_r, sigmam_r, linestyle='None')
         plt.errorbar(phase_r+1, m_r, sigmam_r, linestyle='None')
         plt.xticks(fontsize=axticksize)
         plt.yticks(fontsize=axticksize)
         plt.gca().invert_yaxis()
         plt.title(r"Lightcurve in the r band folded at $P_{\mathrm{best}}/2$: "+telescope.split("/")[0]+" data",fontsize=subsize)
         plt.xlabel(r"Phase ($\phi$)",fontsize=axlabelsize)
         plt.ylabel(r"Apparent magnitude",fontsize=axlabelsize)

      if r_check == 1:
         ztf_phase_r = ((ztf_mjd_r-t_0)/(2*best_periods[i])) % 1
         plt.subplot(5, 2, 7)
         plt.plot(ztf_phase_r, ztf_mag_r, 'bo', linestyle='None')
         plt.plot(ztf_phase_r+1, ztf_mag_r, 'ro', linestyle='None')
         plt.errorbar(ztf_phase_r, ztf_mag_r, ztf_magerr_r, linestyle='None')
         plt.errorbar(ztf_phase_r+1, ztf_mag_r, ztf_magerr_r, linestyle='None')
         plt.xticks(fontsize=axticksize)
         plt.yticks(fontsize=axticksize)
         plt.gca().invert_yaxis()
         plt.title(r"Lightcurve in the r band folded at $2P_{\mathrm{best}}$: ZTF data",fontsize=subsize)
         plt.xlabel(r"Phase ($\phi$)",fontsize=axlabelsize)
         plt.ylabel(r"Apparent magnitude",fontsize=axlabelsize)

      if os.path.exists(path+'/rfilter/periods/V'+str(index)+'_LombScargle_N1_Lightcurve.csv'):
         phase_r = ((mjd_r-t_0)/(2*best_periods[i])) % 1
         plt.subplot(5, 2, 8)
         plt.plot(phase_r, m_r, 'bo', linestyle='None')
         plt.plot(phase_r+1, m_r, 'ro', linestyle='None')
         plt.errorbar(phase_r, m_r, sigmam_r, linestyle='None')
         plt.errorbar(phase_r+1, m_r, sigmam_r, linestyle='None')
         plt.xticks(fontsize=axticksize)
         plt.yticks(fontsize=axticksize)
         plt.gca().invert_yaxis()
         plt.title(r"Lightcurve in the r band folded at $2P_{\mathrm{best}}$: "+telescope.split("/")[0]+" data",fontsize=subsize)
         plt.xlabel(r"Phase ($\phi$)",fontsize=axlabelsize)
         plt.ylabel(r"Apparent magnitude",fontsize=axlabelsize)
   
      ax1 = plt.subplot(5, 2, 9)
      plt.title(r'Median r FoV bkg subtracted in the '+str(2*bound)+'$\sigma$ square around the 4FGL coordinates',fontsize=subsize)
      plt.xlim(max(0,xtar-bound*max(ell_a_95,ell_b_95)), min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)))
      plt.ylim(max(0,ytar-bound*max(ell_a_95,ell_b_95)), min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)))
      plt.imshow(data, interpolation='nearest', cmap='gray', vmin=m-colscale*s, vmax=m+colscale*s, origin='lower')
      plt.text(xtar,ytar,FoVname,horizontalalignment='center',verticalalignment='center',color='cyan',fontsize=subsize,alpha=0.8)
      target = Ellipse(xy=(xtar, ytar), width=2*ell_a_95, height=2*ell_b_95, angle=ell_theta_95)
      target.set_facecolor('none')
      target.set_edgecolor('cyan')
      ax1.add_artist(target)
      if np.size(sources_3FGL)!=0:
         for k in range(0,np.shape(sources_3FGL)[0]):
            if x_3FGL[k]>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_3FGL[k]<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_3FGL[k]>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_3FGL[k]<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
               plt.text(x_3FGL[k],y_3FGL[k],sources_3FGL[k].replace("L ","L"),horizontalalignment='center',verticalalignment='center',color='yellow',fontsize=subsize,alpha=0.8)
               a = Ellipse(xy=(x_3FGL[k], y_3FGL[k]), width=2*3600./float(platescale)*coord_3FGL[k,3], height=2*3600./float(platescale)*coord_3FGL[k,4], angle=coord_3FGL[k,5])
               a.set_facecolor('none')
               a.set_edgecolor('yellow')
               ax1.add_artist(a)
      #if np.size(sources_4FGLe)!=0:
      #   for k in range(0,np.shape(sources_4FGLe)[0]):
      #      if x_4FGLe[k]>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_4FGLe[k]<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_4FGLe[k]>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_4FGLe[k]<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
      #         plt.text(x_4FGLe[k],y_4FGLe[k],sources_4FGLe[k].replace("L ","L"),horizontalalignment='center',verticalalignment='center',color='cyan',fontsize=subsize,alpha=0.8)
      #         a = Ellipse(xy=(x_4FGLe[k], y_4FGLe[k]), width=2*3600./float(platescale)*np.asarray(coord_4FGLe)[k,2], height=2*3600./float(platescale)*np.asarray(coord_4FGLe)[k,3], angle=np.asarray(coord_4FGLe)[k,4])
      #         a.set_facecolor('none')
      #         a.set_edgecolor('cyan')
      #         ax1.add_artist(a)
         
      x, y = w.all_world2pix(ra, dec, 1)
      plt.text(x,y+50,str(i+1),horizontalalignment='center',verticalalignment='center',color='r',fontsize=subsize)
      a = Circle(xy=(x, y), radius=fext*fwhmavg)
      a.set_linewidth(0.5)
      a.set_facecolor('none')
      a.set_edgecolor('red')
      ax1.add_artist(a)
      #if hdr['INSTRUME']=='WIFSIP V2.0':
      #   if np.size(obs_Rcparts)!=0:
      #      for k in range(0,np.size(obs_Rcparts)):
      #         if x_Rcparts[k]>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Rcparts[k]<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Rcparts[k]>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Rcparts[k]<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
      #            a = Circle(xy=(x_Rcparts[k], y_Rcparts[k]), radius=8)
      #            a.set_linewidth(0.5)
      #            a.set_facecolor('none')
      #            plt.text(x_Rcparts[k],y_Rcparts[k]+50,obs_Rcparts[k],horizontalalignment='center',verticalalignment='center',color='limegreen',fontsize=1.2*subsize,alpha=0.8)
      #            a.set_edgecolor('limegreen')
      #            ax1.add_artist(a)

      if os.path.exists(path+'/2sxps.cat'):
         if np.size(x_Xcparts_Swift)!=1:
            for j in range(0,np.size(x_Xcparts_Swift)):
               if x_Xcparts_Swift[j]>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Xcparts_Swift[j]<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Xcparts_Swift[j]>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Xcparts_Swift[j]<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
                  a = Circle(xy=(x_Xcparts_Swift[j], y_Xcparts_Swift[j]), radius=err_Xcparts_Swift[j])
                  a.set_linewidth(0.5)
                  a.set_facecolor('none')
                  plt.text(x_Xcparts_Swift[j],y_Xcparts_Swift[j]+50,'swift',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.2*subsize,alpha=1.0)
                  a.set_edgecolor('fuchsia')
                  ax1.add_artist(a)
         else:
            if x_Xcparts_Swift>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Xcparts_Swift<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Xcparts_Swift>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Xcparts_Swift<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
               a = Circle(xy=(x_Xcparts_Swift, y_Xcparts_Swift), radius=err_Xcparts_Swift)
               a.set_linewidth(0.5)
               a.set_facecolor('none')
               plt.text(x_Xcparts_Swift,y_Xcparts_Swift+50,'swift',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.2*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax1.add_artist(a)
      if os.path.exists(path+'/csc2.cat'):
         if np.size(x_Xcparts_Chandra)!=1:
            for j in range(0,np.size(x_Xcparts_Chandra)):
               if x_Xcparts_Chandra[j]>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Xcparts_Chandra[j]<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Xcparts_Chandra[j]>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Xcparts_Chandra[j]<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
                  a = Ellipse(xy=(x_Xcparts_Chandra[j], y_Xcparts_Chandra[j]), width=2*ell_a_Chandra[j], height=2*ell_b_Chandra[j], angle=ell_theta_Chandra[j])
                  a.set_linewidth(0.5)
                  a.set_facecolor('none')
                  plt.text(x_Xcparts_Chandra[j],y_Xcparts_Chandra[j]+50,'chandra',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.2*subsize,alpha=1.0)
                  a.set_edgecolor('fuchsia')
                  ax1.add_artist(a)
         else:
            if x_Xcparts_Chandra>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Xcparts_Chandra<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Xcparts_Chandra>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Xcparts_Chandra<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
               a = Ellipse(xy=(x_Xcparts_Chandra, y_Xcparts_Chandra), width=2*ell_a_Chandra, height=2*ell_b_Chandra, angle=ell_theta_Chandra)
               a.set_linewidth(0.5)
               a.set_facecolor('none')
               plt.text(x_Xcparts_Chandra,y_Xcparts_Chandra+50,'chandra',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.2*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax1.add_artist(a)
      if os.path.exists(path+'/XMM4d13s.cat'):
         if np.size(x_Xcparts_XMM)!=1:
            for j in range(0,np.size(x_Xcparts_XMM)):
               if x_Xcparts_XMM[j]>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Xcparts_XMM[j]<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Xcparts_XMM[j]>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Xcparts_XMM[j]<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
                  a = Circle(xy=(x_Xcparts_XMM[j], y_Xcparts_XMM[j]), radius=err_Xcparts_XMM[j])
                  a.set_linewidth(0.5)
                  a.set_facecolor('none')
                  plt.text(x_Xcparts_XMM[j],y_Xcparts_XMM[j]+50,'XMM',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.2*subsize,alpha=1.0)
                  a.set_edgecolor('fuchsia')
                  ax1.add_artist(a)
         else:
            if x_Xcparts_XMM>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Xcparts_XMM<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Xcparts_XMM>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Xcparts_XMM<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
               a = Circle(xy=(x_Xcparts_XMM, y_Xcparts_XMM), radius=err_Xcparts_XMM)
               a.set_linewidth(0.5)
               a.set_facecolor('none')
               plt.text(x_Xcparts_XMM,y_Xcparts_XMM+50,'XMM',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.2*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax1.add_artist(a)
      if os.path.exists(path+'/eRASS1d1.cat'):
         if np.size(x_Xcparts_eRASS)!=1:
            for j in range(0,np.size(x_Xcparts_eRASS)):
               if x_Xcparts_eRASS[j]>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Xcparts_eRASS[j]<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Xcparts_eRASS[j]>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Xcparts_eRASS[j]<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
                  a = Circle(xy=(x_Xcparts_eRASS[j], y_Xcparts_eRASS[j]), radius=err_Xcparts_eRASS[j])
                  a.set_linewidth(0.5)
                  a.set_facecolor('none')
                  plt.text(x_Xcparts_eRASS[j],y_Xcparts_eRASS[j]+50,'eRASS',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.2*subsize,alpha=1.0)
                  a.set_edgecolor('fuchsia')
                  ax1.add_artist(a)
         else:
            if x_Xcparts_eRASS>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Xcparts_eRASS<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Xcparts_eRASS>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Xcparts_eRASS<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
               a = Circle(xy=(x_Xcparts_eRASS, y_Xcparts_eRASS), radius=err_Xcparts_eRASS)
               a.set_linewidth(0.5)
               a.set_facecolor('none')
               plt.text(x_Xcparts_eRASS,y_Xcparts_eRASS+50,'eRASS',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.2*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax1.add_artist(a)
      if np.size(x_Rcparts) > 0:
         if np.size(x_Rcparts) > 1:
            for j in range(0,np.size(x_Rcparts)):
               if x_Rcparts[j]>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Rcparts[j]<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Rcparts[j]>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Rcparts[j]<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
                  if not np.isnan(errmaj_Rcparts[j]) and not np.isnan(errmin_Rcparts[j]) and not np.isnan(errangle_Rcparts[j]):
                     a = Ellipse(xy=(x_Rcparts[j], y_Rcparts[j]), width=2*errmaj_Rcparts[j], height=2*errmin_Rcparts[j], angle=errangle_Rcparts[j])
                  elif not np.isnan(errmaj_Rcparts[j]) and not np.isnan(errmin_Rcparts[j]) and np.isnan(errangle_Rcparts[j]):
                     a = Circle(xy=(x_Rcparts[j], y_Rcparts[j]),radius=errmaj_Rcparts[j])
                  else:
                     a = Circle(xy=(x_Rcparts[j], y_Rcparts[j]),radius=2.)
                  a.set_linewidth(0.5)
                  a.set_facecolor('none')
                  plt.text(x_Rcparts[j],y_Rcparts[j]+50,cat_R[j],horizontalalignment='center',verticalalignment='center',color='limegreen',fontsize=1.2*subsize,alpha=1.0)
                  a.set_edgecolor('limegreen')
                  ax1.add_artist(a)
         else:
            if x_Rcparts>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Rcparts<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Rcparts>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Rcparts<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
               if not np.isnan(errmaj_Rcparts) and not np.isnan(errmin_Rcparts) and not np.isnan(errangle_Rcparts):
                  a = Ellipse(xy=(x_Rcparts, y_Rcparts), width=2*errmaj_Rcparts, height=2*errmin_Rcparts, angle=errangle_Rcparts)
               elif not np.isnan(errmaj_Rcparts) and not np.isnan(errmin_Rcparts) and np.isnan(errangle_Rcparts):
                  a = Circle(xy=(x_Rcparts, y_Rcparts),radius=errmaj_Rcparts)
               else:
                  a = Circle(xy=(x_Rcparts, y_Rcparts),radius=2.)
               a.set_linewidth(0.5)
               a.set_facecolor('none')
               plt.text(x_Rcparts,y_Rcparts+50,cat_R,horizontalalignment='center',verticalalignment='center',color='limegreen',fontsize=1.2*subsize,alpha=1.0)
               a.set_edgecolor('limegreen')
               ax1.add_artist(a)

      if 'ra_query' in locals() and 'dec_query' in locals():
         x_query, y_query = w.all_world2pix(ra_query, dec_query, 1)
         plt.text(x_query,y_query+50,"query",horizontalalignment='center',verticalalignment='center',color='darkorange',fontsize=1.2*subsize,alpha=1.0)
         a = Circle(xy=(x_query, y_query), radius=8)
         a.set_linewidth(0.5)
         a.set_facecolor('none')
         a.set_edgecolor('darkorange')
         ax1.add_artist(a)
         

      ax2 = plt.subplot(5, 2, 10)
      plt.title('Zoom on the variable source: (RA, DEC) = (%.3f'%ra+', %.3f'%dec+') deg',fontsize=subsize)
      plt.xlim(x-zoom,x+zoom)
      plt.ylim(y-zoom,y+zoom)
      plt.imshow(data, interpolation='nearest', cmap='gray', vmin=m-colscale*s, vmax=m+colscale*s, origin='lower')
      rxnew = [rx[j] for j in range(0,len(rx)) if rx[j]>=(x-zoom) and rx[j]<=(x+zoom) and ry[j]>=(y-zoom) and ry[j]<=(y+zoom)]
      rynew = [ry[j] for j in range(0,len(rx)) if rx[j]>=(x-zoom) and rx[j]<=(x+zoom) and ry[j]>=(y-zoom) and ry[j]<=(y+zoom)]
      for j in range(len(rxnew)):
         if not(np.round(x,ndig)==np.round(rxnew[j],ndig) and np.round(y,ndig)==np.round(rynew[j],ndig)):
            a = Circle(xy=(rxnew[j], rynew[j]), radius=fext*fwhmavg)
            a.set_linewidth(1.0)
            a.set_facecolor('none')
            a.set_edgecolor('blue')
            ax2.add_artist(a)
      plt.text(x,y+20,str(i+1),horizontalalignment='center',verticalalignment='center',color='r',fontsize=2*subsize)
      a = Circle(xy=(x, y), radius=fext*fwhmavg)
      a.set_linewidth(1.0)
      a.set_facecolor('none')
      a.set_edgecolor('red')
      ax2.add_artist(a)
      #if hdr['INSTRUME']=='WIFSIP V2.0':
      #   if np.size(obs_Rcparts)!=0:
      #      xnew_Rcparts = [x_Rcparts[j] for j in range(0,len(x_Rcparts)) if x_Rcparts[j]>=(x-zoom) and x_Rcparts[j]<=(x+zoom) and y_Rcparts[j]>=(y-zoom) and y_Rcparts[j]<=(y+zoom)]
      #      ynew_Rcparts = [y_Rcparts[j] for j in range(0,len(x_Rcparts)) if x_Rcparts[j]>=(x-zoom) and x_Rcparts[j]<=(x+zoom) and y_Rcparts[j]>=(y-zoom) and y_Rcparts[j]<=(y+zoom)]
      #      obsnew_Rcparts = [obs_Rcparts[j] for j in range(0,len(x_Rcparts)) if x_Rcparts[j]>=(x-zoom) and x_Rcparts[j]<=(x+zoom) and y_Rcparts[j]>=(y-zoom) and y_Rcparts[j]<=(y+zoom)]
      #      for j in range(0,np.size(xnew_Rcparts)):
      #         a = Circle(xy=(xnew_Rcparts[j], ynew_Rcparts[j]), radius=8)
      #         a.set_linewidth(1.0)
      #         a.set_facecolor('none')
      #         plt.text(xnew_Rcparts[j],ynew_Rcparts[j],obsnew_Rcparts[j],horizontalalignment='center',verticalalignment='center',color='limegreen',fontsize=1.5*subsize,alpha=1.0)
      #         a.set_edgecolor('limegreen')
      #         ax2.add_artist(a)
      if os.path.exists(path+'/2sxps.cat'):
         if np.size(x_Xcparts_Swift)!=1:
            xnew_Xcparts = [x_Xcparts_Swift[j] for j in range(0,len(x_Xcparts_Swift)) if x_Xcparts_Swift[j]>=(x-zoom) and x_Xcparts_Swift[j]<=(x+zoom) and y_Xcparts_Swift[j]>=(y-zoom) and y_Xcparts_Swift[j]<=(y+zoom)]
            ynew_Xcparts = [y_Xcparts_Swift[j] for j in range(0,len(x_Xcparts_Swift)) if x_Xcparts_Swift[j]>=(x-zoom) and x_Xcparts_Swift[j]<=(x+zoom) and y_Xcparts_Swift[j]>=(y-zoom) and y_Xcparts_Swift[j]<=(y+zoom)]
            errnew_Xcparts = [err_Xcparts_Swift[j] for j in range(0,len(x_Xcparts_Swift)) if x_Xcparts_Swift[j]>=(x-zoom) and x_Xcparts_Swift[j]<=(x+zoom) and y_Xcparts_Swift[j]>=(y-zoom) and y_Xcparts_Swift[j]<=(y+zoom)]
            for j in range(0,np.size(xnew_Xcparts)):
               a = Circle(xy=(xnew_Xcparts[j], ynew_Xcparts[j]), radius=errnew_Xcparts[j])
               a.set_linewidth(1.0)
               a.set_facecolor('none')
               plt.text(xnew_Xcparts[j],ynew_Xcparts[j],'swift',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.5*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax2.add_artist(a)
         else:
            if x_Xcparts_Swift>=(x-zoom) and x_Xcparts_Swift<=(x+zoom) and y_Xcparts_Swift>=(y-zoom) and y_Xcparts_Swift<=(y+zoom):
               xnew_Xcparts = x_Xcparts_Swift
               ynew_Xcparts = y_Xcparts_Swift
               errnew_Xcparts = err_Xcparts_Swift
               a = Circle(xy=(xnew_Xcparts, ynew_Xcparts), radius=errnew_Xcparts)
               a.set_linewidth(1.0)
               a.set_facecolor('none')
               plt.text(xnew_Xcparts,ynew_Xcparts,'swift',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.5*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax2.add_artist(a)
      if os.path.exists(path+'/csc2.cat'):
         if np.size(x_Xcparts_Chandra)!=1:
            xnew_Xcparts = [x_Xcparts_Chandra[j] for j in range(0,len(x_Xcparts_Chandra)) if x_Xcparts_Chandra[j]>=(x-zoom) and x_Xcparts_Chandra[j]<=(x+zoom) and y_Xcparts_Chandra[j]>=(y-zoom) and y_Xcparts_Chandra[j]<=(y+zoom)]
            ynew_Xcparts = [y_Xcparts_Chandra[j] for j in range(0,len(x_Xcparts_Chandra)) if x_Xcparts_Chandra[j]>=(x-zoom) and x_Xcparts_Chandra[j]<=(x+zoom) and y_Xcparts_Chandra[j]>=(y-zoom) and y_Xcparts_Chandra[j]<=(y+zoom)]
            ellanew_Xcparts = [ell_a_Chandra[j] for j in range(0,len(x_Xcparts_Chandra)) if x_Xcparts_Chandra[j]>=(x-zoom) and x_Xcparts_Chandra[j]<=(x+zoom) and y_Xcparts_Chandra[j]>=(y-zoom) and y_Xcparts_Chandra[j]<=(y+zoom)]
            ellbnew_Xcparts = [ell_b_Chandra[j] for j in range(0,len(x_Xcparts_Chandra)) if x_Xcparts_Chandra[j]>=(x-zoom) and x_Xcparts_Chandra[j]<=(x+zoom) and y_Xcparts_Chandra[j]>=(y-zoom) and y_Xcparts_Chandra[j]<=(y+zoom)]
            ellthetanew_Xcparts = [ell_theta_Chandra[j] for j in range(0,len(x_Xcparts_Chandra)) if x_Xcparts_Chandra[j]>=(x-zoom) and x_Xcparts_Chandra[j]<=(x+zoom) and y_Xcparts_Chandra[j]>=(y-zoom) and y_Xcparts_Chandra[j]<=(y+zoom)]
            for j in range(0,np.size(xnew_Xcparts)):
               a = Ellipse(xy=(xnew_Xcparts[j], ynew_Xcparts[j]), width=2*ellanew_Xcparts[j], height=2*ellbnew_Xcparts[j], angle=ellthetanew_Xcparts[j])
               a.set_linewidth(1.0)
               a.set_facecolor('none')
               plt.text(xnew_Xcparts[j],ynew_Xcparts[j],'chandra',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.5*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax2.add_artist(a)
         else:
            if x_Xcparts_Chandra>=(x-zoom) and x_Xcparts_Chandra<=(x+zoom) and y_Xcparts_Chandra>=(y-zoom) and y_Xcparts_Chandra<=(y+zoom):
               xnew_Xcparts = x_Xcparts_Chandra
               ynew_Xcparts = y_Xcparts_Chandra
               ellanew_Xcparts = ell_a_Chandra
               ellbnew_Xcparts = ell_b_Chandra
               ellthetanew_Xcparts = ell_theta_Chandra
               a = Ellipse(xy=(xnew_Xcparts, ynew_Xcparts), width=2*ellanew_Xcparts, height=2*ellbnew_Xcparts, angle=ellthetanew_Xcparts)
               a.set_linewidth(1.0)
               a.set_facecolor('none')
               plt.text(xnew_Xcparts,ynew_Xcparts,'chandra',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.5*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax2.add_artist(a)
      if os.path.exists(path+'/XMM4d13s.cat'):
         if np.size(x_Xcparts_XMM)!=1:
            xnew_Xcparts = [x_Xcparts_XMM[j] for j in range(0,len(x_Xcparts_XMM)) if x_Xcparts_XMM[j]>=(x-zoom) and x_Xcparts_XMM[j]<=(x+zoom) and y_Xcparts_XMM[j]>=(y-zoom) and y_Xcparts_XMM[j]<=(y+zoom)]
            ynew_Xcparts = [y_Xcparts_XMM[j] for j in range(0,len(x_Xcparts_XMM)) if x_Xcparts_XMM[j]>=(x-zoom) and x_Xcparts_XMM[j]<=(x+zoom) and y_Xcparts_XMM[j]>=(y-zoom) and y_Xcparts_XMM[j]<=(y+zoom)]
            errnew_Xcparts = [err_Xcparts_XMM[j] for j in range(0,len(x_Xcparts_XMM)) if x_Xcparts_XMM[j]>=(x-zoom) and x_Xcparts_XMM[j]<=(x+zoom) and y_Xcparts_XMM[j]>=(y-zoom) and y_Xcparts_XMM[j]<=(y+zoom)]
            for j in range(0,np.size(xnew_Xcparts)):
               a = Circle(xy=(xnew_Xcparts[j], ynew_Xcparts[j]), radius=errnew_Xcparts[j])
               a.set_linewidth(1.0)
               a.set_facecolor('none')
               plt.text(xnew_Xcparts[j],ynew_Xcparts[j],'XMM',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.5*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax2.add_artist(a)
         else:
            if x_Xcparts_XMM>=(x-zoom) and x_Xcparts_XMM<=(x+zoom) and y_Xcparts_XMM>=(y-zoom) and y_Xcparts_XMM<=(y+zoom):
               xnew_Xcparts = x_Xcparts_XMM
               ynew_Xcparts = y_Xcparts_XMM
               errnew_Xcparts = err_Xcparts_XMM
               a = Circle(xy=(xnew_Xcparts, ynew_Xcparts), radius=errnew_Xcparts)
               a.set_linewidth(1.0)
               a.set_facecolor('none')
               plt.text(xnew_Xcparts,ynew_Xcparts,'XMM',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.5*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax2.add_artist(a)
      if os.path.exists(path+'/eRASS1d1.cat'):
         if np.size(x_Xcparts_eRASS)!=1:
            xnew_Xcparts = [x_Xcparts_eRASS[j] for j in range(0,len(x_Xcparts_eRASS)) if x_Xcparts_eRASS[j]>=(x-zoom) and x_Xcparts_eRASS[j]<=(x+zoom) and y_Xcparts_eRASS[j]>=(y-zoom) and y_Xcparts_eRASS[j]<=(y+zoom)]
            ynew_Xcparts = [y_Xcparts_eRASS[j] for j in range(0,len(x_Xcparts_eRASS)) if x_Xcparts_eRASS[j]>=(x-zoom) and x_Xcparts_eRASS[j]<=(x+zoom) and y_Xcparts_eRASS[j]>=(y-zoom) and y_Xcparts_eRASS[j]<=(y+zoom)]
            errnew_Xcparts = [err_Xcparts_eRASS[j] for j in range(0,len(x_Xcparts_eRASS)) if x_Xcparts_eRASS[j]>=(x-zoom) and x_Xcparts_eRASS[j]<=(x+zoom) and y_Xcparts_eRASS[j]>=(y-zoom) and y_Xcparts_eRASS[j]<=(y+zoom)]
            for j in range(0,np.size(xnew_Xcparts)):
               a = Circle(xy=(xnew_Xcparts[j], ynew_Xcparts[j]), radius=errnew_Xcparts[j])
               a.set_linewidth(1.0)
               a.set_facecolor('none')
               plt.text(xnew_Xcparts[j],ynew_Xcparts[j],'eRASS',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.5*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax2.add_artist(a)
         else:
            if x_Xcparts_eRASS>=(x-zoom) and x_Xcparts_eRASS<=(x+zoom) and y_Xcparts_eRASS>=(y-zoom) and y_Xcparts_eRASS<=(y+zoom):
               xnew_Xcparts = x_Xcparts_eRASS
               ynew_Xcparts = y_Xcparts_eRASS
               errnew_Xcparts = err_Xcparts_eRASS
               a = Circle(xy=(xnew_Xcparts, ynew_Xcparts), radius=errnew_Xcparts)
               a.set_linewidth(1.0)
               a.set_facecolor('none')
               plt.text(xnew_Xcparts,ynew_Xcparts,'eRASS',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.5*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax2.add_artist(a)
      if np.size(x_Rcparts) > 0:
         if np.size(x_Rcparts) > 1:
            xnew_Rcparts = [x_Rcparts[j] for j in range(0,len(x_Rcparts)) if x_Rcparts[j]>=(x-zoom) and x_Rcparts[j]<=(x+zoom) and y_Rcparts[j]>=(y-zoom) and y_Rcparts[j]<=(y+zoom)]
            ynew_Rcparts = [y_Rcparts[j] for j in range(0,len(x_Rcparts)) if x_Rcparts[j]>=(x-zoom) and x_Rcparts[j]<=(x+zoom) and y_Rcparts[j]>=(y-zoom) and y_Rcparts[j]<=(y+zoom)]
            errmajnew_Rcparts = [errmaj_Rcparts[j] for j in range(0,len(x_Rcparts)) if x_Rcparts[j]>=(x-zoom) and x_Rcparts[j]<=(x+zoom) and y_Rcparts[j]>=(y-zoom) and y_Rcparts[j]<=(y+zoom)]
            errminnew_Rcparts = [errmin_Rcparts[j] for j in range(0,len(x_Rcparts)) if x_Rcparts[j]>=(x-zoom) and x_Rcparts[j]<=(x+zoom) and y_Rcparts[j]>=(y-zoom) and y_Rcparts[j]<=(y+zoom)]
            erranglenew_Rcparts = [errangle_Rcparts[j] for j in range(0,len(x_Rcparts)) if x_Rcparts[j]>=(x-zoom) and x_Rcparts[j]<=(x+zoom) and y_Rcparts[j]>=(y-zoom) and y_Rcparts[j]<=(y+zoom)]
            catnew_R = [cat_R[j] for j in range(0,len(x_Rcparts)) if x_Rcparts[j]>=(x-zoom) and x_Rcparts[j]<=(x+zoom) and y_Rcparts[j]>=(y-zoom) and y_Rcparts[j]<=(y+zoom)]
            for j in range(0,np.size(xnew_Rcparts)):
               if not np.isnan(errmajnew_Rcparts[j]) and not np.isnan(errminnew_Rcparts[j]) and not np.isnan(erranglenew_Rcparts[j]):
                  a = Ellipse(xy=(xnew_Rcparts[j], ynew_Rcparts[j]), width=2*errmajnew_Rcparts[j], height=2*errminnew_Rcparts[j], angle=erranglenew_Rcparts[j])
               elif not np.isnan(errmajnew_Rcparts[j]) and not np.isnan(errminnew_Rcparts[j]) and np.isnan(erranglenew_Rcparts[j]):
                  a = Circle(xy=(xnew_Rcparts[j], ynew_Rcparts[j]),radius=errmajnew_Rcparts[j])
               else:
                  a = Circle(xy=(xnew_Rcparts[j], ynew_Rcparts[j]),radius=2.)
               a.set_linewidth(1.0)
               a.set_facecolor('none')
               plt.text(xnew_Rcparts[j],ynew_Rcparts[j],catnew_R[j],horizontalalignment='center',verticalalignment='center',color='limegreen',fontsize=1.5*subsize,alpha=1.0)
               a.set_edgecolor('limegreen')
               ax2.add_artist(a)
         else:
            if x_Rcparts>=(x-zoom) and x_Rcparts<=(x+zoom) and y_Rcparts>=(y-zoom) and y_Rcparts<=(y+zoom):
               xnew_Rcparts = x_Rcparts
               ynew_Rcparts = y_Rcparts
               errmajnew_Rcparts = errmaj_Rcparts
               errminnew_Rcparts = errmin_Rcparts
               erranglenew_Rcparts = errangle_Rcparts
               catnew_R = cat_R
               if not np.isnan(errmajnew_Rcparts) and not np.isnan(errminnew_Rcparts) and not np.isnan(erranglenew_Rcparts):
                  a = Ellipse(xy=(xnew_Rcparts, ynew_Rcparts), width=2*errmajnew_Rcparts, height=2*errminnew_Rcparts, angle=erranglenew_Rcparts)
               elif not np.isnan(errmajnew_Rcparts) and not np.isnan(errminnew_Rcparts) and np.isnan(erranglenew_Rcparts):
                  a = Circle(xy=(xnew_Rcparts, ynew_Rcparts),radius=errmajnew_Rcparts)
               else:
                  a = Circle(xy=(xnew_Rcparts, ynew_Rcparts),radius=2.)
               a.set_linewidth(1.0)
               a.set_facecolor('none')
               plt.text(xnew_Rcparts,ynew_Rcparts,catnew_R,horizontalalignment='center',verticalalignment='center',color='limegreen',fontsize=1.5*subsize,alpha=1.0)
               a.set_edgecolor('limegreen')
               ax2.add_artist(a)

      if 'ra_query' in locals() and 'dec_query' in locals():
         plt.text(x_query,y_query,"query",horizontalalignment='center',verticalalignment='center',color='darkorange',fontsize=1.5*subsize,alpha=0.8)
         a = Circle(xy=(x_query, y_query), radius=8)
         a.set_linewidth(0.5)
         a.set_facecolor('none')
         a.set_edgecolor('darkorange')
         ax2.add_artist(a)

      if j==0:
         if 'kic_ID' in locals():
            plt.text(-1.85,6.0,"Identified from Kepler as "+kic_ID+r" with $P_{\mathrm{Kepler}}$="+str(kic_period)+" d",transform=ax2.transAxes,fontsize=tsize)
         elif 'ATO_ID' in locals():
            plt.text(-1.85,6.0,"Identified from ATLAS as "+ATO_ID+r" with $P_{\mathrm{ATLAS}}$="+str(ATO_period)+" d",transform=ax2.transAxes,fontsize=tsize)
         elif 'SIMBAD_ID' in locals():
            plt.text(-1.7,6.0,"Identified from SIMBAD as "+SIMBAD_ID,transform=ax2.transAxes,fontsize=tsize)
         plt.savefig(path+'/finalperiodplots'+str(nlist)+'_'+per_search+'/finalperplot_'+str(i+1)+'_'+per_search+'_'+analysis+'_nper'+str(int(n_periods))+'_ZTF.png',bbox_inches='tight')
      else:
         if 'kic_ID' in locals():
            plt.text(-1.85,6.0,"Identified from Kepler as "+kic_ID+r" with $P_{\mathrm{Kepler}}$="+str(kic_period)+" d",transform=ax2.transAxes,fontsize=tsize)
            plt.text(-0.87,5.85,"SIMBAD query",url="https://simbad.cds.unistra.fr/simbad/sim-coo?Coord="+str(ra)+"+"+str(dec)+"&CooFrame=FK5&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=none&Radius=3&Radius.unit=arcsec&submit=submit+query&CoordList=",transform=ax2.transAxes,fontsize=tsize,color='blue')
         elif 'ATO_ID' in locals():
            plt.text(-1.85,6.0,"Identified from ATLAS as "+ATO_ID+r" with $P_{\mathrm{ATLAS}}$="+str(ATO_period)+" d",transform=ax2.transAxes,fontsize=tsize)
            plt.text(-0.87,5.85,"SIMBAD query",url="https://simbad.cds.unistra.fr/simbad/sim-coo?Coord="+str(ra)+"+"+str(dec)+"&CooFrame=FK5&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=none&Radius=3&Radius.unit=arcsec&submit=submit+query&CoordList=",transform=ax2.transAxes,fontsize=tsize,color='blue')
         elif 'SIMBAD_ID' in locals():
            plt.text(-1.7,6.0,"Identified from SIMBAD as "+SIMBAD_ID,transform=ax2.transAxes,fontsize=tsize)
            plt.text(-0.87,5.85,"SIMBAD query",url="https://simbad.cds.unistra.fr/simbad/sim-coo?Coord="+str(ra)+"+"+str(dec)+"&CooFrame=FK5&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=none&Radius=3&Radius.unit=arcsec&submit=submit+query&CoordList=",transform=ax2.transAxes,fontsize=tsize,color='blue')
         else:
            plt.text(-0.87,6.0,"SIMBAD query",url="https://simbad.cds.unistra.fr/simbad/sim-coo?Coord="+str(ra)+"+"+str(dec)+"&CooFrame=FK5&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=none&Radius=3&Radius.unit=arcsec&submit=submit+query&CoordList=",transform=ax2.transAxes,fontsize=tsize,color='blue')
         plt.savefig(path+'/finalperiodplots'+str(nlist)+'_'+per_search+'/finalperplot_'+str(i+1)+'_'+per_search+'_'+analysis+'_nper'+str(int(n_periods))+'_ZTF.pdf',bbox_inches='tight')
      plt.close(figure)

   #Perform Lomb-Scargle/PDM periodograms on ZTF data to get better estimates of the period (better sampling), and fold light curves at the period, half the period and twice the period and t_0 in g and r bands, with the period from the eventual corresponding ATLAS variable source
   if g_check == 1 or r_check==1:
      
      figure=plt.figure(figsize=(figxsize,figysize), dpi=150)
      if np.size(sources_3FGL)!=0:
         if 'flux_Xcpart' in locals():
            figure.suptitle(r"Fermi source: "+FoVname+", Periodic variable #%d"%(i+1)+", $\sigma_{\mathrm{ZTF}}$: %.2f"%ZTF_std+r" mag, $\sigma_{\mathrm{"+telescope.split("/")[0]+"}}$: %.2f"%sigma+" mag\n\nRange: [%.2f"%P_low+",%.1f"%min(P_up,(ZTF_mjd_end-ZTF_mjd_start))+"] d, Steps: %d"%n_periods+", Period folding with "+r"$T_{0}=$%.4f"%t_0+" MJD\n\n%.1f"%(sep_4FGL/ell_a_68)+r"$\sigma$ far from the 4FGL source and %.1f"%(sep_3FGL/ell_a_68_3FGL)+r"$\sigma$ far from the 3FGL source, X-ray flux: "+str(np.round(float(str(flux_Xcpart).split('e')[0]),1))+"$\cdot 10^{"+str(flux_Xcpart).split('e')[1]+"}$ erg s$^{-1}$ cm$^{-2}$", fontsize=tsize)
         else:
            if 'UL_X' in locals():
               figure.suptitle(r"Fermi source: "+FoVname+", Periodic variable #%d"%(i+1)+", $\sigma_{\mathrm{ZTF}}$: %.2f"%ZTF_std+r" mag, $\sigma_{\mathrm{"+telescope.split("/")[0]+"}}$: %.2f"%sigma+" mag\n\nRange: [%.2f"%P_low+",%.1f"%min(P_up,(ZTF_mjd_end-ZTF_mjd_start))+"] d, Steps: %d"%n_periods+", Period folding with "+r"$T_{0}=$%.4f"%t_0+" MJD\n\n%.1f"%(sep_4FGL/ell_a_68)+r"$\sigma$ far from the 4FGL source and %.1f"%(sep_3FGL/ell_a_68_3FGL)+r"$\sigma$ far from the 3FGL source, X-ray flux 3$\sigma$ UL: "+str(np.round(float(str(UL_X).split('e')[0]),1))+"$\cdot 10^{"+str(UL_X).split('e')[1]+"}$ erg s$^{-1}$ cm$^{-2}$", fontsize=tsize)
            else:
               figure.suptitle(r"Fermi source: "+FoVname+", Periodic variable #%d"%(i+1)+", $\sigma_{\mathrm{ZTF}}$: %.2f"%ZTF_std+r" mag, $\sigma_{\mathrm{"+telescope.split("/")[0]+"}}$: %.2f"%sigma+" mag\n\nRange: [%.2f"%P_low+",%.1f"%min(P_up,(ZTF_mjd_end-ZTF_mjd_start))+"] d, Steps: %d"%n_periods+", Period folding with "+r"$T_{0}=$%.4f"%t_0+" MJD\n\n%.1f"%(sep_4FGL/ell_a_68)+r"$\sigma$ far from the 4FGL source and %.1f"%(sep_3FGL/ell_a_68_3FGL)+r"$\sigma$ far from the 3FGL source", fontsize=tsize)
      else:
         if 'flux_Xcpart' in locals():
            figure.suptitle(r"Fermi source: "+FoVname+", Periodic variable #%d"%(i+1)+", $\sigma_{\mathrm{ZTF}}$: %.2f"%ZTF_std+r" mag, $\sigma_{\mathrm{"+telescope.split("/")[0]+"}}$: %.2f"%sigma+" mag\n\nRange: [%.2f"%P_low+",%.1f"%min(P_up,(ZTF_mjd_end-ZTF_mjd_start))+"] d, Steps: %d"%n_periods+", Period folding with "+r"$T_{0}=$%.4f"%t_0+" MJD\n\n%.1f"%(sep_4FGL/ell_a_68)+r"$\sigma$ far from the 4FGL source, X-ray flux: "+str(np.round(float(str(flux_Xcpart).split('e')[0]),1))+"$\cdot 10^{"+str(flux_Xcpart).split('e')[1]+"}$ erg s$^{-1}$ cm$^{-2}$", fontsize=tsize)
         else:
            if 'UL_X' in locals():
               figure.suptitle(r"Fermi source: "+FoVname+", Periodic variable #%d"%(i+1)+", $\sigma_{\mathrm{ZTF}}$: %.2f"%ZTF_std+r" mag, $\sigma_{\mathrm{"+telescope.split("/")[0]+"}}$: %.2f"%sigma+" mag\n\nRange: [%.2f"%P_low+",%.1f"%min(P_up,(ZTF_mjd_end-ZTF_mjd_start))+"] d, Steps: %d"%n_periods+", Period folding with "+r"$T_{0}=$%.4f"%t_0+" MJD\n\n%.1f"%(sep_4FGL/ell_a_68)+r"$\sigma$ far from the 4FGL source, X-ray flux 3$\sigma$ UL: "+str(np.round(float(str(UL_X).split('e')[0]),1))+"$\cdot 10^{"+str(UL_X).split('e')[1]+"}$ erg s$^{-1}$ cm$^{-2}$", fontsize=tsize)
            else:
               figure.suptitle(r"Fermi source: "+FoVname+", Periodic variable #%d"%(i+1)+", $\sigma_{\mathrm{ZTF}}$: %.2f"%ZTF_std+r" mag, $\sigma_{\mathrm{"+telescope.split("/")[0]+"}}$: %.2f"%sigma+" mag\n\nRange: [%.2f"%P_low+",%.1f"%min(P_up,(ZTF_mjd_end-ZTF_mjd_start))+"] d, Steps: %d"%n_periods+", Period folding with "+r"$T_{0}=$%.4f"%t_0+" MJD\n\n%.1f"%(sep_4FGL/ell_a_68)+r"$\sigma$ far from the 4FGL source", fontsize=tsize)

      if g_check == 1:
         if per_search == 'LS':
            periods, periodogram = LombScargle_search(ztf_mjd_g, ztf_mag_g, ztf_magerr_g, P_low, P_up, n_periods, 1)
            ztfg_best_period = periods[~np.isnan(periodogram)][np.argmax(periodogram[~np.isnan(periodogram)])]
         else:
            if per_search == 'PDM':
               periods, periodogram = PDM_search(ztf_mjd_g, ztf_mag_g, P_low, P_up, 0., n_periods, n_bin, min_bin)
            else:
               if g_check == 1 and r_check == 1:
                  periods, periodogram = PDM_bisearch(ztf_mjd_g, ztf_mag_g, ztf_mjd_r, ztf_mag_r, P_low, P_up, 0., n_periods, n_bin, min_bin)
               else:
                  periods, periodogram = PDM_search(ztf_mjd_g, ztf_mag_g, P_low, P_up, 0., n_periods, n_bin, min_bin)
            if np.size(periodogram[~np.isnan(periodogram)])!=0:
               ztfg_best_period = periods[~np.isnan(periodogram)][np.argmin(periodogram[~np.isnan(periodogram)])]
         
         ztf_phase_g = ((ztf_mjd_g-t_0)/ztfg_best_period) % 1
         plt.subplot(5, 2, 1)
         plt.plot(ztf_phase_g, ztf_mag_g, 'bo', linestyle='None')
         plt.plot(ztf_phase_g+1, ztf_mag_g, 'ro', linestyle='None')
         plt.errorbar(ztf_phase_g, ztf_mag_g, ztf_magerr_g, linestyle='None')
         plt.errorbar(ztf_phase_g+1, ztf_mag_g, ztf_magerr_g, linestyle='None')
         plt.xticks(fontsize=axticksize)
         plt.yticks(fontsize=axticksize)
         plt.gca().invert_yaxis()
         plt.title(r"Folded ZTF lightcurve in the g band: $P_{\mathrm{ZTF,g}}$="+"{0} d".format(ztfg_best_period),fontsize=subsize)
         plt.xlabel(r"Phase ($\phi$)",fontsize=axlabelsize)
         plt.ylabel(r"Apparent magnitude",fontsize=axlabelsize)

         plt.subplot(5, 2, 2)
         plt.plot(periods,periodogram)
         if not np.isnan(max(periodogram)) and not np.isinf(max(periodogram)) and not np.isnan(min(periodogram)) and not np.isinf(min(periodogram)) and not np.isnan(ztfg_best_period) and not np.isinf(ztfg_best_period):
            plt.axvline(x=ztfg_best_period, color='r', ls='--')
         plt.xscale("log")
         plt.xticks(fontsize=axticksize)
         plt.yticks(fontsize=axticksize)
         if per_search == 'LS':
            plt.title("Lomb-Scargle ZTF data periodogram in the g band",fontsize=subsize)
            plt.ylabel(r"Power",fontsize=axlabelsize)
            if period_FAPl == 'yes':
               tic_1 = time.perf_counter()
               FAP_periods, FAP_levels = LombScargle_FAP_level(ztf_mjd_g, ztf_mag_g, ztf_magerr_g, P_low, P_up, n_periods, 1, FAPl_factor, conf_level)
               toc_1 = time.perf_counter()
               print(f"FAP levels computed in {toc_1 - tic_1:0.4f} seconds")
               plt.plot(FAP_periods,FAP_levels,color='darkorange',ls='--')
               plt.text(plt.xlim()[0]*1.4, np.mean(FAP_levels)*1.2, r'0.1% FAP', color='darkorange', ha='center',fontsize=axticksize)
         elif per_search == 'PDM':
            plt.title("PDM ZTF data periodogram in the g band",fontsize=subsize)
            plt.ylabel(r"Normalized variance",fontsize=axlabelsize)
            if period_FAPl == 'yes':
               tic_1 = time.perf_counter()
               FAP_periods, FAP_levels = PDM_FAP_level(ztf_mjd_g, ztf_mag_g, P_low, P_up, 0., n_periods, n_bin, min_bin, FAPl_factor, conf_level)
               toc_1 = time.perf_counter()
               print(f"FAP levels computed in {toc_1 - tic_1:0.4f} seconds")
               plt.plot(FAP_periods,FAP_levels,color='darkorange',ls='--')
               plt.text(plt.xlim()[0]*1.4, np.mean(FAP_levels)*1.02, r'0.1% FAP', color='darkorange', ha='center',fontsize=axticksize)
         else:
            plt.title("PDM multi-band ZTF data periodogram",fontsize=subsize)
            plt.ylabel(r"Sum of normalized variances",fontsize=axlabelsize)
         plt.xlabel(r"Trial period (days)",fontsize=axlabelsize)

      if r_check == 1:
         if per_search == 'LS':
            periods, periodogram = LombScargle_search(ztf_mjd_r, ztf_mag_r, ztf_magerr_r, P_low, P_up, n_periods, 1)
            ztfr_best_period = periods[~np.isnan(periodogram)][np.argmax(periodogram[~np.isnan(periodogram)])]
         else:
            if per_search == 'PDM':
               periods, periodogram = PDM_search(ztf_mjd_r, ztf_mag_r, P_low, P_up, 0., n_periods, n_bin, min_bin)
            else:
               if g_check == 1 and r_check == 1:
                  periods, periodogram = PDM_bisearch(ztf_mjd_g, ztf_mag_g, ztf_mjd_r, ztf_mag_r, P_low, P_up, 0., n_periods, n_bin, min_bin)
               else:
                  periods, periodogram = PDM_search(ztf_mjd_r, ztf_mag_r, P_low, P_up, 0., n_periods, n_bin, min_bin)
            if np.size(periodogram[~np.isnan(periodogram)])!=0:
               ztfr_best_period = periods[~np.isnan(periodogram)][np.argmin(periodogram[~np.isnan(periodogram)])]
         
         ztf_phase_r = ((ztf_mjd_r-t_0)/ztfr_best_period) % 1
         plt.subplot(5, 2, 3)
         plt.plot(ztf_phase_r, ztf_mag_r, 'bo', linestyle='None')
         plt.plot(ztf_phase_r+1, ztf_mag_r, 'ro', linestyle='None')
         plt.errorbar(ztf_phase_r, ztf_mag_r, ztf_magerr_r, linestyle='None')
         plt.errorbar(ztf_phase_r+1, ztf_mag_r, ztf_magerr_r, linestyle='None')
         plt.xticks(fontsize=axticksize)
         plt.yticks(fontsize=axticksize)
         plt.gca().invert_yaxis()
         plt.title(r"Folded ZTF lightcurve in the r band: $P_{\mathrm{ZTF,r}}$="+"{0} d".format(ztfr_best_period),fontsize=subsize)
         plt.xlabel(r"Phase ($\phi$)",fontsize=axlabelsize)
         plt.ylabel(r"Apparent magnitude",fontsize=axlabelsize)

         plt.subplot(5, 2, 4)
         plt.plot(periods,periodogram)
         if not np.isnan(max(periodogram)) and not np.isinf(max(periodogram)) and not np.isnan(min(periodogram)) and not np.isinf(min(periodogram)) and not np.isnan(ztfr_best_period) and not np.isinf(ztfr_best_period):
            plt.axvline(x=ztfr_best_period, color='r', ls='--')
         plt.xscale("log")
         plt.xticks(fontsize=axticksize)
         plt.yticks(fontsize=axticksize)
         if per_search == 'LS':
            plt.title("Lomb-Scargle ZTF data periodogram in the r band",fontsize=subsize)
            plt.ylabel(r"Power",fontsize=axlabelsize)
            if period_FAPl == 'yes':
               tic_1 = time.perf_counter()
               FAP_periods, FAP_levels = LombScargle_FAP_level(ztf_mjd_r, ztf_mag_r, ztf_magerr_r, P_low, P_up, n_periods, 1, FAPl_factor, conf_level)
               toc_1 = time.perf_counter()
               print(f"FAP levels computed in {toc_1 - tic_1:0.4f} seconds")
               plt.plot(FAP_periods,FAP_levels,color='darkorange',ls='--')
               plt.text(plt.xlim()[0]*1.4, np.mean(FAP_levels)*1.2, r'0.1% FAP', color='darkorange', ha='center',fontsize=axticksize)
         elif per_search == 'PDM':
            plt.title("PDM ZTF data periodogram in the r band",fontsize=subsize)
            plt.ylabel(r"Normalized variance",fontsize=axlabelsize)
            if period_FAPl == 'yes':
               tic_1 = time.perf_counter()
               FAP_periods, FAP_levels = PDM_FAP_level(ztf_mjd_r, ztf_mag_r, P_low, P_up, 0., n_periods, n_bin, min_bin, FAPl_factor, conf_level)
               toc_1 = time.perf_counter()
               print(f"FAP levels computed in {toc_1 - tic_1:0.4f} seconds")
               plt.plot(FAP_periods,FAP_levels,color='darkorange',ls='--')
               plt.text(plt.xlim()[0]*1.4, np.mean(FAP_levels)*1.02, r'0.1% FAP', color='darkorange', ha='center',fontsize=axticksize)
         else:
            plt.title("PDM multi-band ZTF data periodogram",fontsize=subsize)
            plt.ylabel(r"Sum of normalized variances",fontsize=axlabelsize)
         plt.xlabel(r"Trial period (days)",fontsize=axlabelsize)

      if g_check == 1:
         ztf_phase_g = ((ztf_mjd_g-t_0)/(0.5*ztfg_best_period)) % 1
         plt.subplot(5, 2, 5)
         plt.plot(ztf_phase_g, ztf_mag_g, 'bo', linestyle='None')
         plt.plot(ztf_phase_g+1, ztf_mag_g, 'ro', linestyle='None')
         plt.errorbar(ztf_phase_g, ztf_mag_g, ztf_magerr_g, linestyle='None')
         plt.errorbar(ztf_phase_g+1, ztf_mag_g, ztf_magerr_g, linestyle='None')
         plt.xticks(fontsize=axticksize)
         plt.yticks(fontsize=axticksize)
         plt.gca().invert_yaxis()
         plt.title(r"ZTF lightcurve in the g band folded at $P_{\mathrm{ZTF,g}}/2$",fontsize=subsize)
         plt.xlabel(r"Phase ($\phi$)",fontsize=axlabelsize)
         plt.ylabel(r"Apparent magnitude",fontsize=axlabelsize)

         ztf_phase_g = ((ztf_mjd_g-t_0)/(2*ztfg_best_period)) % 1
         plt.subplot(5, 2, 6)
         plt.plot(ztf_phase_g, ztf_mag_g, 'bo', linestyle='None')
         plt.plot(ztf_phase_g+1, ztf_mag_g, 'ro', linestyle='None')
         plt.errorbar(ztf_phase_g, ztf_mag_g, ztf_magerr_g, linestyle='None')
         plt.errorbar(ztf_phase_g+1, ztf_mag_g, ztf_magerr_g, linestyle='None')
         plt.xticks(fontsize=axticksize)
         plt.yticks(fontsize=axticksize)
         plt.gca().invert_yaxis()
         plt.title(r"ZTF lightcurve in the g band folded at $2P_{\mathrm{ZTF,g}}$",fontsize=subsize)
         plt.xlabel(r"Phase ($\phi$)",fontsize=axlabelsize)
         plt.ylabel(r"Apparent magnitude",fontsize=axlabelsize)

      if r_check == 1:
         ztf_phase_r = ((ztf_mjd_r-t_0)/(0.5*ztfr_best_period)) % 1
         plt.subplot(5, 2, 7)
         plt.plot(ztf_phase_r, ztf_mag_r, 'bo', linestyle='None')
         plt.plot(ztf_phase_r+1, ztf_mag_r, 'ro', linestyle='None')
         plt.errorbar(ztf_phase_r, ztf_mag_r, ztf_magerr_r, linestyle='None')
         plt.errorbar(ztf_phase_r+1, ztf_mag_r, ztf_magerr_r, linestyle='None')
         plt.xticks(fontsize=axticksize)
         plt.yticks(fontsize=axticksize)
         plt.gca().invert_yaxis()
         plt.title(r"ZTF lightcurve in the r band folded at $P_{\mathrm{ZTF,r}}/2$",fontsize=subsize)
         plt.xlabel(r"Phase ($\phi$)",fontsize=axlabelsize)
         plt.ylabel(r"Apparent magnitude",fontsize=axlabelsize)

         ztf_phase_r = ((ztf_mjd_r-t_0)/(2*ztfr_best_period)) % 1
         plt.subplot(5, 2, 8)
         plt.plot(ztf_phase_r, ztf_mag_r, 'bo', linestyle='None')
         plt.plot(ztf_phase_r+1, ztf_mag_r, 'ro', linestyle='None')
         plt.errorbar(ztf_phase_r, ztf_mag_r, ztf_magerr_r, linestyle='None')
         plt.errorbar(ztf_phase_r+1, ztf_mag_r, ztf_magerr_r, linestyle='None')
         plt.xticks(fontsize=axticksize)
         plt.yticks(fontsize=axticksize)
         plt.gca().invert_yaxis()
         plt.title(r"ZTF lightcurve in the r band folded at $2P_{\mathrm{ZTF,r}}$",fontsize=subsize)
         plt.xlabel(r"Phase ($\phi$)",fontsize=axlabelsize)
         plt.ylabel(r"Apparent magnitude",fontsize=axlabelsize)
   
      ax1 = plt.subplot(5, 2, 9)
      plt.title(r'Median r FoV bkg subtracted in the '+str(2*bound)+'$\sigma$ square around the 4FGL coordinates',fontsize=subsize)
      plt.xlim(max(0,xtar-bound*max(ell_a_95,ell_b_95)), min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)))
      plt.ylim(max(0,ytar-bound*max(ell_a_95,ell_b_95)), min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)))
      plt.imshow(data, interpolation='nearest', cmap='gray', vmin=m-colscale*s, vmax=m+colscale*s, origin='lower')
      plt.text(xtar,ytar,FoVname,horizontalalignment='center',verticalalignment='center',color='cyan',fontsize=subsize,alpha=0.8)
      target = Ellipse(xy=(xtar, ytar), width=2*ell_a_95, height=2*ell_b_95, angle=ell_theta_95)
      target.set_facecolor('none')
      target.set_edgecolor('cyan')
      ax1.add_artist(target)
      if np.size(sources_3FGL)!=0:
         for j in range(0,np.shape(sources_3FGL)[0]):
            if x_3FGL[j]>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_3FGL[j]<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_3FGL[j]>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_3FGL[j]<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
               plt.text(x_3FGL[j],y_3FGL[j],sources_3FGL[j].replace("L ","L"),horizontalalignment='center',verticalalignment='center',color='yellow',fontsize=subsize,alpha=0.8)
               a = Ellipse(xy=(x_3FGL[j], y_3FGL[j]), width=2*3600./float(platescale)*coord_3FGL[j,3], height=2*3600./float(platescale)*coord_3FGL[j,4], angle=coord_3FGL[j,5])
               a.set_facecolor('none')
               a.set_edgecolor('yellow')
               ax1.add_artist(a)
      #if np.size(sources_4FGLe)!=0:
      #   for j in range(0,np.shape(sources_4FGLe)[0]):
      #      if x_4FGLe[j]>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_4FGLe[j]<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_4FGLe[j]>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_4FGLe[j]<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
      #         plt.text(x_4FGLe[j],y_4FGLe[j],sources_4FGLe[j].replace("L ","L"),horizontalalignment='center',verticalalignment='center',color='cyan',fontsize=subsize,alpha=0.8)
      #         a = Ellipse(xy=(x_4FGLe[j], y_4FGLe[j]), width=2*3600./float(platescale)*np.asarray(coord_4FGLe)[j,2], height=2*3600./float(platescale)*np.asarray(coord_4FGLe)[j,3], angle=np.asarray(coord_4FGLe)[j,4])
      #         a.set_facecolor('none')
      #         a.set_edgecolor('cyan')
      #         ax1.add_artist(a)
         
      x, y = w.all_world2pix(ra, dec, 1)
      plt.text(x,y+50,str(i+1),horizontalalignment='center',verticalalignment='center',color='r',fontsize=subsize)
      a = Circle(xy=(x, y), radius=fext*fwhmavg)
      a.set_linewidth(0.5)
      a.set_facecolor('none')
      a.set_edgecolor('red')
      ax1.add_artist(a)
      #if hdr['INSTRUME']=='WIFSIP V2.0':
      #   if np.size(obs_Rcparts)!=0:
      #      for j in range(0,np.size(obs_Rcparts)):
      #         if x_Rcparts[j]>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Rcparts[j]<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Rcparts[j]>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Rcparts[j]<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
      #            a = Circle(xy=(x_Rcparts[j], y_Rcparts[j]), radius=8)
      #            a.set_linewidth(0.5)
      #            a.set_facecolor('none')
      #            plt.text(x_Rcparts[j],y_Rcparts[j]+50,obs_Rcparts[j],horizontalalignment='center',verticalalignment='center',color='limegreen',fontsize=1.2*subsize,alpha=0.8)
      #            a.set_edgecolor('limegreen')
      #            ax1.add_artist(a)
      if os.path.exists(path+'/2sxps.cat'):
         if np.size(x_Xcparts_Swift)!=1:
            for j in range(0,np.size(x_Xcparts_Swift)):
               if x_Xcparts_Swift[j]>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Xcparts_Swift[j]<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Xcparts_Swift[j]>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Xcparts_Swift[j]<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
                  a = Circle(xy=(x_Xcparts_Swift[j], y_Xcparts_Swift[j]), radius=err_Xcparts_Swift[j])
                  a.set_linewidth(0.5)
                  a.set_facecolor('none')
                  plt.text(x_Xcparts_Swift[j],y_Xcparts_Swift[j]+50,'swift',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.2*subsize,alpha=1.0)
                  a.set_edgecolor('fuchsia')
                  ax1.add_artist(a)
         else:
            if x_Xcparts_Swift>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Xcparts_Swift<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Xcparts_Swift>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Xcparts_Swift<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
               a = Circle(xy=(x_Xcparts_Swift, y_Xcparts_Swift), radius=err_Xcparts_Swift)
               a.set_linewidth(0.5)
               a.set_facecolor('none')
               plt.text(x_Xcparts_Swift,y_Xcparts_Swift+50,'swift',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.2*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax1.add_artist(a)
      if os.path.exists(path+'/csc2.cat'):
         if np.size(x_Xcparts_Chandra)!=1:
            for j in range(0,np.size(x_Xcparts_Chandra)):
               if x_Xcparts_Chandra[j]>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Xcparts_Chandra[j]<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Xcparts_Chandra[j]>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Xcparts_Chandra[j]<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
                  a = Ellipse(xy=(x_Xcparts_Chandra[j], y_Xcparts_Chandra[j]), width=2*ell_a_Chandra[j], height=2*ell_b_Chandra[j], angle=ell_theta_Chandra[j])
                  a.set_linewidth(0.5)
                  a.set_facecolor('none')
                  plt.text(x_Xcparts_Chandra[j],y_Xcparts_Chandra[j]+50,'chandra',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.2*subsize,alpha=1.0)
                  a.set_edgecolor('fuchsia')
                  ax1.add_artist(a)
         else:
            if x_Xcparts_Chandra>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Xcparts_Chandra<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Xcparts_Chandra>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Xcparts_Chandra<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
               a = Ellipse(xy=(x_Xcparts_Chandra, y_Xcparts_Chandra), width=2*ell_a_Chandra, height=2*ell_b_Chandra, angle=ell_theta_Chandra)
               a.set_linewidth(0.5)
               a.set_facecolor('none')
               plt.text(x_Xcparts_Chandra,y_Xcparts_Chandra+50,'chandra',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.2*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax1.add_artist(a)
      if os.path.exists(path+'/XMM4d13s.cat'):
         if np.size(x_Xcparts_XMM)!=1:
            for j in range(0,np.size(x_Xcparts_XMM)):
               if x_Xcparts_XMM[j]>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Xcparts_XMM[j]<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Xcparts_XMM[j]>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Xcparts_XMM[j]<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
                  a = Circle(xy=(x_Xcparts_XMM[j], y_Xcparts_XMM[j]), radius=err_Xcparts_XMM[j])
                  a.set_linewidth(0.5)
                  a.set_facecolor('none')
                  plt.text(x_Xcparts_XMM[j],y_Xcparts_XMM[j]+50,'XMM',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.2*subsize,alpha=1.0)
                  a.set_edgecolor('fuchsia')
                  ax1.add_artist(a)
         else:
            if x_Xcparts_XMM>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Xcparts_XMM<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Xcparts_XMM>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Xcparts_XMM<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
               a = Circle(xy=(x_Xcparts_XMM, y_Xcparts_XMM), radius=err_Xcparts_XMM)
               a.set_linewidth(0.5)
               a.set_facecolor('none')
               plt.text(x_Xcparts_XMM,y_Xcparts_XMM+50,'XMM',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.2*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax1.add_artist(a)
      if os.path.exists(path+'/eRASS1d1.cat'):
         if np.size(x_Xcparts_eRASS)!=1:
            for j in range(0,np.size(x_Xcparts_eRASS)):
               if x_Xcparts_eRASS[j]>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Xcparts_eRASS[j]<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Xcparts_eRASS[j]>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Xcparts_eRASS[j]<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
                  a = Circle(xy=(x_Xcparts_eRASS[j], y_Xcparts_eRASS[j]), radius=err_Xcparts_eRASS[j])
                  a.set_linewidth(0.5)
                  a.set_facecolor('none')
                  plt.text(x_Xcparts_eRASS[j],y_Xcparts_eRASS[j]+50,'eRASS',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.2*subsize,alpha=1.0)
                  a.set_edgecolor('fuchsia')
                  ax1.add_artist(a)
         else:
            if x_Xcparts_eRASS>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Xcparts_eRASS<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Xcparts_eRASS>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Xcparts_eRASS<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
               a = Circle(xy=(x_Xcparts_eRASS, y_Xcparts_eRASS), radius=err_Xcparts_eRASS)
               a.set_linewidth(0.5)
               a.set_facecolor('none')
               plt.text(x_Xcparts_eRASS,y_Xcparts_eRASS+50,'eRASS',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.2*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax1.add_artist(a)
      if np.size(x_Rcparts) > 0:
         if np.size(x_Rcparts) > 1:
            for j in range(0,np.size(x_Rcparts)):
               if x_Rcparts[j]>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Rcparts[j]<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Rcparts[j]>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Rcparts[j]<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
                  if not np.isnan(errmaj_Rcparts[j]) and not np.isnan(errmin_Rcparts[j]) and not np.isnan(errangle_Rcparts[j]):
                     a = Ellipse(xy=(x_Rcparts[j], y_Rcparts[j]), width=2*errmaj_Rcparts[j], height=2*errmin_Rcparts[j], angle=errangle_Rcparts[j])
                  elif not np.isnan(errmaj_Rcparts[j]) and not np.isnan(errmin_Rcparts[j]) and np.isnan(errangle_Rcparts[j]):
                     a = Circle(xy=(x_Rcparts[j], y_Rcparts[j]),radius=errmaj_Rcparts[j])
                  else:
                     a = Circle(xy=(x_Rcparts[j], y_Rcparts[j]),radius=2.)
                  a.set_linewidth(0.5)
                  a.set_facecolor('none')
                  plt.text(x_Rcparts[j],y_Rcparts[j]+50,cat_R[j],horizontalalignment='center',verticalalignment='center',color='limegreen',fontsize=1.2*subsize,alpha=1.0)
                  a.set_edgecolor('limegreen')
                  ax1.add_artist(a)
         else:
            if x_Rcparts>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Rcparts<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Rcparts>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Rcparts<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
               if not np.isnan(errmaj_Rcparts) and not np.isnan(errmin_Rcparts) and not np.isnan(errangle_Rcparts):
                  a = Ellipse(xy=(x_Rcparts, y_Rcparts), width=2*errmaj_Rcparts, height=2*errmin_Rcparts, angle=errangle_Rcparts)
               elif not np.isnan(errmaj_Rcparts) and not np.isnan(errmin_Rcparts) and np.isnan(errangle_Rcparts):
                  a = Circle(xy=(x_Rcparts, y_Rcparts),radius=errmaj_Rcparts)
               else:
                  a = Circle(xy=(x_Rcparts, y_Rcparts),radius=2.)
               a.set_linewidth(0.5)
               a.set_facecolor('none')
               plt.text(x_Rcparts,y_Rcparts+50,cat_R,horizontalalignment='center',verticalalignment='center',color='limegreen',fontsize=1.2*subsize,alpha=1.0)
               a.set_edgecolor('limegreen')
               ax1.add_artist(a)

      if 'ra_query' in locals() and 'dec_query' in locals():
         x_query, y_query = w.all_world2pix(ra_query, dec_query, 1)
         plt.text(x_query,y_query+50,"query",horizontalalignment='center',verticalalignment='center',color='darkorange',fontsize=1.2*subsize,alpha=1.0)
         a = Circle(xy=(x_query, y_query), radius=8)
         a.set_linewidth(0.5)
         a.set_facecolor('none')
         a.set_edgecolor('darkorange')
         ax1.add_artist(a)
         

      ax2 = plt.subplot(5, 2, 10)
      plt.title('Zoom on the variable source: (RA, DEC) = (%.3f'%ra+', %.3f'%dec+') deg',fontsize=subsize)
      plt.xlim(x-zoom,x+zoom)
      plt.ylim(y-zoom,y+zoom)
      plt.imshow(data, interpolation='nearest', cmap='gray', vmin=m-colscale*s, vmax=m+colscale*s, origin='lower')
      rxnew = [rx[j] for j in range(0,len(rx)) if rx[j]>=(x-zoom) and rx[j]<=(x+zoom) and ry[j]>=(y-zoom) and ry[j]<=(y+zoom)]
      rynew = [ry[j] for j in range(0,len(rx)) if rx[j]>=(x-zoom) and rx[j]<=(x+zoom) and ry[j]>=(y-zoom) and ry[j]<=(y+zoom)]
      for j in range(len(rxnew)):
         if not(np.round(x,ndig)==np.round(rxnew[j],ndig) and np.round(y,ndig)==np.round(rynew[j],ndig)):
            a = Circle(xy=(rxnew[j], rynew[j]), radius=fext*fwhmavg)
            a.set_linewidth(1.0)
            a.set_facecolor('none')
            a.set_edgecolor('blue')
            ax2.add_artist(a)
      plt.text(x,y+20,str(i+1),horizontalalignment='center',verticalalignment='center',color='r',fontsize=2*subsize)
      a = Circle(xy=(x, y), radius=fext*fwhmavg)
      a.set_linewidth(1.0)
      a.set_facecolor('none')
      a.set_edgecolor('red')
      ax2.add_artist(a)
      #if hdr['INSTRUME']=='WIFSIP V2.0':
      #   if np.size(obs_Rcparts)!=0:
      #      xnew_Rcparts = [x_Rcparts[j] for j in range(0,len(x_Rcparts)) if x_Rcparts[j]>=(x-zoom) and x_Rcparts[j]<=(x+zoom) and y_Rcparts[j]>=(y-zoom) and y_Rcparts[j]<=(y+zoom)]
      #      ynew_Rcparts = [y_Rcparts[j] for j in range(0,len(x_Rcparts)) if x_Rcparts[j]>=(x-zoom) and x_Rcparts[j]<=(x+zoom) and y_Rcparts[j]>=(y-zoom) and y_Rcparts[j]<=(y+zoom)]
      #      obsnew_Rcparts = [obs_Rcparts[j] for j in range(0,len(x_Rcparts)) if x_Rcparts[j]>=(x-zoom) and x_Rcparts[j]<=(x+zoom) and y_Rcparts[j]>=(y-zoom) and y_Rcparts[j]<=(y+zoom)]
      #      for j in range(0,np.size(xnew_Rcparts)):
      #         a = Circle(xy=(xnew_Rcparts[j], ynew_Rcparts[j]), radius=8)
      #         a.set_linewidth(1.0)
      #         a.set_facecolor('none')
      #         plt.text(xnew_Rcparts[j],ynew_Rcparts[j],obsnew_Rcparts[j],horizontalalignment='center',verticalalignment='center',color='limegreen',fontsize=1.5*subsize,alpha=1.0)
      #         a.set_edgecolor('limegreen')
      #         ax2.add_artist(a)
      if os.path.exists(path+'/2sxps.cat'):
         if np.size(x_Xcparts_Swift)!=1:
            xnew_Xcparts = [x_Xcparts_Swift[j] for j in range(0,len(x_Xcparts_Swift)) if x_Xcparts_Swift[j]>=(x-zoom) and x_Xcparts_Swift[j]<=(x+zoom) and y_Xcparts_Swift[j]>=(y-zoom) and y_Xcparts_Swift[j]<=(y+zoom)]
            ynew_Xcparts = [y_Xcparts_Swift[j] for j in range(0,len(x_Xcparts_Swift)) if x_Xcparts_Swift[j]>=(x-zoom) and x_Xcparts_Swift[j]<=(x+zoom) and y_Xcparts_Swift[j]>=(y-zoom) and y_Xcparts_Swift[j]<=(y+zoom)]
            errnew_Xcparts = [err_Xcparts_Swift[j] for j in range(0,len(x_Xcparts_Swift)) if x_Xcparts_Swift[j]>=(x-zoom) and x_Xcparts_Swift[j]<=(x+zoom) and y_Xcparts_Swift[j]>=(y-zoom) and y_Xcparts_Swift[j]<=(y+zoom)]
            for j in range(0,np.size(xnew_Xcparts)):
               a = Circle(xy=(xnew_Xcparts[j], ynew_Xcparts[j]), radius=errnew_Xcparts[j])
               a.set_linewidth(1.0)
               a.set_facecolor('none')
               plt.text(xnew_Xcparts[j],ynew_Xcparts[j],'swift',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.5*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax2.add_artist(a)
         else:
            if x_Xcparts_Swift>=(x-zoom) and x_Xcparts_Swift<=(x+zoom) and y_Xcparts_Swift>=(y-zoom) and y_Xcparts_Swift<=(y+zoom):
               xnew_Xcparts = x_Xcparts_Swift
               ynew_Xcparts = y_Xcparts_Swift
               errnew_Xcparts = err_Xcparts_Swift
               a = Circle(xy=(xnew_Xcparts, ynew_Xcparts), radius=errnew_Xcparts)
               a.set_linewidth(1.0)
               a.set_facecolor('none')
               plt.text(xnew_Xcparts,ynew_Xcparts,'swift',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.5*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax2.add_artist(a)
      if os.path.exists(path+'/csc2.cat'):
         if np.size(x_Xcparts_Chandra)!=1:
            xnew_Xcparts = [x_Xcparts_Chandra[j] for j in range(0,len(x_Xcparts_Chandra)) if x_Xcparts_Chandra[j]>=(x-zoom) and x_Xcparts_Chandra[j]<=(x+zoom) and y_Xcparts_Chandra[j]>=(y-zoom) and y_Xcparts_Chandra[j]<=(y+zoom)]
            ynew_Xcparts = [y_Xcparts_Chandra[j] for j in range(0,len(x_Xcparts_Chandra)) if x_Xcparts_Chandra[j]>=(x-zoom) and x_Xcparts_Chandra[j]<=(x+zoom) and y_Xcparts_Chandra[j]>=(y-zoom) and y_Xcparts_Chandra[j]<=(y+zoom)]
            ellanew_Xcparts = [ell_a_Chandra[j] for j in range(0,len(x_Xcparts_Chandra)) if x_Xcparts_Chandra[j]>=(x-zoom) and x_Xcparts_Chandra[j]<=(x+zoom) and y_Xcparts_Chandra[j]>=(y-zoom) and y_Xcparts_Chandra[j]<=(y+zoom)]
            ellbnew_Xcparts = [ell_b_Chandra[j] for j in range(0,len(x_Xcparts_Chandra)) if x_Xcparts_Chandra[j]>=(x-zoom) and x_Xcparts_Chandra[j]<=(x+zoom) and y_Xcparts_Chandra[j]>=(y-zoom) and y_Xcparts_Chandra[j]<=(y+zoom)]
            ellthetanew_Xcparts = [ell_theta_Chandra[j] for j in range(0,len(x_Xcparts_Chandra)) if x_Xcparts_Chandra[j]>=(x-zoom) and x_Xcparts_Chandra[j]<=(x+zoom) and y_Xcparts_Chandra[j]>=(y-zoom) and y_Xcparts_Chandra[j]<=(y+zoom)]
            for j in range(0,np.size(xnew_Xcparts)):
               a = Ellipse(xy=(xnew_Xcparts[j], ynew_Xcparts[j]), width=2*ellanew_Xcparts[j], height=2*ellbnew_Xcparts[j], angle=ellthetanew_Xcparts[j])
               a.set_linewidth(1.0)
               a.set_facecolor('none')
               plt.text(xnew_Xcparts[j],ynew_Xcparts[j],'chandra',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.5*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax2.add_artist(a)
         else:
            if x_Xcparts_Chandra>=(x-zoom) and x_Xcparts_Chandra<=(x+zoom) and y_Xcparts_Chandra>=(y-zoom) and y_Xcparts_Chandra<=(y+zoom):
               xnew_Xcparts = x_Xcparts_Chandra
               ynew_Xcparts = y_Xcparts_Chandra
               ellanew_Xcparts = ell_a_Chandra
               ellbnew_Xcparts = ell_b_Chandra
               ellthetanew_Xcparts = ell_theta_Chandra
               a = Ellipse(xy=(xnew_Xcparts, ynew_Xcparts), width=2*ellanew_Xcparts, height=2*ellbnew_Xcparts, angle=ellthetanew_Xcparts)
               a.set_linewidth(1.0)
               a.set_facecolor('none')
               plt.text(xnew_Xcparts,ynew_Xcparts,'chandra',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.5*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax2.add_artist(a)
      if os.path.exists(path+'/XMM4d13s.cat'):
         if np.size(x_Xcparts_XMM)!=1:
            xnew_Xcparts = [x_Xcparts_XMM[j] for j in range(0,len(x_Xcparts_XMM)) if x_Xcparts_XMM[j]>=(x-zoom) and x_Xcparts_XMM[j]<=(x+zoom) and y_Xcparts_XMM[j]>=(y-zoom) and y_Xcparts_XMM[j]<=(y+zoom)]
            ynew_Xcparts = [y_Xcparts_XMM[j] for j in range(0,len(x_Xcparts_XMM)) if x_Xcparts_XMM[j]>=(x-zoom) and x_Xcparts_XMM[j]<=(x+zoom) and y_Xcparts_XMM[j]>=(y-zoom) and y_Xcparts_XMM[j]<=(y+zoom)]
            errnew_Xcparts = [err_Xcparts_XMM[j] for j in range(0,len(x_Xcparts_XMM)) if x_Xcparts_XMM[j]>=(x-zoom) and x_Xcparts_XMM[j]<=(x+zoom) and y_Xcparts_XMM[j]>=(y-zoom) and y_Xcparts_XMM[j]<=(y+zoom)]
            for j in range(0,np.size(xnew_Xcparts)):
               a = Circle(xy=(xnew_Xcparts[j], ynew_Xcparts[j]), radius=errnew_Xcparts[j])
               a.set_linewidth(1.0)
               a.set_facecolor('none')
               plt.text(xnew_Xcparts[j],ynew_Xcparts[j],'XMM',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.5*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax2.add_artist(a)
         else:
            if x_Xcparts_XMM>=(x-zoom) and x_Xcparts_XMM<=(x+zoom) and y_Xcparts_XMM>=(y-zoom) and y_Xcparts_XMM<=(y+zoom):
               xnew_Xcparts = x_Xcparts_XMM
               ynew_Xcparts = y_Xcparts_XMM
               errnew_Xcparts = err_Xcparts_XMM
               a = Circle(xy=(xnew_Xcparts, ynew_Xcparts), radius=errnew_Xcparts)
               a.set_linewidth(1.0)
               a.set_facecolor('none')
               plt.text(xnew_Xcparts,ynew_Xcparts,'XMM',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.5*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax2.add_artist(a)
      if os.path.exists(path+'/eRASS1d1.cat'):
         if np.size(x_Xcparts_eRASS)!=1:
            xnew_Xcparts = [x_Xcparts_eRASS[j] for j in range(0,len(x_Xcparts_eRASS)) if x_Xcparts_eRASS[j]>=(x-zoom) and x_Xcparts_eRASS[j]<=(x+zoom) and y_Xcparts_eRASS[j]>=(y-zoom) and y_Xcparts_eRASS[j]<=(y+zoom)]
            ynew_Xcparts = [y_Xcparts_eRASS[j] for j in range(0,len(x_Xcparts_eRASS)) if x_Xcparts_eRASS[j]>=(x-zoom) and x_Xcparts_eRASS[j]<=(x+zoom) and y_Xcparts_eRASS[j]>=(y-zoom) and y_Xcparts_eRASS[j]<=(y+zoom)]
            errnew_Xcparts = [err_Xcparts_eRASS[j] for j in range(0,len(x_Xcparts_eRASS)) if x_Xcparts_eRASS[j]>=(x-zoom) and x_Xcparts_eRASS[j]<=(x+zoom) and y_Xcparts_eRASS[j]>=(y-zoom) and y_Xcparts_eRASS[j]<=(y+zoom)]
            for j in range(0,np.size(xnew_Xcparts)):
               a = Circle(xy=(xnew_Xcparts[j], ynew_Xcparts[j]), radius=errnew_Xcparts[j])
               a.set_linewidth(1.0)
               a.set_facecolor('none')
               plt.text(xnew_Xcparts[j],ynew_Xcparts[j],'eRASS',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.5*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax2.add_artist(a)
         else:
            if x_Xcparts_eRASS>=(x-zoom) and x_Xcparts_eRASS<=(x+zoom) and y_Xcparts_eRASS>=(y-zoom) and y_Xcparts_eRASS<=(y+zoom):
               xnew_Xcparts = x_Xcparts_eRASS
               ynew_Xcparts = y_Xcparts_eRASS
               errnew_Xcparts = err_Xcparts_eRASS
               a = Circle(xy=(xnew_Xcparts, ynew_Xcparts), radius=errnew_Xcparts)
               a.set_linewidth(1.0)
               a.set_facecolor('none')
               plt.text(xnew_Xcparts,ynew_Xcparts,'eRASS',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.5*subsize,alpha=1.0)
               a.set_edgecolor('fuchsia')
               ax2.add_artist(a)
      if np.size(x_Rcparts) > 0:
         if np.size(x_Rcparts) > 1:
            xnew_Rcparts = [x_Rcparts[j] for j in range(0,len(x_Rcparts)) if x_Rcparts[j]>=(x-zoom) and x_Rcparts[j]<=(x+zoom) and y_Rcparts[j]>=(y-zoom) and y_Rcparts[j]<=(y+zoom)]
            ynew_Rcparts = [y_Rcparts[j] for j in range(0,len(x_Rcparts)) if x_Rcparts[j]>=(x-zoom) and x_Rcparts[j]<=(x+zoom) and y_Rcparts[j]>=(y-zoom) and y_Rcparts[j]<=(y+zoom)]
            errmajnew_Rcparts = [errmaj_Rcparts[j] for j in range(0,len(x_Rcparts)) if x_Rcparts[j]>=(x-zoom) and x_Rcparts[j]<=(x+zoom) and y_Rcparts[j]>=(y-zoom) and y_Rcparts[j]<=(y+zoom)]
            errminnew_Rcparts = [errmin_Rcparts[j] for j in range(0,len(x_Rcparts)) if x_Rcparts[j]>=(x-zoom) and x_Rcparts[j]<=(x+zoom) and y_Rcparts[j]>=(y-zoom) and y_Rcparts[j]<=(y+zoom)]
            erranglenew_Rcparts = [errangle_Rcparts[j] for j in range(0,len(x_Rcparts)) if x_Rcparts[j]>=(x-zoom) and x_Rcparts[j]<=(x+zoom) and y_Rcparts[j]>=(y-zoom) and y_Rcparts[j]<=(y+zoom)]
            catnew_R = [cat_R[j] for j in range(0,len(x_Rcparts)) if x_Rcparts[j]>=(x-zoom) and x_Rcparts[j]<=(x+zoom) and y_Rcparts[j]>=(y-zoom) and y_Rcparts[j]<=(y+zoom)]
            for j in range(0,np.size(xnew_Rcparts)):
               if not np.isnan(errmajnew_Rcparts[j]) and not np.isnan(errminnew_Rcparts[j]) and not np.isnan(erranglenew_Rcparts[j]):
                  a = Ellipse(xy=(xnew_Rcparts[j], ynew_Rcparts[j]), width=2*errmajnew_Rcparts[j], height=2*errminnew_Rcparts[j], angle=erranglenew_Rcparts[j])
               elif not np.isnan(errmajnew_Rcparts[j]) and not np.isnan(errminnew_Rcparts[j]) and np.isnan(erranglenew_Rcparts[j]):
                  a = Circle(xy=(xnew_Rcparts[j], ynew_Rcparts[j]),radius=errmajnew_Rcparts[j])
               else:
                  a = Circle(xy=(xnew_Rcparts[j], ynew_Rcparts[j]),radius=2.)
               a.set_linewidth(1.0)
               a.set_facecolor('none')
               plt.text(xnew_Rcparts[j],ynew_Rcparts[j],catnew_R[j],horizontalalignment='center',verticalalignment='center',color='limegreen',fontsize=1.5*subsize,alpha=1.0)
               a.set_edgecolor('limegreen')
               ax2.add_artist(a)
         else:
            if x_Rcparts>=(x-zoom) and x_Rcparts<=(x+zoom) and y_Rcparts>=(y-zoom) and y_Rcparts<=(y+zoom):
               xnew_Rcparts = x_Rcparts
               ynew_Rcparts = y_Rcparts
               errmajnew_Rcparts = errmaj_Rcparts
               errminnew_Rcparts = errmin_Rcparts
               erranglenew_Rcparts = errangle_Rcparts
               catnew_R = cat_R
               if not np.isnan(errmajnew_Rcparts) and not np.isnan(errminnew_Rcparts) and not np.isnan(erranglenew_Rcparts):
                  a = Ellipse(xy=(xnew_Rcparts, ynew_Rcparts), width=2*errmajnew_Rcparts, height=2*errminnew_Rcparts, angle=erranglenew_Rcparts)
               elif not np.isnan(errmajnew_Rcparts) and not np.isnan(errminnew_Rcparts) and np.isnan(erranglenew_Rcparts):
                  a = Circle(xy=(xnew_Rcparts, ynew_Rcparts),radius=errmajnew_Rcparts)
               else:
                  a = Circle(xy=(xnew_Rcparts, ynew_Rcparts),radius=2.)
               a.set_linewidth(1.0)
               a.set_facecolor('none')
               plt.text(xnew_Rcparts,ynew_Rcparts,catnew_R,horizontalalignment='center',verticalalignment='center',color='limegreen',fontsize=1.5*subsize,alpha=1.0)
               a.set_edgecolor('limegreen')
               ax2.add_artist(a)
      
      if 'ra_query' in locals() and 'dec_query' in locals():
         plt.text(x_query,y_query,"query",horizontalalignment='center',verticalalignment='center',color='darkorange',fontsize=1.5*subsize,alpha=0.8)
         a = Circle(xy=(x_query, y_query), radius=8)
         a.set_linewidth(0.5)
         a.set_facecolor('none')
         a.set_edgecolor('darkorange')
         ax2.add_artist(a)
      
      if 'kic_ID' in locals():
         plt.text(-1.85,6.0,"Identified from Kepler as "+kic_ID+r" with $P_{\mathrm{Kepler}}$="+str(kic_period)+" d",transform=ax2.transAxes,fontsize=tsize)
      elif 'ATO_ID' in locals():
         plt.text(-1.85,6.0,"Identified from ATLAS as "+ATO_ID+r" with $P_{\mathrm{ATLAS}}$="+str(ATO_period)+" d",transform=ax2.transAxes,fontsize=tsize)
      elif 'SIMBAD_ID' in locals():
         plt.text(-1.7,6.0,"Identified from SIMBAD as "+SIMBAD_ID,transform=ax2.transAxes,fontsize=tsize)
      plt.savefig(path+'/finalperiodplots'+str(nlist)+'_'+per_search+'/finalperplot_'+str(i+1)+'_'+per_search+'_ZTF_nper'+str(int(n_periods))+'.png',bbox_inches='tight')
      plt.close(figure)

   if ATLAS=='yes':
      #Plot ATLAS lightcurves folded at best_periods[i] and t_0 in c and o ATLAS filters together with Lomb-Scargle periodograms obtained with forced photometry (https://fallingstar-data.com/forcedphot)
      BASEURL = "https://fallingstar-data.com/forcedphot"
      resp = requests.post(url=f"{BASEURL}/api-token-auth/", data={'username': "marcotu", 'password': "Marco6813!"})
      if resp.status_code == 200:
         token = resp.json()['token']
         print(f'Your token is {token}')
         headers = {'Authorization': f'Token {token}', 'Accept': 'application/json'}
      else:
         print(f'ERROR {resp.status_code}')
         print(resp.json())

      task_url = None
      while not task_url:
         with requests.Session() as ses:
            resp = ses.post(f"{BASEURL}/queue/", headers=headers, data={
               'ra': ra, 'dec': dec, 'mjd_min': 59215., 'send_email': False})

            if resp.status_code == 201:  # successfully queued
               task_url = resp.json()['url']
               print(f'The task URL is {task_url}')
            elif resp.status_code == 429:  # throttled
               message = resp.json()["detail"]
               print(f'{resp.status_code} {message}')
               t_sec = re.findall(r'available in (\d+) seconds', message)
               t_min = re.findall(r'available in (\d+) minutes', message)
               if t_sec:
                  waittime = int(t_sec[0])
               elif t_min:
                  waittime = int(t_min[0]) * 60
               else:
                  waittime = 10
                  print(f'Waiting {waittime} seconds')
                  time.sleep(waittime)
            else:
               print(f'ERROR {resp.status_code}')
               print(resp.json())

      result_url = None
      while not result_url:
         with requests.Session() as ses:
            resp = ses.get(task_url, headers=headers)

            if resp.status_code == 200:  # HTTP OK
               if resp.json()['finishtimestamp']:
                  result_url = resp.json()['result_url']
                  print(f"Task is complete with results available at {result_url}")
                  break
               elif resp.json()['starttimestamp']:
                  print(f"Task is running (started at {resp.json()['starttimestamp']})")
               else:
                  print("Waiting for job to start. Checking again in 10 seconds...")
               time.sleep(10)
            else:
               print(f'ERROR {resp.status_code}')
               print(resp.json())

      with requests.Session() as ses:
         textdata = ses.get(result_url, headers=headers).text
      dfresult = pd.read_csv(io.StringIO(textdata.replace("###", "")), delim_whitespace=True)
      mjd_ATLAS = dfresult['MJD'].to_numpy()
      uJy_ATLAS = dfresult['uJy'].to_numpy()
      duJy_ATLAS = dfresult['duJy'].to_numpy()
      filter_ATLAS = dfresult['F'].to_numpy()

      #Filter out bad data points with negative values or errors above 100 muJy
      mjd_ATLAS = [mjd_ATLAS[i] for i in range(0,len(uJy_ATLAS)) if uJy_ATLAS[i]>0.]
      filter_ATLAS = [filter_ATLAS[i] for i in range(0,len(uJy_ATLAS)) if uJy_ATLAS[i]>0.]
      duJy_ATLAS = [duJy_ATLAS[i] for i in range(0,len(uJy_ATLAS)) if uJy_ATLAS[i]>0.]
      uJy_ATLAS = [uJy_ATLAS[i] for i in range(0,len(uJy_ATLAS)) if uJy_ATLAS[i]>0.]
      mjd_ATLAS = [mjd_ATLAS[i] for i in range(0,len(duJy_ATLAS)) if duJy_ATLAS[i]<100.]
      filter_ATLAS = [filter_ATLAS[i] for i in range(0,len(duJy_ATLAS)) if duJy_ATLAS[i]<100.]
      uJy_ATLAS = [uJy_ATLAS[i] for i in range(0,len(duJy_ATLAS)) if duJy_ATLAS[i]<100.]
      duJy_ATLAS = [duJy_ATLAS[i] for i in range(0,len(duJy_ATLAS)) if duJy_ATLAS[i]<100.]

      mjd_ATLASc = np.array([mjd_ATLAS[i] for i in range(0,len(mjd_ATLAS)) if filter_ATLAS[i] == 'c'])
      uJy_ATLASc = np.array([uJy_ATLAS[i] for i in range(0,len(mjd_ATLAS)) if filter_ATLAS[i] == 'c'])
      duJy_ATLASc = np.array([duJy_ATLAS[i] for i in range(0,len(mjd_ATLAS)) if filter_ATLAS[i] == 'c'])

      mjd_ATLASo = np.array([mjd_ATLAS[i] for i in range(0,len(mjd_ATLAS)) if filter_ATLAS[i] == 'o'])
      uJy_ATLASo = np.array([uJy_ATLAS[i] for i in range(0,len(mjd_ATLAS)) if filter_ATLAS[i] == 'o'])
      duJy_ATLASo = np.array([duJy_ATLAS[i] for i in range(0,len(mjd_ATLAS)) if filter_ATLAS[i] == 'o'])

      del textdata, dfresult, mjd_ATLAS, uJy_ATLAS, duJy_ATLAS, filter_ATLAS

      if np.size(mjd_ATLASc)!=0 or np.size(mjd_ATLASo)!=0:
         figure=plt.figure(figsize=(figxsize,figysize))

         if 'flux_Xcpart' in locals():
            figure.suptitle(r"Fermi source: "+FoVname+", Periodic variable #%d"%(i+1)+r", Period folding with $P_{\mathrm{best}}=$%.4f"%best_periods[i]+" d and $T_{0}=$%.4f"%t_0+" MJD\n\nRange of the Lomb-Scargle periodogram: [%.2f"%P_low+",%.1f"%min(P_up,(mjd_ATLASo[-1]-mjd_ATLASo[0]))+"] d, Steps: %d"%n_periods+"\n\nX-ray flux: "+str(np.round(float(str(flux_Xcpart).split('e')[0]),1))+"$\cdot 10^{"+str(flux_Xcpart).split('e')[1]+"}$ erg s$^{-1}$ cm$^{-2}$", fontsize=tsize)
         else:
            if 'UL_X' in locals():
               figure.suptitle("Fermi source: "+FoVname+", Periodic variable #%d"%(i+1)+r", Period folding with $P_{\mathrm{best}}=$%.4f"%best_periods[i]+" d and $T_{0}=$%.4f"%t_0+" MJD\n\nRange of the Lomb-Scargle periodogram: [%.2f"%P_low+",%.1f"%min(P_up,(mjd_ATLASo[-1]-mjd_ATLASo[0]))+"] d, Steps: %d"%n_periods+"\n\nX-ray flux 3$\sigma$ UL: "+str(np.round(float(str(UL_X).split('e')[0]),1))+"$\cdot 10^{"+str(UL_X).split('e')[1]+"}$ erg s$^{-1}$ cm$^{-2}$", fontsize=tsize)
            else:
               figure.suptitle("Fermi source: "+FoVname+", Periodic variable #%d"%(i+1)+r", Period folding with $P_{\mathrm{best}}=$%.4f"%best_periods[i]+" d and $T_{0}=$%.4f"%t_0+" MJD\n\nRange of the Lomb-Scargle periodogram: [%.2f"%P_low+",%.1f"%min(P_up,(mjd_ATLASo[-1]-mjd_ATLASo[0]))+"] d, Steps: %d"%n_periods, fontsize=tsize)

         if np.size(mjd_ATLASc)!=0:
            phase_ATLASc = ((mjd_ATLASc-t_0)/best_periods[i]) % 1
            plt.subplot(4, 2, 1)
            plt.plot(phase_ATLASc, uJy_ATLASc, 'bo', linestyle='None')
            plt.plot(phase_ATLASc+1, uJy_ATLASc, 'ro', linestyle='None')
            plt.errorbar(phase_ATLASc, uJy_ATLASc, duJy_ATLASc, linestyle='None')
            plt.errorbar(phase_ATLASc+1, uJy_ATLASc, duJy_ATLASc, linestyle='None')
            plt.xticks(fontsize=axticksize)
            plt.yticks(fontsize=axticksize)
            plt.title(r"Folded ATLAS lightcurve in the c band: $P$="+"{0} d".format(best_periods[i]),fontsize=subsize)
            plt.xlabel(r"Phase ($\phi$)",fontsize=axlabelsize)
            plt.ylabel(r"Flux ($\mu$Jy)",fontsize=axlabelsize)

         if np.size(mjd_ATLASo)!=0:
            phase_ATLASo = ((mjd_ATLASo-t_0)/best_periods[i]) % 1
            plt.subplot(4, 2, 2)
            plt.plot(phase_ATLASo, uJy_ATLASo, 'bo', linestyle='None')
            plt.plot(phase_ATLASo+1, uJy_ATLASo, 'ro', linestyle='None')
            plt.errorbar(phase_ATLASo, uJy_ATLASo, duJy_ATLASo, linestyle='None')
            plt.errorbar(phase_ATLASo+1, uJy_ATLASo, duJy_ATLASo, linestyle='None')
            plt.xticks(fontsize=axticksize)
            plt.yticks(fontsize=axticksize)
            plt.title(r"Folded ATLAS lightcurve in the o band: $P$="+"{0} d".format(best_periods[i]),fontsize=subsize)
            plt.xlabel(r"Phase ($\phi$)",fontsize=axlabelsize)
            plt.ylabel(r"Flux ($\mu$Jy)",fontsize=axlabelsize)

         if np.size(mjd_ATLASc)!=0:
            if per_search == 'LS':
               periods, periodogram = LombScargle_search(mjd_ATLASc, uJy_ATLASc, duJy_ATLASc, P_low, P_up, n_periods, 1)
               best_period_ATLASc = periods[~np.isnan(periodogram)][np.argmax(periodogram[~np.isnan(periodogram)])]
            else:
               if per_search == 'PDM':
                  periods, periodogram = PDM_search(mjd_ATLASc, uJy_ATLASc, P_low, P_up, 0., n_periods, n_bin, min_bin)
               else:
                  if np.size(mjd_ATLASc)!=0 and np.size(mjd_ATLASo)!=0:
                     periods, periodogram = PDM_bisearch(mjd_ATLASc, uJy_ATLASc, mjd_ATLASo, uJy_ATLASo, P_low, P_up, 0., n_periods, n_bin, min_bin)
                  else:
                     periods, periodogram = PDM_search(mjd_ATLASc, uJy_ATLASc, P_low, P_up, 0., n_periods, n_bin, min_bin)
               best_period_ATLASc = periods[~np.isnan(periodogram)][np.argmin(periodogram[~np.isnan(periodogram)])]
         
            phase_ATLASc = ((mjd_ATLASc-t_0)/best_period_ATLASc) % 1
            plt.subplot(4, 2, 3)
            plt.plot(phase_ATLASc, uJy_ATLASc, 'bo', linestyle='None')
            plt.plot(phase_ATLASc+1, uJy_ATLASc, 'ro', linestyle='None')
            plt.errorbar(phase_ATLASc, uJy_ATLASc, duJy_ATLASc, linestyle='None')
            plt.errorbar(phase_ATLASc+1, uJy_ATLASc, duJy_ATLASc, linestyle='None')
            plt.xticks(fontsize=axticksize)
            plt.yticks(fontsize=axticksize)
            plt.title(r"Folded ATLAS lightcurve in the c band: $P_{\mathrm{ATLAS,c}}$="+"{0} d".format(best_period_ATLASc),fontsize=subsize)
            plt.xlabel(r"Phase ($\phi$)",fontsize=axlabelsize)
            plt.ylabel(r"Flux ($\mu$Jy)",fontsize=axlabelsize)

            plt.subplot(4, 2, 4)
            plt.plot(periods,periodogram)
            if not np.isnan(max(periodogram)) and not np.isinf(max(periodogram)) and not np.isnan(min(periodogram)) and not np.isinf(min(periodogram)) and not np.isnan(best_period_ATLASc) and not np.isinf(best_period_ATLASc):
               plt.axvline(x=best_period_ATLASc, color='g', ls='--')
            plt.axvline(x=best_periods[i], color='r', ls='--')
            plt.xscale("log")
            plt.xticks(fontsize=axticksize)
            plt.yticks(fontsize=axticksize)
            if per_search == 'LS':
               plt.title("Lomb-Scargle ATLAS data periodogram in the c band",fontsize=subsize)
               plt.ylabel(r"Power",fontsize=axlabelsize)
            elif per_search == 'PDM':
               plt.title("PDM ATLAS data periodogram in the c band",fontsize=subsize)
               plt.ylabel(r"Normalized variance",fontsize=axlabelsize)
            else:
               plt.title("PDM multi-band ATLAS data periodogram in the c band",fontsize=subsize)
               plt.ylabel(r"Sum of normalized variances",fontsize=axlabelsize)
            plt.xlabel(r"Trial period (days)",fontsize=axlabelsize)

         if np.size(mjd_ATLASo)!=0:
            if per_search == 'LS':
               periods, periodogram = LombScargle_search(mjd_ATLASo, uJy_ATLASo, duJy_ATLASo, P_low, P_up, n_periods, 1)
               best_period_ATLASo = periods[~np.isnan(periodogram)][np.argmax(periodogram[~np.isnan(periodogram)])]
            else:
               if per_search == 'PDM':
                  periods, periodogram = PDM_search(mjd_ATLASo, uJy_ATLASo, P_low, P_up, 0., n_periods, n_bin, min_bin)
               else:
                  if np.size(mjd_ATLASc)!=0 and np.size(mjd_ATLASo)!=0:
                     periods, periodogram = PDM_bisearch(mjd_ATLASc, uJy_ATLASc, mjd_ATLASo, uJy_ATLASo, P_low, P_up, 0., n_periods, n_bin, min_bin)
                  else:
                     periods, periodogram = PDM_search(mjd_ATLASo, uJy_ATLASo, P_low, P_up, 0., n_periods, n_bin, min_bin)
               best_period_ATLASo = periods[~np.isnan(periodogram)][np.argmin(periodogram[~np.isnan(periodogram)])]
         
            phase_ATLASo = ((mjd_ATLASo-t_0)/best_period_ATLASo) % 1
            plt.subplot(4, 2, 5)
            plt.plot(phase_ATLASo, uJy_ATLASo, 'bo', linestyle='None')
            plt.plot(phase_ATLASo+1, uJy_ATLASo, 'ro', linestyle='None')
            plt.errorbar(phase_ATLASo, uJy_ATLASo, duJy_ATLASo, linestyle='None')
            plt.errorbar(phase_ATLASo+1, uJy_ATLASo, duJy_ATLASo, linestyle='None')
            plt.xticks(fontsize=axticksize)
            plt.yticks(fontsize=axticksize)
            plt.title(r"Folded ATLAS lightcurve in the o band: $P_{\mathrm{ATLAS,o}}$="+"{0} d".format(best_period_ATLASo),fontsize=subsize)
            plt.xlabel(r"Phase ($\phi$)",fontsize=axlabelsize)
            plt.ylabel(r"Flux ($\mu$Jy)",fontsize=axlabelsize)

            plt.subplot(4, 2, 6)
            plt.plot(periods,periodogram)
            if not np.isnan(max(periodogram)) and not np.isinf(max(periodogram)) and not np.isnan(min(periodogram)) and not np.isinf(min(periodogram)) and not np.isnan(best_period_ATLASo) and not np.isinf(best_period_ATLASo):
               plt.axvline(x=best_period_ATLASo, color='g', ls='--')
            plt.axvline(x=best_periods[i], color='r', ls='--')
            plt.xscale("log")
            plt.xticks(fontsize=axticksize)
            plt.yticks(fontsize=axticksize)
            if per_search == 'LS':
               plt.title("Lomb-Scargle ATLAS data periodogram in the o band",fontsize=subsize)
               plt.ylabel(r"Power",fontsize=axlabelsize)
            elif per_search == 'PDM':
               plt.title("PDM ATLAS data periodogram in the o band",fontsize=subsize)
               plt.ylabel(r"Normalized variance",fontsize=axlabelsize)
            else:
               plt.title("PDM multi-band ATLAS data periodogram in the o band",fontsize=subsize)
               plt.ylabel(r"Sum of normalized variances",fontsize=axlabelsize)
            plt.xlabel(r"Trial period (days)",fontsize=axlabelsize)

         ax1 = plt.subplot(4, 2, 7)
         plt.title(r'Median r FoV bkg subtracted in the '+str(2*bound)+'$\sigma$ square around the 4FGL coordinates',fontsize=subsize)
         plt.xlim(max(0,xtar-bound*max(ell_a_95,ell_b_95)), min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)))
         plt.ylim(max(0,ytar-bound*max(ell_a_95,ell_b_95)), min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)))
         plt.imshow(data, interpolation='nearest', cmap='gray', vmin=m-colscale*s, vmax=m+colscale*s, origin='lower')
         plt.text(xtar,ytar,FoVname,horizontalalignment='center',verticalalignment='center',color='cyan',fontsize=subsize,alpha=0.8)
         target = Ellipse(xy=(xtar, ytar), width=2*ell_a_95, height=2*ell_b_95, angle=ell_theta_95)
         target.set_facecolor('none')
         target.set_edgecolor('cyan')
         ax1.add_artist(target)
         if np.size(sources_3FGL)!=0:
            for j in range(0,np.shape(sources_3FGL)[0]):
               if x_3FGL[j]>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_3FGL[j]<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_3FGL[j]>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_3FGL[j]<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
                  plt.text(x_3FGL[j],y_3FGL[j],sources_3FGL[j].replace("L ","L"),horizontalalignment='center',verticalalignment='center',color='yellow',fontsize=subsize,alpha=0.8)
                  a = Ellipse(xy=(x_3FGL[j], y_3FGL[j]), width=2*3600./float(platescale)*coord_3FGL[j,3], height=2*3600./float(platescale)*coord_3FGL[j,4], angle=coord_3FGL[j,5])
                  a.set_facecolor('none')
                  a.set_edgecolor('yellow')
                  ax1.add_artist(a)
         #if np.size(sources_4FGLe)!=0:
         #   for j in range(0,np.shape(sources_4FGLe)[0]):
         #      if x_4FGLe[j]>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_4FGLe[j]<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_4FGLe[j]>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_4FGLe[j]<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
         #         plt.text(x_4FGLe[j],y_4FGLe[j],sources_4FGLe[j].replace("L ","L"),horizontalalignment='center',verticalalignment='center',color='cyan',fontsize=subsize,alpha=0.8)
         #         a = Ellipse(xy=(x_4FGLe[j], y_4FGLe[j]), width=2*3600./float(platescale)*np.asarray(coord_4FGLe)[j,2], height=2*3600./float(platescale)*np.asarray(coord_4FGLe)[j,3], angle=np.asarray(coord_4FGLe)[j,4])
         #         a.set_facecolor('none')
         #         a.set_edgecolor('cyan')
         #         ax1.add_artist(a)
         
         x, y = w.all_world2pix(ra, dec, 1)
         plt.text(x,y+50,str(i+1),horizontalalignment='center',verticalalignment='center',color='r',fontsize=subsize)
         a = Circle(xy=(x, y), radius=fext*fwhmavg)
         a.set_linewidth(0.5)
         a.set_facecolor('none')
         a.set_edgecolor('red')
         ax1.add_artist(a)
         if hdr['INSTRUME']=='WIFSIP V2.0':
            if np.size(obs_Rcparts)!=0:
               for j in range(0,np.size(obs_Rcparts)):
                  if x_Rcparts[j]>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Rcparts[j]<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Rcparts[j]>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Rcparts[j]<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
                     a = Circle(xy=(x_Rcparts[j], y_Rcparts[j]), radius=8)
                     a.set_linewidth(0.5)
                     a.set_facecolor('none')
                     plt.text(x_Rcparts[j],y_Rcparts[j]+50,obs_Rcparts[j],horizontalalignment='center',verticalalignment='center',color='limegreen',fontsize=1.2*subsize,alpha=0.8)
                     a.set_edgecolor('limegreen')
                     ax1.add_artist(a)
         if os.path.exists(path+'/2sxps.cat'):
            if np.size(x_Xcparts_Swift)!=1:
               for j in range(0,np.size(x_Xcparts_Swift)):
                  if x_Xcparts_Swift[j]>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Xcparts_Swift[j]<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Xcparts_Swift[j]>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Xcparts_Swift[j]<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
                     a = Circle(xy=(x_Xcparts_Swift[j], y_Xcparts_Swift[j]), radius=err_Xcparts_Swift[j])
                     a.set_linewidth(0.5)
                     a.set_facecolor('none')
                     plt.text(x_Xcparts_Swift[j],y_Xcparts_Swift[j]+50,'swift',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.2*subsize,alpha=0.8)
                     a.set_edgecolor('fuchsia')
                     ax1.add_artist(a)
            else:
               if x_Xcparts_Swift>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Xcparts_Swift<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Xcparts_Swift>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Xcparts_Swift<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
                  a = Circle(xy=(x_Xcparts_Swift, y_Xcparts_Swift), radius=err_Xcparts_Swift)
                  a.set_linewidth(0.5)
                  a.set_facecolor('none')
                  plt.text(x_Xcparts_Swift,y_Xcparts_Swift+50,'swift',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.2*subsize,alpha=0.8)
                  a.set_edgecolor('fuchsia')
                  ax1.add_artist(a)
         if os.path.exists(path+'/csc2.cat'):
            if np.size(x_Xcparts_Chandra)!=1:
               for j in range(0,np.size(x_Xcparts_Chandra)):
                  if x_Xcparts_Chandra[j]>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Xcparts_Chandra[j]<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Xcparts_Chandra[j]>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Xcparts_Chandra[j]<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
                     a = Ellipse(xy=(x_Xcparts_Chandra[j], y_Xcparts_Chandra[j]), width=2*ell_a_Chandra[j], height=2*ell_b_Chandra[j], angle=ell_theta_Chandra[j])
                     a.set_linewidth(0.5)
                     a.set_facecolor('none')
                     plt.text(x_Xcparts_Chandra[j],y_Xcparts_Chandra[j]+50,'chandra',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.2*subsize,alpha=0.8)
                     a.set_edgecolor('fuchsia')
                     ax1.add_artist(a)
            else:
               if x_Xcparts_Chandra>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Xcparts_Chandra<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Xcparts_Chandra>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Xcparts_Chandra<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
                  a = Ellipse(xy=(x_Xcparts_Chandra, y_Xcparts_Chandra), width=2*ell_a_Chandra, height=2*ell_b_Chandra, angle=ell_theta_Chandra)
                  a.set_linewidth(0.5)
                  a.set_facecolor('none')
                  plt.text(x_Xcparts_Chandra,y_Xcparts_Chandra+50,'chandra',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.2*subsize,alpha=0.8)
                  a.set_edgecolor('fuchsia')
                  ax1.add_artist(a)
         if os.path.exists(path+'/XMM4d13s.cat'):
            if np.size(x_Xcparts_XMM)!=1:
               for j in range(0,np.size(x_Xcparts_XMM)):
                  if x_Xcparts_XMM[j]>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Xcparts_XMM[j]<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Xcparts_XMM[j]>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Xcparts_XMM[j]<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
                     a = Circle(xy=(x_Xcparts_XMM[j], y_Xcparts_XMM[j]), radius=err_Xcparts_XMM[j])
                     a.set_linewidth(0.5)
                     a.set_facecolor('none')
                     plt.text(x_Xcparts_XMM[j],y_Xcparts_XMM[j]+50,'XMM',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.2*subsize,alpha=0.8)
                     a.set_edgecolor('fuchsia')
                     ax1.add_artist(a)
            else:
               if x_Xcparts_XMM>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_Xcparts_XMM<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_Xcparts_XMM>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_Xcparts_XMM<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
                  a = Circle(xy=(x_Xcparts_XMM, y_Xcparts_XMM), radius=err_Xcparts_XMM)
                  a.set_linewidth(0.5)
                  a.set_facecolor('none')
                  plt.text(x_Xcparts_XMM,y_Xcparts_XMM+50,'XMM',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.2*subsize,alpha=0.8)
                  a.set_edgecolor('fuchsia')
                  ax1.add_artist(a)

         if 'ra_query' in locals() and 'dec_query' in locals():
            x_query, y_query = w.all_world2pix(ra_query, dec_query, 1)
            plt.text(x_query,y_query+50,"query",horizontalalignment='center',verticalalignment='center',color='darkorange',fontsize=1.2*subsize,alpha=0.8)
            a = Circle(xy=(x_query, y_query), radius=8)
            a.set_linewidth(0.5)
            a.set_facecolor('none')
            a.set_edgecolor('darkorange')
            ax1.add_artist(a)
         

         ax2 = plt.subplot(int(4, 2, 8))
         plt.title('Zoom on the variable source: (RA, DEC) = (%.3f'%ra+', %.3f'%dec+') deg',fontsize=subsize)
         plt.xlim(x-zoom,x+zoom)
         plt.ylim(y-zoom,y+zoom)
         plt.imshow(data, interpolation='nearest', cmap='gray', vmin=m-colscale*s, vmax=m+colscale*s, origin='lower')
         rxnew = [rx[j] for j in range(0,len(rx)) if rx[j]>=(x-zoom) and rx[j]<=(x+zoom) and ry[j]>=(y-zoom) and ry[j]<=(y+zoom)]
         rynew = [ry[j] for j in range(0,len(rx)) if rx[j]>=(x-zoom) and rx[j]<=(x+zoom) and ry[j]>=(y-zoom) and ry[j]<=(y+zoom)]
         for j in range(len(rxnew)):
            if not(np.round(x,ndig)==np.round(rxnew[j],ndig) and np.round(y,ndig)==np.round(rynew[j],ndig)):
               a = Circle(xy=(rxnew[j], rynew[j]), radius=fext*fwhmavg)
               a.set_linewidth(1.0)
               a.set_facecolor('none')
               a.set_edgecolor('blue')
               ax2.add_artist(a)
         plt.text(x,y+20,str(i+1),horizontalalignment='center',verticalalignment='center',color='r',fontsize=2*subsize)
         a = Circle(xy=(x, y), radius=fext*fwhmavg)
         a.set_linewidth(1.0)
         a.set_facecolor('none')
         a.set_edgecolor('red')
         ax2.add_artist(a)
         if hdr['INSTRUME']=='WIFSIP V2.0':
            if np.size(obs_Rcparts)!=0:
               xnew_Rcparts = [x_Rcparts[j] for j in range(0,len(x_Rcparts)) if x_Rcparts[j]>=(x-zoom) and x_Rcparts[j]<=(x+zoom) and y_Rcparts[j]>=(y-zoom) and y_Rcparts[j]<=(y+zoom)]
               ynew_Rcparts = [y_Rcparts[j] for j in range(0,len(x_Rcparts)) if x_Rcparts[j]>=(x-zoom) and x_Rcparts[j]<=(x+zoom) and y_Rcparts[j]>=(y-zoom) and y_Rcparts[j]<=(y+zoom)]
               obsnew_Rcparts = [obs_Rcparts[j] for j in range(0,len(x_Rcparts)) if x_Rcparts[j]>=(x-zoom) and x_Rcparts[j]<=(x+zoom) and y_Rcparts[j]>=(y-zoom) and y_Rcparts[j]<=(y+zoom)]
               for j in range(0,np.size(xnew_Rcparts)):
                  a = Circle(xy=(xnew_Rcparts[j], ynew_Rcparts[j]), radius=8)
                  a.set_linewidth(1.0)
                  a.set_facecolor('none')
                  plt.text(xnew_Rcparts[j],ynew_Rcparts[j],obsnew_Rcparts[j],horizontalalignment='center',verticalalignment='center',color='limegreen',fontsize=1.5*subsize,alpha=1.0)
                  a.set_edgecolor('limegreen')
                  ax2.add_artist(a)
         if os.path.exists(path+'/2sxps.cat'):
            if np.size(x_Xcparts_Swift)!=1:
               xnew_Xcparts = [x_Xcparts_Swift[j] for j in range(0,len(x_Xcparts_Swift)) if x_Xcparts_Swift[j]>=(x-zoom) and x_Xcparts_Swift[j]<=(x+zoom) and y_Xcparts_Swift[j]>=(y-zoom) and y_Xcparts_Swift[j]<=(y+zoom)]
               ynew_Xcparts = [y_Xcparts_Swift[j] for j in range(0,len(x_Xcparts_Swift)) if x_Xcparts_Swift[j]>=(x-zoom) and x_Xcparts_Swift[j]<=(x+zoom) and y_Xcparts_Swift[j]>=(y-zoom) and y_Xcparts_Swift[j]<=(y+zoom)]
               errnew_Xcparts = [err_Xcparts_Swift[j] for j in range(0,len(x_Xcparts_Swift)) if x_Xcparts_Swift[j]>=(x-zoom) and x_Xcparts_Swift[j]<=(x+zoom) and y_Xcparts_Swift[j]>=(y-zoom) and y_Xcparts_Swift[j]<=(y+zoom)]
               for j in range(0,np.size(xnew_Xcparts)):
                  a = Circle(xy=(xnew_Xcparts[j], ynew_Xcparts[j]), radius=errnew_Xcparts[j])
                  a.set_linewidth(1.0)
                  a.set_facecolor('none')
                  plt.text(xnew_Xcparts[j],ynew_Xcparts[j],'swift',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.5*subsize,alpha=1.0)
                  a.set_edgecolor('fuchsia')
                  ax2.add_artist(a)
            else:
               if x_Xcparts_Swift>=(x-zoom) and x_Xcparts_Swift<=(x+zoom) and y_Xcparts_Swift>=(y-zoom) and y_Xcparts_Swift<=(y+zoom):
                  xnew_Xcparts = x_Xcparts_Swift
                  ynew_Xcparts = y_Xcparts_Swift
                  errnew_Xcparts = err_Xcparts_Swift
                  a = Circle(xy=(xnew_Xcparts, ynew_Xcparts), radius=errnew_Xcparts)
                  a.set_linewidth(1.0)
                  a.set_facecolor('none')
                  plt.text(xnew_Xcparts,ynew_Xcparts,'swift',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.5*subsize,alpha=1.0)
                  a.set_edgecolor('fuchsia')
                  ax2.add_artist(a)
         if os.path.exists(path+'/csc2.cat'):
            if np.size(x_Xcparts_Chandra)!=1:
               xnew_Xcparts = [x_Xcparts_Chandra[j] for j in range(0,len(x_Xcparts_Chandra)) if x_Xcparts_Chandra[j]>=(x-zoom) and x_Xcparts_Chandra[j]<=(x+zoom) and y_Xcparts_Chandra[j]>=(y-zoom) and y_Xcparts_Chandra[j]<=(y+zoom)]
               ynew_Xcparts = [y_Xcparts_Chandra[j] for j in range(0,len(x_Xcparts_Chandra)) if x_Xcparts_Chandra[j]>=(x-zoom) and x_Xcparts_Chandra[j]<=(x+zoom) and y_Xcparts_Chandra[j]>=(y-zoom) and y_Xcparts_Chandra[j]<=(y+zoom)]
               ellanew_Xcparts = [ell_a_Chandra[j] for j in range(0,len(x_Xcparts_Chandra)) if x_Xcparts_Chandra[j]>=(x-zoom) and x_Xcparts_Chandra[j]<=(x+zoom) and y_Xcparts_Chandra[j]>=(y-zoom) and y_Xcparts_Chandra[j]<=(y+zoom)]
               ellbnew_Xcparts = [ell_b_Chandra[j] for j in range(0,len(x_Xcparts_Chandra)) if x_Xcparts_Chandra[j]>=(x-zoom) and x_Xcparts_Chandra[j]<=(x+zoom) and y_Xcparts_Chandra[j]>=(y-zoom) and y_Xcparts_Chandra[j]<=(y+zoom)]
               ellthetanew_Xcparts = [ell_theta_Chandra[j] for j in range(0,len(x_Xcparts_Chandra)) if x_Xcparts_Chandra[j]>=(x-zoom) and x_Xcparts_Chandra[j]<=(x+zoom) and y_Xcparts_Chandra[j]>=(y-zoom) and y_Xcparts_Chandra[j]<=(y+zoom)]
               for j in range(0,np.size(xnew_Xcparts)):
                  a = Ellipse(xy=(xnew_Xcparts[j], ynew_Xcparts[j]), width=2*ellanew_Xcparts[j], height=2*ellbnew_Xcparts[j], angle=ellthetanew_Xcparts[j])
                  a.set_linewidth(1.0)
                  a.set_facecolor('none')
                  plt.text(xnew_Xcparts[j],ynew_Xcparts[j],'chandra',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.5*subsize,alpha=1.0)
                  a.set_edgecolor('fuchsia')
                  ax2.add_artist(a)
            else:
               if x_Xcparts_Chandra>=(x-zoom) and x_Xcparts_Chandra<=(x+zoom) and y_Xcparts_Chandra>=(y-zoom) and y_Xcparts_Chandra<=(y+zoom):
                  xnew_Xcparts = x_Xcparts_Chandra
                  ynew_Xcparts = y_Xcparts_Chandra
                  ellanew_Xcparts = ell_a_Chandra
                  ellbnew_Xcparts = ell_b_Chandra
                  ellthetanew_Xcparts = ell_theta_Chandra
                  a = Ellipse(xy=(xnew_Xcparts, ynew_Xcparts), width=2*ellanew_Xcparts, height=2*ellbnew_Xcparts, angle=ellthetanew_Xcparts)
                  a.set_linewidth(1.0)
                  a.set_facecolor('none')
                  plt.text(xnew_Xcparts,ynew_Xcparts,'chandra',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.5*subsize,alpha=1.0)
                  a.set_edgecolor('fuchsia')
                  ax2.add_artist(a)
         if os.path.exists(path+'/XMM4d13s.cat'):
            if np.size(x_Xcparts_XMM)!=1:
               xnew_Xcparts = [x_Xcparts_XMM[j] for j in range(0,len(x_Xcparts_XMM)) if x_Xcparts_XMM[j]>=(x-zoom) and x_Xcparts_XMM[j]<=(x+zoom) and y_Xcparts_XMM[j]>=(y-zoom) and y_Xcparts_XMM[j]<=(y+zoom)]
               ynew_Xcparts = [y_Xcparts_XMM[j] for j in range(0,len(x_Xcparts_XMM)) if x_Xcparts_XMM[j]>=(x-zoom) and x_Xcparts_XMM[j]<=(x+zoom) and y_Xcparts_XMM[j]>=(y-zoom) and y_Xcparts_XMM[j]<=(y+zoom)]
               errnew_Xcparts = [err_Xcparts_XMM[j] for j in range(0,len(x_Xcparts_XMM)) if x_Xcparts_XMM[j]>=(x-zoom) and x_Xcparts_XMM[j]<=(x+zoom) and y_Xcparts_XMM[j]>=(y-zoom) and y_Xcparts_XMM[j]<=(y+zoom)]
               for j in range(0,np.size(xnew_Xcparts)):
                  a = Circle(xy=(xnew_Xcparts[j], ynew_Xcparts[j]), radius=errnew_Xcparts[j])
                  a.set_linewidth(1.0)
                  a.set_facecolor('none')
                  plt.text(xnew_Xcparts[j],ynew_Xcparts[j],'XMM',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.5*subsize,alpha=1.0)
                  a.set_edgecolor('fuchsia')
                  ax2.add_artist(a)
            else:
               if x_Xcparts_XMM>=(x-zoom) and x_Xcparts_XMM<=(x+zoom) and y_Xcparts_XMM>=(y-zoom) and y_Xcparts_XMM<=(y+zoom):
                  xnew_Xcparts = x_Xcparts_XMM
                  ynew_Xcparts = y_Xcparts_XMM
                  errnew_Xcparts = err_Xcparts_XMM
                  a = Circle(xy=(xnew_Xcparts, ynew_Xcparts), radius=errnew_Xcparts)
                  a.set_linewidth(1.0)
                  a.set_facecolor('none')
                  plt.text(xnew_Xcparts,ynew_Xcparts,'XMM',horizontalalignment='center',verticalalignment='center',color='fuchsia',fontsize=1.5*subsize,alpha=1.0)
                  a.set_edgecolor('fuchsia')
                  ax2.add_artist(a)

         if 'ra_query' in locals() and 'dec_query' in locals():
            plt.text(x_query,y_query,"query",horizontalalignment='center',verticalalignment='center',color='darkorange',fontsize=1.5*subsize,alpha=0.8)
            a = Circle(xy=(x_query, y_query), radius=8)
            a.set_linewidth(0.5)
            a.set_facecolor('none')
            a.set_edgecolor('darkorange')
            ax2.add_artist(a)
   
         if 'kic_ID' in locals():
            plt.text(-1.85,6.0,"Identified from Kepler as "+kic_ID+r" with $P_{\mathrm{Kepler}}$="+str(kic_period)+" d",transform=ax2.transAxes,fontsize=tsize)
         elif 'ATO_ID' in locals():
            plt.text(-1.85,6.0,"Identified from ATLAS as "+ATO_ID+r" with $P_{\mathrm{ATLAS}}$="+str(ATO_period)+" d",transform=ax2.transAxes,fontsize=tsize)
         elif 'SIMBAD_ID' in locals():
            plt.text(-1.7,6.0,"Identified from SIMBAD as "+SIMBAD_ID,transform=ax2.transAxes,fontsize=tsize)
         plt.savefig(path+'/finalperiodplots'+str(nlist)+'_'+per_search+'/finalperplot_'+str(i+1)+'_'+per_search+'_ATLAS_nper'+str(int(n_periods))+'.png',bbox_inches='tight')
         plt.close(figure)

      else:
         print("No ATLAS data acquired for the candidate #"+str(i+1)+".")

   if 'flux_Xcpart' in locals():
      del flux_Xcpart
   if 'flux_Xcpart_Swift' in locals():
      del flux_Xcpart_Swift
   if 'flux_Xcpart_Chandra' in locals():
      del flux_Xcpart_Chandra
   if 'flux_Xcpart_XMM' in locals():
      del flux_Xcpart_XMM
   if 'UL_X' in locals():
      del UL_X
   if 'ATO_ID' in locals():
      del ATO_ID
   if 'ATO_period' in locals():
      del ATO_period
   if 'kic_ID' in locals():
      del kic_ID
   if 'kic_period' in locals():
      del kic_period
   if 'SIMBAD_ID' in locals():
      del SIMBAD_ID
   if 'ra_query' in locals():
      del ra_query
   if 'dec_query' in locals():
      del dec_query
         
if nvar!=1:
   #Follow Python index enumeration
   perSources_fin = np.concatenate((perSources_fin, np.array([list_targ-1]).T), axis=1)
   perSources_fin = np.concatenate((perSources_fin, np.array([best_periods]).T), axis=1)
   np.savetxt(path+'/perSources_'+str(nlist)+'.csv', perSources_fin, delimiter=",", fmt='%0.8f')

toc = time.perf_counter()
print(f"Final periodic plots produced in {toc - tic:0.4f} seconds")
