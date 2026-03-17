#!/bin/python

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rcParams
from matplotlib.patches import Circle
from matplotlib.patches import Ellipse
from astropy.io import fits
from astropy import wcs
from astropy.coordinates import SkyCoord
from astropy.stats import biweight_midcorrelation
import os
import shutil
import time
from astropy.timeseries import LombScargle

from periodicity_search import LombScargle_search, LombScargle_search_test, LombScargle_FAP_level, PDM_search, PDM_FAP_level
   
analysis = input("Enter 'STELLA', 'INT' or 'LCO' for the type of data (default value 'LCO'):\n") or 'LCO'
while analysis!='STELLA' and analysis!='INT' and analysis!='LCO':
   analysis = input("Invalid value, please insert 'STELLA', 'INT' or 'LCO':\n") or 'LCO'
fext = float(input("Enter the scale factor for the FWHM used for the extraction radius (default value 1.2):\n") or 1.2)
while fext<=0.:
   fext = float(input("Invalid value, please insert a positive extraction factor:\n") or 1.2)
nlist = int(input("Enter the list of candidate variables you want to inspect (1 or 2) (default value 1):\n") or 1)
while nlist!=1 and nlist!=2:
   nlist = int(input("Invalid value, please insert 1 or 2:\n") or 1)
P_low = float(input("Enter a positive number for periodlower (in days) to be used in the final Astrosource run (default value 0.02):\n") or "0.02")
while P_low<=0.:
   P_low = float(input("Invalid value, please insert a positive number:\n") or "0.02")

P_up = float(input("Enter a positive number for periodupper (in days) to be used in the final Astrosource run (default value 2.5):\n") or "2.5")
while P_up<=0.:
   P_up = float(input("Invalid value, please insert a positive number:\n") or "2.5")
   
n_periods = float(input("Enter an integer positive number for periodtests to be used in the final Astrosource run (default value 2000):\n") or "2000")
while n_periods<0. or n_periods%1!=0.:
   n_periods = float(input("Invalid value, please insert a positive integer number:\n") or "2000")
n_periods = int(n_periods)
period_FAPl = input("Enter 'yes' if you want to estimate the False-Alarm-Probability level at 0.1% (only for LS method), better than 3sigma (default value 'no', it can increase the computation times):\n")  or "no"
while period_FAPl!='yes' and period_FAPl!='no':
   period_FAPl = input("Invalid value, please insert 'yes' or 'no':\n") or "no"
if period_FAPl == 'yes':
   #Probability of periodic signal in the data
   conf_level = 99.9
   FAPl_factor = float(input("Enter an integer positive number for the reduction factor of the period grid (default value 6):\n")  or "6")
   while FAPl_factor<1 or FAPl_factor>100:
      FAPl_factor = float(input("Invalid value, please insert a integer number between 1 and 100:\n") or "6")
   FAPl_factor = int(FAPl_factor)
g_an = input("Do you want to produce plots also in the g filter? (default value 'no'):\n") or 'no'
while g_an!='yes' and g_an!='no':
   g_an = input("Invalid value, please insert 'yes' or 'no':\n") or 'no'
i_an = input("Do you want to produce plots also in the i filter? (default value 'no'):\n") or 'no'
while i_an!='yes' and i_an!='no':
   i_an = input("Invalid value, please insert 'yes' or 'no':\n") or 'no'

#Defining scale factor (in sigma) to use as squared region around the target, with side=2*bound*semimaj_ax
bound = 2
#Defining radius size (in arcmin) to use for plotting 3FGL sources nearby to our candidate target
radius = 11
#Defining scale factor (in pixel) to use as squared zoom region around the variable source, with side=2*zoom
zoom = 100
#Color range factor to be multiplied to data stdev in pyplot.imshow
colscale = 0.15
#Number of decimal digits to which round the coordinates for not plotting twice the aperture radius of the variable candidate in the last panel (zoomed FoV)
ndig=5
path=os.getcwd()
markersize=8
markersize_main=10
alpha=0.25
#Filter on biweight mid-correlation (robust correlation coefficient) for identifying fake variable sources and exclude them (we perform this check only in the r band)
rthresh1 = 0.8
#Filter on biweight mid-correlation (robust correlation coefficient) for identifying fake variable sources and exclude them if they also have more (or equal) than half of the images with SEP flags (we perform this check only in the r band)
rthresh2 = 0.5
tsize = 35
subsize = 22
axlabelsize = 18
axticksize = 16
#figxsize1 = 25
#figysize1 = 35
#figxsize2 = 30
#figysize2 = 42
figxsize = 35
figysize = 50
#figxsize = 40
#figysize = 65

#We get as reference epoch for period folding the time (in MJD) of the first image in the gband
t_0 = float(np.genfromtxt(path+'/gfilter/usedImages.txt',dtype='str')[0].split('_')[2])

tic = time.perf_counter()
sourcesr = np.genfromtxt(path+'/rfilter/starVariability.csv', delimiter=',')
varSources_fin1=np.genfromtxt(path+'/rfilter/varSources_1.csv', delimiter=',')
#varSources_fin2=np.genfromtxt(path+'/rfilter/varSources_2.csv', delimiter=',')

#Listing the optical filters
filters = ['g', 'r', 'i']

if nlist==1:
   nvar=len(varSources_fin1)
   #nvar=1
   for f in filters:
      if(os.path.isdir(path+'/'+f+'filter/finalplots1')):
         shutil.rmtree(path+'/'+f+'filter/finalplots1')
      os.system('mkdir '+path+'/'+f+'filter/finalplots1')
      if(os.path.isdir(path+'/'+f+'filter/finalcheckplots1')):
         shutil.rmtree(path+'/'+f+'filter/finalcheckplots1')
      os.system('mkdir '+path+'/'+f+'filter/finalcheckplots1')
else:
   nvar=len(varSources_fin2)
   #nvar=1
   for f in filters:
      if(os.path.isdir(path+'/'+f+'filter/finalplots2')):
         shutil.rmtree(path+'/'+f+'filter/finalplots2')
      os.system('mkdir '+path+'/'+f+'filter/finalplots2')
      if(os.path.isdir(path+'/'+f+'filter/finalcheckplots2')):
         shutil.rmtree(path+'/'+f+'filter/finalcheckplots2')
      os.system('mkdir '+path+'/'+f+'filter/finalcheckplots2')

hdul = fits.open(path+'/median_r_after.fits')
hdr = hdul[0].header
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
   hdul = fits.open(path+'/rfilter/'+np.genfromtxt(path+'/rfilter.txt',dtype='str')[0])
   if hdul[0].header['INSTRUME']=='kb84':
      telescope = 'LCO/SBIG'
   else:
      telescope = 'LCO/Sinistro'
   hdul.close()

#Load the corresponding Fermi coordinates, 2sigma ellipses (in degrees) and 1 sigma semi-major axis (in degrees) of the target
hdul = fits.open('/export/work/marcotu/gll_psc_v28.fit')
[[ratar, dectar, ell_a_68, ell_a_95, ell_b_95, ell_theta_95]] = [[hdul[1].data['RAJ2000'][i],hdul[1].data['DEJ2000'][i],hdul[1].data['Conf_68_SemiMajor'][i],3600./float(platescale)*hdul[1].data['Conf_95_SemiMajor'][i],3600./float(platescale)*hdul[1].data['Conf_95_SemiMinor'][i],hdul[1].data['Conf_95_PosAng'][i]] for i in range(0,len(hdul[1].data)) if FoVname.replace("L","L ")==hdul[1].data['Source_Name'][i]]
hdul.close()

#Load the 3FGL Fermi updated coordinates for our target and eventual closeby sources (within radius arcmin from our target), also extended sources
hdul = fits.open('/export/work/marcotu/gll_psc_v16.fit')
sources_3FGL = [hdul[1].data['Source_Name'][i] for i in range(0,len(hdul[1].data))]
coord_3FGL = [[hdul[1].data['RAJ2000'][i],hdul[1].data['DEJ2000'][i],hdul[1].data['Conf_68_SemiMajor'][i],hdul[1].data['Conf_95_SemiMajor'][i],hdul[1].data['Conf_95_SemiMinor'][i],hdul[1].data['Conf_95_PosAng'][i]] for i in range(0,len(hdul[1].data))]
sources_3FGL = np.asarray(sources_3FGL)
coord_3FGL = np.asarray(coord_3FGL)

sep = SkyCoord(ratar, dectar, frame='fk5', unit='deg').separation(SkyCoord(coord_3FGL[:,0],coord_3FGL[:,1],frame='fk5',unit='deg')).arcminute
sources_3FGL = [sources_3FGL[i] for i in range(0,len(sources_3FGL)) if sep[i]<=radius]
sources_3FGL = np.asarray(sources_3FGL)
coord_3FGL = [coord_3FGL[i,:] for i in range(0,len(coord_3FGL)) if sep[i]<=radius]
coord_3FGL = np.asarray(coord_3FGL)
sep = [sep[i] for i in range(0,len(sep)) if sep[i]<=radius]

#Load eventual 4FGL extended closeby sources (within radius arcmin from our target)
#sources_4FGLe = [hdul[2].data['Source_Name'][i] for i in range(0,len(hdul[2].data))]
#coord_4FGLe = [[hdul[2].data['RAJ2000'][i],hdul[2].data['DEJ2000'][i],hdul[2].data['Model_SemiMajor'][i],hdul[2].data['Model_SemiMinor'][i],hdul[2].data['Model_PosAng'][i]] for i in range(0,len(hdul[2].data))]
#sources_4FGLe = np.asarray(sources_4FGLe)
#coord_4FGLe = np.asarray(coord_4FGLe)
#sep = SkyCoord(ratar, dectar, frame='fk5', unit='deg').separation(SkyCoord(coord_4FGLe[:,0],coord_4FGLe[:,1],frame='fk5',unit='deg')).arcminute
#sources_4FGLe = [sources_4FGLe[i] for i in range(0,len(sources_4FGLe)) if sep[i]<=radius]
#coord_4FGLe = [coord_4FGLe[i] for i in range(0,len(coord_4FGLe)) if sep[i]<=radius]
#sep = [sep[i] for i in range(0,len(sep)) if sep[i]<=radius]
#sources_4FGLe = np.asarray(sources_4FGLe)
#coord_4FGLe = np.asarray(coord_4FGLe)
hdul.close()

#We get as reference epoch for period folding the time (in MJD) of the first image in the gband
t_0 = float(np.genfromtxt(path+'/gfilter/usedImages.txt',dtype='str')[0].split('_')[2])

#For loop in the three optical filters
for f in filters:
   if (f=='g' and g_an=='yes') or (f=='i' and i_an=='yes') or f=='r':
      hdul = fits.open(path+'/median_'+f+'_after.fits')
      data= hdul[0].data
      hdr = hdul[0].header
      xsources = hdul[1].data['x']
      ysources = hdul[1].data['y']
      w = wcs.WCS(hdr)
      xtar, ytar = w.all_world2pix(ratar, dectar, 1)
      if np.size(sources_3FGL)!=0:
         x_3FGL, y_3FGL = w.all_world2pix(coord_3FGL[:,0], coord_3FGL[:,1], 1)
      #if np.size(sources_4FGLe)!=0:
      #   x_4FGLe, y_4FGLe = w.all_world2pix(np.asarray(sources_4FGLe)[:,1], np.asarray(sources_4FGLe)[:,2], 1)
      data_crop = data[max(0,round(xtar-bound*max(ell_a_95,ell_b_95))):min(np.shape(data)[1]-1,round(xtar+bound*max(ell_a_95,ell_b_95))),max(0,round(ytar-bound*max(ell_a_95,ell_b_95))):min(np.shape(data)[0]-1,round(ytar+bound*max(ell_a_95,ell_b_95)))]
      hdul.close()

      #Checking nights dates
      images = np.genfromtxt(path+'/'+f+'filter.txt',dtype='str')
      if analysis=='STELLA':
         dates = [s.split("-")[0].split("science")[1].split("A")[0] for s in images]
      elif analysis=='INT':
         dates = []
         for i in range(0,len(images)):
            hdul = fits.open(path+'/'+f+'filter/'+images[i])
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

      #if N_nights==1:
      #   nights_print = nights[0]
      #else:
      #   nights_print = nights[0]+", "
      #   for i in range(1,N_nights):
      #      if i!=N_nights-1:
      #         nights_print = nights_print+nights[i]+", "
      #      else:
      #         nights_print = nights_print+nights[i]

      #Checking nights mjd_ref
      images = np.genfromtxt(path+'/'+f+'filter/usedImages.txt',dtype='str')
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

      m, s = np.mean(data_crop), np.std(data_crop)

      #We get as reference epoch for the light curves the integer part of the time (in MJD) of the first image in the corresponding filter
      ref_ep = np.trunc(float(np.genfromtxt(path+'/'+f+'filter/usedImages.txt',dtype='str')[0].split('_')[2]))
               
      #For loop in the candidate variables
      for i in range(nvar):
         if nlist==1:
            sigma=varSources_fin1[i,3]
            dm=varSources_fin1[i,2]
            ra=varSources_fin1[i,0]
            dec=varSources_fin1[i,1]
         else:
            sigma=varSources_fin2[i,3]
            dm=varSources_fin2[i,2]
            ra=varSources_fin2[i,0]
            dec=varSources_fin2[i,1]

         sep_4FGL = SkyCoord(ra, dec, frame='fk5', unit='deg').separation(SkyCoord(ratar,dectar,frame='fk5',unit='deg')).degree
         if np.size(sources_3FGL)!=0:
            sep_3FGL = SkyCoord(ra, dec, frame='fk5', unit='deg').separation(SkyCoord(coord_3FGL[:,0],coord_3FGL[:,1],frame='fk5',unit='deg')).degree
            ell_a_68_3FGL = coord_3FGL[np.argmin(sep_3FGL),2]
            sep_3FGL = sep_3FGL[np.argmin(sep_3FGL)]

         if os.path.exists(path+'/'+f+'filter/outputcats/doerPhot_V'+str(i+1)+'.csv') and np.shape(np.genfromtxt(path+'/'+f+'filter/outputcats/doerPhot_V'+str(i+1)+'.csv', delimiter=','))[0]!=np.size(np.genfromtxt(path+'/'+f+'filter/outputcats/doerPhot_V'+str(i+1)+'.csv', delimiter=',')) and np.shape(np.genfromtxt(path+'/'+f+'filter/outputcats/doerPhot_V'+str(i+1)+'.csv', delimiter=','))[0]>2:

            xpos, ypos, diffm, errdiffm, fwhm, flag, mjd, airmass, mag, sigma_mag = np.genfromtxt(path+'/'+f+'filter/outputcats/doerPhot_V'+str(i+1)+'.csv', delimiter=',', usecols=(2,3,12,13,6,7,8,9,-2,-1),unpack=True)

            fwhmavg = np.mean(fwhm)

            nightdim = np.zeros(N_nights,int)
            for j in range(0,(np.size(mjd))):
               if j!=0:
                  if (mjd[j]-mjd[j-1]) > 0.8:
                     n_night = np.argmin(np.fabs(mjd[j-1]-mjd_ref))
                     nightdim[n_night] = j-nightdim[0]
                     for k in range(1,n_night):
                        nightdim[n_night] = nightdim[n_night]-nightdim[k]
                  if j==np.size(mjd)-1:
                     n_night = np.argmin(np.fabs(mjd[j]-mjd_ref))
                     nightdim[n_night] = j+1-nightdim[0]
                     for k in range(1,n_night):
                        nightdim[n_night] = nightdim[n_night]-nightdim[k]      

            xpos_single = []
            ypos_single = []
            diffm_single = []
            errdiffm_single = []
            fwhm_single = []
            flag_single = []
            mjd_single = []
            airmass_single = []
            mag_single = []
            sigmamag_single = []
            n_lc = 0
            for j in range(0,N_nights):
               xpos_single.append(xpos[np.sum(nightdim[0:j]):np.sum(nightdim[0:j])+nightdim[j]])
               ypos_single.append(ypos[np.sum(nightdim[0:j]):np.sum(nightdim[0:j])+nightdim[j]])
               diffm_single.append(diffm[np.sum(nightdim[0:j]):np.sum(nightdim[0:j])+nightdim[j]])
               errdiffm_single.append(errdiffm[np.sum(nightdim[0:j]):np.sum(nightdim[0:j])+nightdim[j]])
               fwhm_single.append(fwhm[np.sum(nightdim[0:j]):np.sum(nightdim[0:j])+nightdim[j]])
               flag_single.append(flag[np.sum(nightdim[0:j]):np.sum(nightdim[0:j])+nightdim[j]])
               mjd_single.append(mjd[np.sum(nightdim[0:j]):np.sum(nightdim[0:j])+nightdim[j]])
               mag_single.append(mag[np.sum(nightdim[0:j]):np.sum(nightdim[0:j])+nightdim[j]])
               sigmamag_single.append(sigma_mag[np.sum(nightdim[0:j]):np.sum(nightdim[0:j])+nightdim[j]])
               airmass_single.append(airmass[np.sum(nightdim[0:j]):np.sum(nightdim[0:j])+nightdim[j]])
               if nightdim[j]!=0:
                  n_lc = n_lc+1
               
            rairmass = np.zeros(N_nights)
            rfwhm = np.zeros(N_nights)
            rxpos = np.zeros(N_nights)
            rypos = np.zeros(N_nights)
            checkflags = []
            for j in range(0,N_nights):
               counter = 0
               for k in range(0,len(flag_single[j])):
                  if flag_single[j][k]!=0:
                     counter+=1
               if nightdim[j]!=0:
                  rairmass[j] = biweight_midcorrelation(airmass_single[j], diffm_single[j])
                  rfwhm[j] = biweight_midcorrelation(fwhm_single[j], diffm_single[j])
                  rxpos[j] = biweight_midcorrelation(xpos_single[j], diffm_single[j])
                  rypos[j] = biweight_midcorrelation(ypos_single[j], diffm_single[j])
                  if (np.fabs(rairmass[j]) > rthresh1 or np.fabs(rfwhm[j]) > rthresh1 or np.fabs(rxpos[j]) > rthresh1 or np.fabs(rypos[j]) > rthresh1) or (counter >= len(flag_single[j])/2 and (np.fabs(rairmass[j]) > rthresh2 or np.fabs(rfwhm[j]) > rthresh2 or np.fabs(rxpos[j]) > rthresh2 or np.fabs(rypos[j]) > rthresh2)):
                     checkflags.append(str(j+1))
            if np.size(checkflags)!=0:
               if np.size(checkflags)==1:
                  checkflags_print = checkflags[0]
               else:
                  checkflags_print = checkflags[0]+", "
               for j in range(1,np.size(checkflags)):
                  if j!=np.size(checkflags)-1:
                     checkflags_print = checkflags_print+checkflags[j]+", "
                  else:
                     checkflags_print = checkflags_print+checkflags[j]

            LS1_periods, LS1_periodogram = LombScargle_search(mjd, mag, sigma_mag, P_low, P_up, n_periods, 1)
            best_LS1 = LS1_periods[~np.isnan(LS1_periodogram)][np.argmax(LS1_periodogram[~np.isnan(LS1_periodogram)])]
            phase_LS1 = ((mjd-t_0)/best_LS1) % 1

            PDM_periods, PDM_periodogram = PDM_search(mjd, mag, P_low, P_up, 0., n_periods, int(np.trunc(np.size(mjd)/3)), 2)
            best_PDM = PDM_periods[~np.isnan(PDM_periodogram)][np.argmin(PDM_periodogram[~np.isnan(PDM_periodogram)])]
            phase_PDM = ((mjd-t_0)/best_PDM) % 1

            #Calculate the Lomb-Scargle and PDM periodogram values for a function gaussian distributed, with mu=mean of our data and sigma=0.01 mag, to compare the peaks we find in the periodogram of our time series
            mag_test = np.random.normal(np.mean(mag), 0.01, np.size(mjd))
            LS1_periods_test, LS1_periodogram_test = LombScargle_search_test(mjd, mag_test, P_low, P_up, n_periods, 1)
            PDM_periods_test, PDM_periodogram_test = PDM_search(mjd, mag_test, P_low, P_up, 0., n_periods, int(np.trunc(np.size(mjd)/3)), 2)

            figure=plt.figure(figsize=(figxsize,figysize), dpi=150)
            if np.size(checkflags)!=0:
               if np.size(sources_3FGL)!=0:
                  figure.suptitle(r"Fermi source: "+FoVname+", Instrument: "+telescope+", Filter: "+f+", Source detected in "+str(n_lc)+" nights out of "+str(N_nights)+"\n\nResolution time: %.1f"%deltat+" s, Plate scale: "+str(platescale)+" arcsec/pix, Extraction radius: "+str(fext)+"*FWHM, Checkflags on nights: "+checkflags_print+"\n\nList #%d"%nlist+", Candidate #%d"%(i+1)+", $\sigma$: %.2f"%sigma+" mag, $\Delta$m: %.2f"%dm+" mag, Range: [%.2f"%P_low+",%.1f"%min(P_up,mjd[-1]-mjd[0])+"] d, Steps: %d"%n_periods+", Period folding with $T_{0}=$%.4f"%t_0+" MJD\n\n%.1f"%(sep_4FGL/ell_a_68)+r"$\sigma$ far from the 4FGL source and %.1f"%(sep_3FGL/ell_a_68_3FGL)+r"$\sigma$ far from the 3FGL source", fontsize=tsize)
               else:
                  figure.suptitle(r"Fermi source: "+FoVname+", Instrument: "+telescope+", Filter: "+f+", Source detected in "+str(n_lc)+" nights out of "+str(N_nights)+"\n\nResolution time: %.1f"%deltat+" s, Plate scale: "+str(platescale)+" arcsec/pix, Extraction radius: "+str(fext)+"*FWHM, Checkflags on nights: "+checkflags_print+"\n\nList #%d"%nlist+", Candidate #%d"%(i+1)+", $\sigma$: %.2f"%sigma+" mag, $\Delta$m: %.2f"%dm+" mag, Range: [%.2f"%P_low+",%.1f"%min(P_up,mjd[-1]-mjd[0])+"] d, Steps: %d"%n_periods+", Period folding with $T_{0}=$%.4f"%t_0+" MJD\n\n%.1f"%(sep_4FGL/ell_a_68)+r"$\sigma$ far from the 4FGL source", fontsize=tsize)
            else:
               if np.size(sources_3FGL)!=0:
                  figure.suptitle(r"Fermi source: "+FoVname+", Instrument: "+telescope+", Filter: "+f+", Source detected in "+str(n_lc)+" nights out of "+str(N_nights)+"\n\nResolution time: %.1f"%deltat+" s, Plate scale: "+str(platescale)+" arcsec/pix, Extraction radius: "+str(fext)+"*FWHM\n\nList #%d"%nlist+", Candidate #%d"%(i+1)+", $\sigma$: %.2f"%sigma+" mag, $\Delta$m: %.2f"%dm+" mag, Range: [%.2f"%P_low+",%.1f"%min(P_up,mjd[-1]-mjd[0])+"] d, Steps: %d"%n_periods+", Period folding with $T_{0}=$%.4f"%t_0+" MJD\n\n%.1f"%(sep_4FGL/ell_a_68)+r"$\sigma$ far from the 4FGL source and %.1f"%(sep_3FGL/ell_a_68_3FGL)+r"$\sigma$ far from the 3FGL source", fontsize=tsize)
               else:
                  figure.suptitle(r"Fermi source: "+FoVname+", Instrument: "+telescope+", Filter: "+f+", Source detected in "+str(n_lc)+" nights out of "+str(N_nights)+"\n\nResolution time: %.1f"%deltat+" s, Plate scale: "+str(platescale)+" arcsec/pix, Extraction radius: "+str(fext)+"*FWHM\n\nList #%d"%nlist+", Candidate #%d"%(i+1)+", $\sigma$: %.2f"%sigma+" mag, $\Delta$m: %.2f"%dm+" mag, Range: [%.2f"%P_low+",%.1f"%min(P_up,mjd[-1]-mjd[0])+"] d, Steps: %d"%n_periods+", Period folding with $T_{0}=$%.4f"%t_0+" MJD\n\n%.1f"%(sep_4FGL/ell_a_68)+r"$\sigma$ far from the 4FGL source", fontsize=tsize)

            counter = 0
            for j in range(0,N_nights):
               if nightdim[j]!=0:
                  plt.subplot(int(np.trunc(5+n_lc/2)), 2, counter+1)
                  plt.plot(mjd_single[j]-ref_ep, mag_single[j], 'bo', linestyle='None')
                  plt.errorbar(mjd_single[j]-ref_ep, mag_single[j], sigmamag_single[j], linestyle='None')
                  plt.xticks(fontsize=axticksize)
                  plt.yticks(fontsize=axticksize)
                  plt.gca().invert_yaxis()
                  plt.title("Lightcurve of "+nights[j],fontsize=subsize,fontweight='bold')
                  plt.xlabel(r"Time (MJD-%d)"%ref_ep,fontsize=axlabelsize)
                  plt.ylabel(r"Apparent "+f+" magnitude",fontsize=axlabelsize)
                  counter = counter+1

            plt.subplot(int(np.trunc(5+n_lc/2)), 2, n_lc+1)
            plt.xlabel(r'$\Delta$m (mag in the r band)',fontsize=axlabelsize)
            plt.ylabel(r'$\sigma$ (mag)',fontsize=axlabelsize)
            plt.title("Dispersion vs differential magnitude",fontsize=subsize,fontweight='bold')
            plt.yscale("log")
            plt.plot(sourcesr[:,2],sourcesr[:,3],'ob',markerfacecolor='0.8',markeredgecolor='black',markersize=markersize,alpha=alpha)
            plt.plot(varSources_fin1[:,2],varSources_fin1[:,3],'ob',markerfacecolor='red',markeredgecolor='black',markersize=markersize,alpha=alpha)
            #plt.plot(varSources_fin2[:,2],varSources_fin2[:,3],'ob',markerfacecolor='green',markeredgecolor='black',markersize=markersize,alpha=alpha)
            if nlist==1:
               plt.plot(varSources_fin1[i,2],varSources_fin1[i,3],'ob',markerfacecolor='red',markeredgecolor='black',markersize=markersize_main)
            else:
               plt.plot(varSources_fin2[i,2],varSources_fin2[i,3],'ob',markerfacecolor='green',markeredgecolor='black',markersize=markersize_main)
            plt.xticks(fontsize=axticksize)
            plt.yticks(fontsize=axticksize)

            #ax3 = plt.subplot(7, 2, 3)
            #plt.plot(LS1_lc[:,0]-np.trunc(LS1_lc[0,0]), LS1_lc[:,2], 'bo', linestyle='None')
            #plt.xticks(fontsize=axticksize)
            #plt.yticks(fontsize=axticksize)
            #plt.gca().invert_yaxis()
            #plt.title("Gaussian test lightcurve with mean "+f+"={:.2f} mag".format(np.mean(LS1_lc[:,3])),fontsize=subsize)
            #plt.xlabel(r"Time (MJD-%d)"%np.trunc(LS1_lc[0,0]),fontsize=axlabelsize)
            #plt.ylabel(r"Apparent "+f+" magnitude",fontsize=axlabelsize)
            #plt.ylim(max(LS1_lc[:,3]),min(LS1_lc[:,3]))

            plt.subplot(int(np.trunc(5+n_lc/2)), 2, n_lc+2)
            plt.plot(LS1_periods_test,LS1_periodogram_test)
            if not np.isnan(max(LS1_periodogram_test)) and not np.isinf(max(LS1_periodogram_test)) and not np.isnan(best_LS1) and not np.isinf(best_LS1):
               #plt.ylim(0, max(max(LS1_periodogram[~np.isnan(LS1_periodogram)] or np.isinf(LS1_periodogram))]), max(LS1_periodogram_test[~(np.isnan(LS1_periodogram_test) or np.isinf(LS1_periodogram_test))]))+0.1)
               plt.axvline(x=best_LS1, color='r', ls='--')
            plt.xscale("log")
            plt.xticks(fontsize=axticksize)
            plt.yticks(fontsize=axticksize)
            plt.title("N=1 Lomb-Scargle test periodogram on Gaussian window",fontsize=subsize)
            plt.xlabel(r"Trial period (days)",fontsize=axlabelsize)
            plt.ylabel(r"Power",fontsize=axlabelsize)

            plt.subplot(int(np.trunc(5+n_lc/2)), 2, n_lc+3)
            plt.plot(phase_LS1, mag, 'bo', linestyle='None')
            plt.plot(phase_LS1+1, mag, 'ro', linestyle='None')
            plt.errorbar(phase_LS1, mag, sigma_mag, linestyle='None')
            plt.errorbar(phase_LS1+1, mag, sigma_mag, linestyle='None')
            plt.xticks(fontsize=axticksize)
            plt.yticks(fontsize=axticksize)
            plt.gca().invert_yaxis()
            plt.title("N=1 Lomb-Scargle folded lightcurve: P={0} d".format(best_LS1),fontsize=subsize,fontweight='bold')
            plt.xlabel(r"Phase ($\phi$)",fontsize=axlabelsize)
            plt.ylabel(r"Apparent "+f+" magnitude",fontsize=axlabelsize)

            plt.subplot(int(np.trunc(5+n_lc/2)), 2, n_lc+4)
            plt.plot(LS1_periods,LS1_periodogram)
            if not np.isnan(max(LS1_periodogram)) and not np.isinf(max(LS1_periodogram)) and not np.isnan(best_LS1) and not np.isinf(best_LS1):
               #plt.ylim(0, max(max(LS1_periodogram[~np.isnan(np.isnan(LS1_periodogram) or np.isinf(LS1_periodogram))]), max(LS1_periodogram_test[~np.isnan(np.isnan(LS1_periodogram_test) or np.isinf(LS1_periodogram_test))]))+0.1)
               plt.axvline(x=best_LS1, color='r', ls='--')
            plt.xscale("log")
            plt.xticks(fontsize=axticksize)
            plt.yticks(fontsize=axticksize)
            if period_FAPl == 'yes':
               tic_1 = time.perf_counter()
               FAP_periods, FAP_levels = LombScargle_FAP_level(mjd, mag, sigma_mag, P_low, P_up, n_periods, 1, FAPl_factor, conf_level)
               toc_1 = time.perf_counter()
               print(f"FAP levels computed in {toc_1 - tic_1:0.4f} seconds")
               plt.plot(FAP_periods,FAP_levels,color='darkorange',ls='--')
               if max(LS1_periodogram[~np.isnan(LS1_periodogram)]) >= FAP_levels[np.argmin(np.fabs(FAP_periods-LS1_periods[~np.isnan(LS1_periodogram)][np.argmax(LS1_periodogram[~np.isnan(LS1_periodogram)])]))]:
                  plt.title(r"N=1 Lomb-Scargle data periodogram, >0.1% FAP",fontsize=subsize,fontweight='bold')
               else:
                  plt.title(r"N=1 Lomb-Scargle data periodogram, <0.1% FAP",fontsize=subsize,fontweight='bold')
            else:
               plt.title("N=1 Lomb-Scargle data periodogram",fontsize=subsize,fontweight='bold')
            plt.xlabel(r"Trial period (days)",fontsize=axlabelsize)
            plt.ylabel(r"Power",fontsize=axlabelsize)

            plt.subplot(int(np.trunc(5+n_lc/2)), 2, n_lc+5)
            plt.plot(PDM_periods_test,PDM_periodogram_test)
            if not np.isnan(max(PDM_periodogram_test)) and not np.isinf(max(PDM_periodogram_test)) and not np.isnan(best_PDM) and not np.isinf(best_PDM):
               #plt.ylim(0, max(max(PDM_periodogram[~np.isnan(np.isnan(PDM_periodogram) or np.isinf(PDM_periodogram))]), max(PDM_periodogram_test[~np.isnan(np.isnan(PDM_periodogram_test) or np.isinf(PDM_periodogram_test))]))+0.1)
               plt.axvline(x=best_PDM, color='r', ls='--')
            plt.xscale("log")
            plt.xticks(fontsize=axticksize)
            plt.yticks(fontsize=axticksize)
            plt.gca().invert_yaxis()
            plt.title("PDM test periodogram on Gaussian window",fontsize=subsize)
            plt.xlabel(r"Trial period (days)",fontsize=axlabelsize)
            plt.ylabel(r"Likelihood of period",fontsize=axlabelsize)

            plt.subplot(int(np.trunc(5+n_lc/2)), 2, n_lc+6)
            plt.plot(phase_PDM, mag, 'bo', linestyle='None')
            plt.plot(phase_PDM+1, mag, 'ro', linestyle='None')
            plt.errorbar(phase_PDM, mag, sigma_mag, linestyle='None')
            plt.errorbar(phase_PDM+1, mag, sigma_mag, linestyle='None')
            plt.xticks(fontsize=axticksize)
            plt.yticks(fontsize=axticksize)
            plt.gca().invert_yaxis()
            plt.title("PDM folded lightcurve: P={0} d".format(best_PDM),fontsize=subsize,fontweight='bold')
            plt.xlabel(r"Phase ($\phi$)",fontsize=axlabelsize)
            plt.ylabel(r"Apparent "+f+" magnitude",fontsize=axlabelsize)

            plt.subplot(int(np.trunc(5+n_lc/2)), 2, n_lc+7)
            plt.plot(PDM_periods,PDM_periodogram)
            plt.axvline(x=best_PDM, color='r', ls='--')
            plt.xscale("log")
            plt.xticks(fontsize=axticksize)
            plt.yticks(fontsize=axticksize)
            plt.gca().invert_yaxis()
            plt.title("PDM data periodogram",fontsize=subsize,fontweight='bold')
            plt.xlabel(r"Trial period (days)",fontsize=axlabelsize)
            plt.ylabel(r"Likelihood of period",fontsize=axlabelsize)

            ax1 = plt.subplot(int(np.trunc(5+n_lc/2)), 2, n_lc+8)
            plt.title(r'Median '+f+' FoV bkg subtracted in the '+str(2*bound)+'$\sigma$ square around the target coordinates',fontsize=subsize,fontweight='bold')
            plt.xlim(max(0,xtar-bound*max(ell_a_95,ell_b_95)), min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)))
            plt.ylim(max(0,ytar-bound*max(ell_a_95,ell_b_95)), min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)))
            plt.imshow(data, interpolation='nearest', cmap='gray', vmin=m-colscale*s, vmax=m+colscale*s, origin='lower')
            plt.text(xtar,ytar,FoVname,horizontalalignment='center',verticalalignment='center',color='cyan',fontsize=subsize)
            target = Ellipse(xy=(xtar, ytar), width=2*ell_a_95, height=2*ell_b_95, angle=ell_theta_95)
            target.set_facecolor('none')
            target.set_edgecolor('cyan')
            ax1.add_artist(target)
            if np.size(sources_3FGL)!=0:
               for j in range(0,np.shape(sources_3FGL)[0]):
                  if x_3FGL[j]>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_3FGL[j]<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_3FGL[j]>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_3FGL[j]<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
                     plt.text(x_3FGL[j],y_3FGL[j],sources_3FGL[j].replace("L ","L"),horizontalalignment='center',verticalalignment='center',color='yellow',fontsize=subsize)
                     a = Ellipse(xy=(x_3FGL[j], y_3FGL[j]), width=2*3600./float(platescale)*coord_3FGL[j,2], height=2*3600./float(platescale)*coord_3FGL[j,3], angle=coord_3FGL[j,4])
                     a.set_facecolor('none')
                     a.set_edgecolor('yellow')
                     ax1.add_artist(a)
            #if np.size(sources_4FGLe)!=0:
            #   for j in range(0,np.shape(sources_4FGLe)[0]):
            #      if x_4FGLe[j]>=max(0,xtar-bound*max(ell_a_95,ell_b_95)) and x_4FGLe[j]<=min(np.shape(data)[1]-1,xtar+bound*max(ell_a_95,ell_b_95)) and y_4FGLe[j]>=max(0,ytar-bound*max(ell_a_95,ell_b_95)) and y_4FGLe[j]<=min(np.shape(data)[0]-1,ytar+bound*max(ell_a_95,ell_b_95)):
            #         plt.text(x_4FGLe[j],y_4FGLe[j],sources_4FGLe[j][0].replace("L ","L"),horizontalalignment='center',verticalalignment='center',color='cyan',fontsize=subsize)
            #         a = Ellipse(xy=(x_4FGLe[j], y_4FGLe[j]), width=2*3600./float(platescale)*np.asarray(sources_4FGLe)[:,3], height=2*3600./float(platescale)*np.asarray(sources_4FGLe)[:,4], angle=np.asarray(sources_4FGLe)[:,5])
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

            ax2 = plt.subplot(int(np.trunc(5+n_lc/2)), 2, n_lc+9)
            plt.title('Zoom on the variable source: (RA, DEC) = (%.3f'%ra+', %.3f'%dec+') deg',fontsize=subsize,fontweight='bold')
            plt.xlim(x-zoom,x+zoom)
            plt.ylim(y-zoom,y+zoom)
            plt.imshow(data, interpolation='nearest', cmap='gray', vmin=m-colscale*s, vmax=m+colscale*s, origin='lower')
            xsourcesnew = [xsources[j] for j in range(0,len(xsources)) if xsources[j]>=(x-zoom) and xsources[j]<=(x+zoom) and ysources[j]>=(y-zoom) and ysources[j]<=(y+zoom)]
            ysourcesnew = [ysources[j] for j in range(0,len(ysources)) if xsources[j]>=(x-zoom) and xsources[j]<=(x+zoom) and ysources[j]>=(y-zoom) and ysources[j]<=(y+zoom)]
            for j in range(len(xsourcesnew)):
               if not(np.round(x,ndig)==np.round(xsourcesnew[j],ndig) and np.round(y,ndig)==np.round(ysourcesnew[j],ndig)):
                  a = Circle(xy=(xsourcesnew[j], ysourcesnew[j]), radius=fext*fwhmavg)
                  a.set_linewidth(1.0)
                  a.set_facecolor('none')
                  a.set_edgecolor('blue')
                  #a.set_alpha(2*alpha)
                  ax2.add_artist(a)
            plt.text(x,y+20,str(i+1),horizontalalignment='center',verticalalignment='center',color='r',fontsize=2*subsize)
            a = Circle(xy=(x, y), radius=fext*fwhmavg)
            a.set_linewidth(1.0)
            a.set_facecolor('none')
            a.set_edgecolor('red')
            ax2.add_artist(a)

            plt.savefig(path+'/'+f+'filter/finalplots'+str(nlist)+'/finalplot_'+str(i+1)+'.png',bbox_inches='tight')
            plt.close(figure)

            for j in range(0,N_nights):
               if nightdim[j]!=0:
                  figure=plt.figure(figsize=(figxsize,figysize), dpi=150)

                  figure.suptitle("List #%d"%nlist+" Candidate #%d"%(i+1)+" Night of "+nights[j]+"\nCheck on SEP flags and dependence of the "+f+" photometry with airmass, seeing (FWHM), x and y positions", fontsize=tsize)

                  ax1 = plt.subplot(6, 1, 1)
                  plt.plot(mjd_single[j]-np.trunc(mjd_single[j][0]), diffm_single[j], 'ko', linestyle='None')
                  plt.errorbar(mjd_single[j]-np.trunc(mjd_single[j][0]), diffm_single[j], errdiffm_single[j], ecolor='k', linestyle='None')
                  plt.xlim(mjd_single[j][0]-np.trunc(mjd_single[j][0])-0.003,mjd_single[j][len(mjd_single[j])-1]-np.trunc(mjd_single[j][0])+0.003)
                  plt.xticks(fontsize=axticksize)
                  plt.yticks(fontsize=axticksize)
                  plt.gca().invert_yaxis()
                  plt.title("Photometry differential magnitudes vs time",fontsize=subsize)
                  plt.xlabel(r"Time (MJD-%d)"%np.trunc(mjd_single[j][0]),fontsize=axlabelsize)
                  plt.ylabel(r"Differential magnitude",fontsize=axlabelsize)

                  flag0 = np.empty(len(mjd_single[j]))
                  flag0[:] = np.nan
                  flag1 = np.empty(len(mjd_single[j]))
                  flag1[:] = np.nan
                  flag2 = np.empty(len(mjd_single[j]))
                  flag2[:] = np.nan
                  for k in range(0,len(flag_single[j])):
                     if int(flag_single[j][k]/4)!=0:
                        flag2[k]=2
                     if flag_single[j][k]%4==2 or flag_single[j][k]%4==3:
                        flag1[k]=1
                     if flag_single[j][k]%2==1:
                        flag0[k]=0

                  ax2 = plt.subplot(6, 1, 2)
                  plt.plot(mjd_single[j]-np.trunc(mjd_single[j][0]), flag0, 'ro', linestyle='None', markersize=markersize_main, label='Flag 0: nearby sources or bad pixels')
                  plt.plot(mjd_single[j]-np.trunc(mjd_single[j][0]), flag1, 'r^', linestyle='None', markersize=markersize_main, label='Flag 1: deblended object')
                  plt.plot(mjd_single[j]-np.trunc(mjd_single[j][0]), flag2, 'rs', linestyle='None', markersize=markersize_main, label='Flag 2: saturated pixels')
                  plt.xlim(mjd_single[j][0]-np.trunc(mjd_single[j][0])-0.003,mjd_single[j][len(mjd_single[j])-1]-np.trunc(mjd_single[j][0])+0.003)
                  plt.xticks(fontsize=axticksize)
                  plt.yticks(ticks=[0,1,2], fontsize=axticksize)
                  plt.title("SEP detection flags vs time",fontsize=subsize)
                  plt.xlabel(r"Time (MJD-%d)"%np.trunc(mjd_single[j][0]),fontsize=axlabelsize)
                  plt.ylabel(r"Number of raised flag",fontsize=axlabelsize)
                  ax2.legend(fontsize=subsize)
      
                  ax3 = plt.subplot(6, 1, 3)
                  if np.fabs(rairmass[j]) > rthresh2:
                     plt.plot(airmass_single[j], diffm_single[j], 'ro', linestyle='None')
                     plt.errorbar(airmass_single[j], diffm_single[j], errdiffm_single[j], ecolor='r', linestyle='None')
                     plt.xticks(fontsize=axticksize)
                     plt.yticks(fontsize=axticksize)
                     plt.gca().invert_yaxis()
                     plt.title("Photometry differential magnitudes vs airmass",fontsize=subsize)
                     plt.xlabel(r"Airmass",fontsize=axlabelsize)
                     plt.ylabel(r"Differential magnitude",fontsize=axlabelsize)
                     plt.text(0.15,0.9,"Correlation coefficient r=%.1f"%rairmass[j],horizontalalignment='center',verticalalignment='center',transform=ax3.transAxes,color='r',fontsize=subsize)
                  else:
                     plt.plot(airmass_single[j], diffm_single[j], 'ko', linestyle='None')
                     plt.errorbar(airmass_single[j], diffm_single[j], errdiffm_single[j], ecolor='k', linestyle='None')
                     plt.xticks(fontsize=axticksize)
                     plt.yticks(fontsize=axticksize)
                     plt.gca().invert_yaxis()
                     plt.title("Photometry differential magnitudes vs airmass",fontsize=subsize)
                     plt.xlabel(r"Airmass",fontsize=axlabelsize)
                     plt.ylabel(r"Differential magnitude",fontsize=axlabelsize)
                     plt.text(0.15,0.9,"Correlation coefficient r=%.1f"%rairmass[j],horizontalalignment='center',verticalalignment='center',transform=ax3.transAxes,color='k',fontsize=subsize)

                  ax4 = plt.subplot(6, 1, 4)
                  if np.fabs(rfwhm[j]) > rthresh2:
                     plt.plot(fwhm_single[j], diffm_single[j], 'ro', linestyle='None')
                     plt.errorbar(fwhm_single[j], diffm_single[j], errdiffm_single[j], ecolor='r', linestyle='None')
                     plt.xticks(fontsize=axticksize)
                     plt.yticks(fontsize=axticksize)
                     plt.gca().invert_yaxis()
                     plt.title("Photometry differential magnitudes vs seeing (FWHM)", fontsize=subsize)
                     plt.xlabel(r"FWHM (pixel)",fontsize=axlabelsize)
                     plt.ylabel(r"Differential magnitude",fontsize=axlabelsize)
                     plt.text(0.15,0.9,"Correlation coefficient r=%.1f"%rfwhm[j],horizontalalignment='center',verticalalignment='center',transform=ax4.transAxes,color='r',fontsize=subsize)
                  else:
                     plt.plot(fwhm_single[j], diffm_single[j], 'ko', linestyle='None')
                     plt.errorbar(fwhm_single[j], diffm_single[j], errdiffm_single[j], ecolor='k', linestyle='None')
                     plt.xticks(fontsize=axticksize)
                     plt.yticks(fontsize=axticksize)
                     plt.gca().invert_yaxis()
                     plt.title("Photometry differential magnitudes vs seeing (FWHM)",fontsize=subsize)
                     plt.xlabel(r"FWHM (pixel)",fontsize=axlabelsize)
                     plt.ylabel(r"Differential magnitude",fontsize=axlabelsize)
                     plt.text(0.15,0.9,"Correlation coefficient r=%.1f"%rfwhm[j],horizontalalignment='center',verticalalignment='center',transform=ax4.transAxes,color='k',fontsize=subsize)

                  ax5 = plt.subplot(6, 1, 5)
                  if np.fabs(rxpos[j]) > rthresh2:
                     plt.plot(xpos_single[j], diffm_single[j], 'ro', linestyle='None')
                     plt.errorbar(xpos_single[j], diffm_single[j], errdiffm_single[j], ecolor='r', linestyle='None')
                     plt.xticks(fontsize=axticksize)
                     plt.yticks(fontsize=axticksize)
                     plt.gca().invert_yaxis()
                     plt.title("Photometry differential magnitudes vs x position", fontsize=subsize)
                     plt.xlabel(r"x (pixel)",fontsize=axlabelsize)
                     plt.ylabel(r"Differential magnitude",fontsize=axlabelsize)
                     plt.text(0.15,0.9,"Correlation coefficient r=%.1f"%rxpos[j],horizontalalignment='center',verticalalignment='center',transform=ax5.transAxes,color='r',fontsize=subsize)
                  else:
                     plt.plot(xpos_single[j], diffm_single[j], 'ko', linestyle='None')
                     plt.errorbar(xpos_single[j], diffm_single[j], errdiffm_single[j], ecolor='r', linestyle='None')
                     plt.xticks(fontsize=axticksize)
                     plt.yticks(fontsize=axticksize)
                     plt.gca().invert_yaxis()
                     plt.title("Photometry differential magnitudes vs x position", fontsize=subsize)
                     plt.xlabel(r"x (pixel)",fontsize=axlabelsize)
                     plt.ylabel(r"Differential magnitude",fontsize=axlabelsize)
                     plt.text(0.15,0.9,"Correlation coefficient r=%.1f"%rxpos[j],horizontalalignment='center',verticalalignment='center',transform=ax5.transAxes,color='k',fontsize=subsize)

                  ax6 = plt.subplot(6, 1, 6)
                  if np.fabs(rypos[j]) > rthresh2:
                     plt.plot(ypos_single[j], diffm_single[j], 'ro', linestyle='None')
                     plt.errorbar(ypos_single[j], diffm_single[j], errdiffm_single[j], ecolor='r', linestyle='None')
                     plt.xticks(fontsize=axticksize)
                     plt.yticks(fontsize=axticksize)
                     plt.gca().invert_yaxis()
                     plt.title("Photometry differential magnitudes vs y position", fontsize=subsize)
                     plt.xlabel(r"y (pixel)",fontsize=axlabelsize)
                     plt.ylabel(r"Differential magnitude",fontsize=axlabelsize)
                     plt.text(0.15,0.9,"Correlation coefficient r=%.1f"%rypos[j],horizontalalignment='center',verticalalignment='center',transform=ax6.transAxes,color='r',fontsize=subsize)
                  else:
                     plt.plot(ypos_single[j], diffm_single[j], 'ko', linestyle='None')
                     plt.errorbar(ypos_single[j], diffm_single[j], errdiffm_single[j], ecolor='r', linestyle='None')
                     plt.xticks(fontsize=axticksize)
                     plt.yticks(fontsize=axticksize)
                     plt.gca().invert_yaxis()
                     plt.title("Photometry differential magnitudes vs y position", fontsize=subsize)
                     plt.xlabel(r"y (pixel)",fontsize=axlabelsize)
                     plt.ylabel(r"Differential magnitude",fontsize=axlabelsize)
                     plt.text(0.15,0.9,"Correlation coefficient r=%.1f"%rypos[j],horizontalalignment='center',verticalalignment='center',transform=ax6.transAxes,color='k',fontsize=subsize)

                  plt.savefig(path+'/'+f+'filter/finalcheckplots'+str(nlist)+'/finalplot_'+str(i+1)+'_check_night'+str(j+1)+'.png',bbox_inches='tight')
                  plt.close(figure)

         else:
            print("No output files for the variable %d"%(i+1)+" in the "+f+" filter, neither Lomb-Scargle nor PDM methods worked.\n")
   
toc = time.perf_counter()
print(f"Final plots produced in {toc - tic:0.4f} seconds")
