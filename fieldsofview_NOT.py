#!/bin/python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.pyplot import cm
#from matplotlib.colors import LogNorm
from pylab import cm
import os
from astropy.io import fits
from astropy import wcs
from astropy.coordinates import SkyCoord
from matplotlib.patches import Ellipse, Circle
#from astropy.utils.data import get_pkg_data_filename
import sep

path = os.getcwd()

#Defining scale factor to use as squared region around the target, with side=2*bound*semi_axis
#bound = 2

#Defining scale factor (in pixel) to use as squared zoom region around the target, with side=2*zoom
zoom = 50

gimages = np.genfromtxt(path+'/gfilter/usedImages.txt',dtype='str')
ng = np.size(gimages)-1
gcomps = np.atleast_2d(np.genfromtxt(path+'/gfilter/compsUsed.csv', delimiter=','))
gsigma = np.mean(gcomps[:,2])
rimages = np.genfromtxt(path+'/rfilter/usedImages.txt',dtype='str')
nr = np.size(rimages)-1
rcomps = np.atleast_2d(np.genfromtxt(path+'/rfilter/compsUsed.csv', delimiter=','))
rsigma = np.mean(rcomps[:,2])
iimages = np.genfromtxt(path+'/ifilter/usedImages.txt',dtype='str')
ni = np.size(iimages)-1
icomps = np.atleast_2d(np.genfromtxt(path+'/ifilter/compsUsed.csv', delimiter=','))
isigma = np.mean(icomps[:,2])

#Checking night dates
files = np.genfromtxt(path+'/rfilter.txt',dtype='str')

dates = []
for i in range(0,len(files)):
   hdul = fits.open(path+'/rfilter/'+files[i])
   hdr = hdul[0].header
   if hdr['DATE-OBS'].split('T')[1][0]=='0':
      date = hdr['DATE-OBS'].split('T')[0][:-1]+str(int(hdr['DATE-OBS'].split('T')[0][-1])-1)
      if date[-1]=='0':
         date = date[:-2]+'31'
   else:
      date = hdr['DATE-OBS'].split('T')[0]
   dates.append(date)
   hdul.close()
            
nights = []
for i in range(0,np.size(dates)):
   if i==0:
      nights.append(dates[i])
   else:
      if dates[i]!=dates[i-1]:
         nights.append(dates[i])
            
N_nights = np.size(nights)                                    

#Checking nights mjd_ref in g, r and i filters
mjd = [float(s.split('_')[2]) for s in gimages if 'median' not in s]
gmjd_ref = np.zeros(N_nights)
counter = 0
for i in range(0,np.size(mjd)):
            if i!=0:
                        if (mjd[i]-mjd[i-1]) > 0.8:
                                    gmjd_ref[counter] = mjd[i-1]
                                    counter+=1
            if i==np.size(mjd)-1:
                        gmjd_ref[counter] = mjd[i]
mjd = [float(s.split('_')[2]) for s in rimages if 'median' not in s]
rmjd_ref = np.zeros(N_nights)
counter = 0
for i in range(0,np.size(mjd)):
            if i!=0:
                        if (mjd[i]-mjd[i-1]) > 0.8:
                                    rmjd_ref[counter] = mjd[i-1]
                                    counter+=1
            if i==np.size(mjd)-1:
                        rmjd_ref[counter] = mjd[i]
mjd = [float(s.split('_')[2]) for s in iimages if 'median' not in s]
imjd_ref = np.zeros(N_nights)
counter = 0
for i in range(0,np.size(mjd)):
            if i!=0:
                        if (mjd[i]-mjd[i-1]) > 0.8:
                                    imjd_ref[counter] = mjd[i-1]
                                    counter+=1
            if i==np.size(mjd)-1:
                        imjd_ref[counter] = mjd[i]

#gfilter field plot

image = fits.open(path+'/gfilter/median_g.fits')
hdr = image[0].header
FoVname = hdr['OBJECT']
pixscale = 0.43
gw = wcs.WCS(hdr)
gdata = image[0].data

#Load the target coordinates from the header of the fits
ra = hdr['OBJRA']*360/24
dec = hdr['OBJDEC']

gxtar, gytar = gw.all_world2pix(ra, dec, 1)

try:
    bkg = sep.Background(gdata, bw=8, bh=8, fw=3, fh=3)
except ValueError:
    gdata = gdata.byteswap(True).newbyteorder()
    bkg = sep.Background(gdata, bw=8, bh=8, fw=3, fh=3)
bkg.subfrom(gdata)
hdu1 = fits.PrimaryHDU(data=gdata,header=hdr)
hdu2 = fits.BinTableHDU(data=image[1].data,header=image[1].header)
hdul = fits.HDUList([hdu1, hdu2])
hdul.writeto(path+'/median_g_after.fits',overwrite=True)

#Cropping the median image in a squared region around the target
gdata_crop = gdata[max(0,round(gxtar-zoom)):min(np.shape(gdata)[1]-1,round(gxtar+zoom)),max(0,round(gytar-zoom)):min(np.shape(gdata)[0]-1,round(gytar+zoom))]
gx = image[1].data['x']
gy = image[1].data['y']
ngsources = len(gx)
mg, sg = np.mean(gdata_crop), np.std(gdata_crop)

rcParams['figure.figsize'] = [9., 11.]
plt.rcParams['figure.dpi'] = 400
plt.rcParams['savefig.dpi'] = 400
fig, (ax1, ax2, ax3) = plt.subplots(3,1)
ax1.set_title(r''+FoVname+' FoV bkg subtracted, g band, '+str(ng)+' images and '+str(ngsources)+' sources', fontsize=12)
target = Circle(xy=(gxtar,gytar), radius=8.0)
target.set_linewidth(0.5)
#ax1.set_xlim(max(0,gxtar-zoom), min(np.shape(gdata)[1]-1,gxtar+zoom))
#ax1.set_ylim(max(0,gytar-zoom), min(np.shape(gdata)[1]-1,gytar+zoom))
ax1.imshow(gdata, interpolation='nearest', cmap='gray', vmin=mg-3*sg, vmax=mg+3*sg, origin='lower')
ax1.text(gxtar,gytar,FoVname,horizontalalignment='center',verticalalignment='center',color='y',fontsize=7)
target.set_facecolor('none')
target.set_edgecolor('yellow')
ax1.add_artist(target)
               
#for i in range(0,ngsources):
#    a = Circle(xy=(gx[i], gy[i]), radius=5.)
#    a.set_linewidth(0.5)
#    a.set_facecolor('none')
#    a.set_edgecolor('red')
#    ax1.add_artist(a)
#if os.path.exists(path+'/Sources_g_sep.reg'):
#      os.remove(path+'/Sources_g_sep.reg')
#reg = open(path+'/Sources_g_sep.reg','w')
#print("# Region file format: DS9 version 4.1",file=reg)
#print("global dashlist=8 3 width=1 font=\"helvetica 10 normal roman\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1",file=reg)
#print("text(2000,2000) # text={g filter} color=green width=1",file=reg)
#print("text(2000,1700) # text={%d images} color=green width=1" % ngimages,file=reg)
#print("text(2000,1400) # text={%d sources} color=green width=1" % ngsources,file=reg)
#for i in range(0,ngsources):
#    print("point(%f,%f) # point=cross size=0.01 # color=red width=1" % (image[1].data['x'][i],image[1].data['y'][i]),file=reg)
#reg.close()

image.close()

#rfilter field plot

image = fits.open(path+'/rfilter/median_r.fits')
hdr = image[0].header
rw = wcs.WCS(hdr)
rdata = image[0].data
rxtar, rytar = rw.all_world2pix(ra, dec, 1)
try:
    bkg = sep.Background(rdata, bw=8, bh=8, fw=3, fh=3)
except ValueError:
    rdata = rdata.byteswap(True).newbyteorder()
    bkg = sep.Background(rdata, bw=8, bh=8, fw=3, fh=3)
bkg.subfrom(rdata)
hdu1 = fits.PrimaryHDU(data=rdata,header=hdr)
hdu2 = fits.BinTableHDU(data=image[1].data,header=image[1].header)
hdul = fits.HDUList([hdu1, hdu2])
hdul.writeto(path+'/median_r_after.fits',overwrite=True)

rdata_crop = gdata[max(0,round(rxtar-zoom)):min(np.shape(rdata)[1]-1,round(rxtar+zoom)),max(0,round(rytar-zoom)):min(np.shape(rdata)[0]-1,round(rytar+zoom))]
rx = image[1].data['x']
ry = image[1].data['y']
nrsources = len(rx)
mr, sr = np.mean(rdata_crop), np.std(rdata_crop)
ax2.set_title(r''+FoVname+' FoV bkg subtracted, r band, '+str(nr)+' images and '+str(nrsources)+' sources', fontsize=12)
target = Circle(xy=(rxtar,rytar), radius=8)
target.set_linewidth(0.5)
ax2.imshow(rdata, interpolation='nearest', cmap='gray', vmin=mr-3.0*sr, vmax=mr+3.0*sr, origin='lower')
ax2.text(rxtar,rytar,FoVname,horizontalalignment='center',verticalalignment='center',color='y',fontsize=7)
target.set_facecolor('none')
target.set_edgecolor('yellow')
ax2.add_artist(target)
                                    
#for i in range(0,nrsources):
#    a = Circle(xy=(rx[i], ry[i]), radius=2.)
#    a.set_linewidth(1.0)
#    a.set_facecolor('none')
#    a.set_edgecolor('red')
#    ax2.add_artist(a)
#if os.path.exists(path+'/Sources_r_sep.reg'):
#      os.remove(path+'/Sources_r_sep.reg')
#reg = open(path+'/Sources_r_sep.reg','w')
#print("# Region file format: DS9 version 4.1",file=reg)
#print("global dashlist=8 3 width=1 font=\"helvetica 10 normal roman\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1",file=reg)
#print("text(2000,2000) # text={r filter} color=green width=1",file=reg)
#print("text(2000,1700) # text={%d images} color=green width=1" % nrimages,file=reg)
#print("text(2000,1400) # text={%d sources} color=green width=1" % nrsources,file=reg)
#for i in range(0,nrsources):
#    print("point(%f,%f) # point=cross size=0.01 # color=red width=1" % (image[1].data['x'][i],image[1].data['y'][i]),file=reg)
#reg.close()
image.close()

#ifilter field plot

image = fits.open(path+'/ifilter/median_i.fits')
hdr = image[0].header
iw = wcs.WCS(hdr)
idata = image[0].data
ixtar, iytar = iw.all_world2pix(ra, dec, 1)
try:
    bkg = sep.Background(idata, bw=8, bh=8, fw=3, fh=3)
except ValueError:
    idata = idata.byteswap(True).newbyteorder()
    bkg = sep.Background(idata, bw=8, bh=8, fw=3, fh=3)
bkg.subfrom(idata)
hdu1 = fits.PrimaryHDU(data=idata,header=hdr)
hdu2 = fits.BinTableHDU(data=image[1].data,header=image[1].header)
hdul = fits.HDUList([hdu1, hdu2])
hdul.writeto(path+'/median_i_after.fits',overwrite=True)
idata_crop = idata[max(0,round(ixtar-zoom)):min(np.shape(idata)[1]-1,round(ixtar+zoom)),max(0,round(iytar-zoom)):min(np.shape(idata)[0]-1,round(iytar+zoom))]
ix = image[1].data['x']
iy = image[1].data['y']
nisources = len(ix)
mi, si = np.mean(idata_crop), np.std(idata_crop)
ax3.set_title(r''+FoVname+' FoV bkg subtracted, i band, '+str(ni)+' images and '+str(nisources)+' sources', fontsize=12)
target = Circle(xy=(ixtar,iytar), radius=8)
target.set_linewidth(0.5)
ax3.imshow(idata, interpolation='nearest', cmap='gray', vmin=mi-3.0*si, vmax=mi+3.0*si, origin='lower')
ax3.text(ixtar,iytar,FoVname,horizontalalignment='center',verticalalignment='center',color='y',fontsize=7)
target.set_facecolor('none')
target.set_edgecolor('yellow')
ax3.add_artist(target)
                                    
#for i in range(0,nisources):
#    a = Circle(xy=(ix[i], iy[i]), radius=2.)
#    a.set_linewidth(1.0)
#    a.set_facecolor('none')
#    a.set_edgecolor('red')
#    ax3.add_artist(a)
#if os.path.exists(path+'/Sources_i_sep.reg'):
#      os.remove(path+'/Sources_i_sep.reg')
#reg = open(path+'/Sources_i_sep.reg','w')
#print("# Region file format: DS9 version 4.1",file=reg)
#print("global dashlist=8 3 width=1 font=\"helvetica 10 normal roman\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1",file=reg)
#print("text(2000,2000) # text={i filter} color=green width=1",file=reg)
#print("text(2000,1700) # text={%d images} color=green width=1" % niimages,file=reg)
#print("text(2000,1400) # text={%d sources} color=green width=1" % nisources,file=reg)
#for i in range(0,nisources):
#    print("point(%f,%f) # point=cross size=0.01 # color=red width=1" % (image[1].data['x'][i],image[1].data['y'][i]),file=reg)
#reg.close()
image.close()

plt.savefig(path+'/'+((FoVname.replace(" ","")).replace("-","m")).replace("+","p")+'_FoV.png')
plt.close(fig)

#gfilter comparison plot

rcParams['figure.figsize'] = [9., 11.]
plt.rcParams['figure.dpi'] = 400
plt.rcParams['savefig.dpi'] = 400
fig, (ax1, ax2, ax3) = plt.subplots(3,1)
if np.shape(gcomps)[0] > 1:
            ax1.set_title(r''+FoVname+' FoV bkg subtracted, g band, '+str(len(gcomps))+' comparison stars with $<\sigma>$=%.5f'%gsigma+' mag', fontsize=12)
            ax1.imshow(gdata, interpolation='nearest', cmap='gray', vmin=mg-3.0*sg, vmax=mg+3.0*sg, origin='lower')
            ax1.text(gxtar,gytar,FoVname,horizontalalignment='center',verticalalignment='center',color='y',fontsize=7)
            target = Circle(xy=(gxtar,gytar), radius=8)
            target.set_linewidth(0.5)
            target.set_facecolor('none')
            target.set_edgecolor('yellow')
            ax1.add_artist(target)
            for i in range(0,len(gcomps)):
                        #[xcomp, ycomp] = np.dot(Rinv,[gcomps[i,0]-raref, gcomps[i,1]-decref])
                        #[xcomp, ycomp] = [xcomp+xref, ycomp+yref]
                        xcomp, ycomp = gw.all_world2pix(gcomps[i,0], gcomps[i,1], 1)
                        a = Circle(xy=(xcomp, ycomp), radius=8.)
                        a.set_linewidth(0.5)
                        a.set_facecolor('none')
                        a.set_edgecolor('red')
                        ax1.add_artist(a)

#rfilter comparison plot
if np.shape(rcomps)[0] > 1:
            ax2.set_title(r''+FoVname+' FoV bkg subtracted, r band, '+str(len(rcomps))+' comparison stars with $<\sigma>$=%.5f'%rsigma+' mag', fontsize=12)
            ax2.imshow(rdata, interpolation='nearest', cmap='gray', vmin=mr-3.0*sr, vmax=mr+3.0*sr, origin='lower')
            ax2.text(rxtar,rytar,FoVname,horizontalalignment='center',verticalalignment='center',color='y',fontsize=7)
            target = Circle(xy=(rxtar,rytar), radius=8)
            target.set_linewidth(0.5)
            target.set_facecolor('none')
            target.set_edgecolor('yellow')
            ax2.add_artist(target)
            for i in range(0,len(rcomps)):
                        xcomp, ycomp = rw.all_world2pix(rcomps[i,0], rcomps[i,1], 1)
                        a = Circle(xy=(xcomp, ycomp), radius=8.)
                        a.set_linewidth(0.5)
                        a.set_facecolor('none')
                        a.set_edgecolor('red')
                        ax2.add_artist(a)

#ifilter comparison plot
if np.shape(icomps)[0] > 1:
            ax3.set_title(r''+FoVname+' FoV bkg subtracted, i band, '+str(len(icomps))+' comparison stars with $<\sigma>$=%.5f'%isigma+' mag', fontsize=12)
            ax3.imshow(idata, interpolation='nearest', cmap='gray', vmin=mi-3.0*si, vmax=mi+3.0*si, origin='lower')
            ax3.text(ixtar,iytar,FoVname,horizontalalignment='center',verticalalignment='center',color='y',fontsize=7)
            target = Circle(xy=(ixtar,iytar), radius=8)
            target.set_linewidth(0.5)
            target.set_facecolor('none')
            target.set_edgecolor('yellow')
            ax3.add_artist(target)
            for i in range(0,len(icomps)):
                        xcomp, ycomp = iw.all_world2pix(icomps[i,0], icomps[i,1], 1)
                        a = Circle(xy=(xcomp, ycomp), radius=8.)
                        a.set_linewidth(0.5)
                        a.set_facecolor('none')
                        a.set_edgecolor('red')
                        ax3.add_artist(a)

plt.savefig(path+'/'+((FoVname.replace(" ","")).replace("-","m")).replace("+","p")+'_compFoV.png')
plt.close(fig)

#For loop on the multiple nights
for i in range(0,N_nights):
            fig=plt.figure(figsize=(7.,14.), dpi=300)
            #gfilter comparison light curves
            if np.shape(gcomps)[0] > 1:
                        color = cm.rainbow(np.linspace(0, 1, len(gcomps)))
                        plt.subplot(3, 1, 1)
                        plt.title(r''+FoVname+' light curves of the '+str(len(gcomps))+' comparison stars in the g band\nNight of '+nights[i], fontsize=12)
                        plt.gca().invert_yaxis()
                        for j, c in zip(range(len(gcomps)), color):
                                    if os.path.exists(path+'/gfilter/outputcats/doerPhot_V'+str(j+1)+'.csv'):
                                                mjd, dm, sigma_dm = np.genfromtxt(path+'/gfilter/outputcats/doerPhot_V'+str(j+1)+'.csv', delimiter=',',usecols=(8,12,13),unpack=True)
                                                nightdim = np.zeros(N_nights,int)
                                                for k in range(0,(np.size(mjd))):
                                                            if k!=0:
                                                                        if (mjd[k]-mjd[k-1]) > 0.8:
                                                                                    n_night = np.argmin(np.fabs(mjd[k-1]-gmjd_ref))
                                                                                    nightdim[n_night] = k-nightdim[0]
                                                                                    for l in range(1,n_night):
                                                                                                nightdim[n_night] = nightdim[n_night]-nightdim[l]
                                                                        if k==np.size(mjd)-1:
                                                                                    n_night = np.argmin(np.fabs(mjd[k]-gmjd_ref))
                                                                                    nightdim[n_night] = k+1-nightdim[0]
                                                                                    for l in range(1,n_night):
                                                                                                nightdim[n_night] = nightdim[n_night]-nightdim[l]
                                                mjd_single = []
                                                dm_single = []
                                                sigma_dm_single = []
                                                for k in range(0,N_nights):
                                                            mjd_single.append(mjd[np.sum(nightdim[0:k]):np.sum(nightdim[0:k])+nightdim[k]])
                                                            dm_single.append(dm[np.sum(nightdim[0:k]):np.sum(nightdim[0:k])+nightdim[k]])
                                                            sigma_dm_single.append(sigma_dm[np.sum(nightdim[0:k]):np.sum(nightdim[0:k])+nightdim[k]])
                                                if j==0:
                                                            #plt.xlabel(r"Time (MJD-%d)"%np.trunc(mjd_single[0]),fontsize=8)
                                                            plt.ylabel(r"Differential g magnitude",fontsize=8)
                                                if np.size(mjd_single[i])!=0:
                                                   plt.plot(mjd_single[i]-np.trunc(mjd_single[i][0]), dm_single[i], c=c, marker='o', markersize=2, linestyle='None')
                                                   plt.errorbar(mjd_single[i]-np.trunc(mjd_single[i][0]), dm_single[i], sigma_dm_single[i], c=c, marker='o', markersize=2, linestyle='None')

            #rfilter comparison light curves

            if np.shape(rcomps)[0] > 1:
                        color = cm.rainbow(np.linspace(0, 1, len(rcomps)))
                        plt.subplot(3, 1, 2)
                        plt.title(r''+FoVname+' light curves of the '+str(len(rcomps))+' comparison stars in the r band\nNight of '+nights[i], fontsize=12)
                        plt.gca().invert_yaxis()
                        for j, c in zip(range(len(rcomps)), color):
                                    if os.path.exists(path+'/rfilter/outputcats/doerPhot_V'+str(j+1)+'.csv'):
                                                mjd, dm, sigma_dm = np.genfromtxt(path+'/rfilter/outputcats/doerPhot_V'+str(j+1)+'.csv', delimiter=',',usecols=(8,12,13),unpack=True)
                                                nightdim = np.zeros(N_nights,int)
                                                for k in range(0,(np.size(mjd))):
                                                            if k!=0:
                                                                        if (mjd[k]-mjd[k-1]) > 0.8:
                                                                                    n_night = np.argmin(np.fabs(mjd[k-1]-rmjd_ref))
                                                                                    nightdim[n_night] = k-nightdim[0]
                                                                                    for l in range(1,n_night):
                                                                                                nightdim[n_night] = nightdim[n_night]-nightdim[l]
                                                                        if k==np.size(mjd)-1:
                                                                                    n_night = np.argmin(np.fabs(mjd[k]-rmjd_ref))
                                                                                    nightdim[n_night] = k+1-nightdim[0]
                                                                                    for l in range(1,n_night):
                                                                                                nightdim[n_night] = nightdim[n_night]-nightdim[l]
                                                mjd_single = []
                                                dm_single = []
                                                sigma_dm_single = []
                                                for k in range(0,N_nights):
                                                            mjd_single.append(mjd[np.sum(nightdim[0:k]):np.sum(nightdim[0:k])+nightdim[k]])
                                                            dm_single.append(dm[np.sum(nightdim[0:k]):np.sum(nightdim[0:k])+nightdim[k]])
                                                            sigma_dm_single.append(sigma_dm[np.sum(nightdim[0:k]):np.sum(nightdim[0:k])+nightdim[k]])
                                                if j==0:
                                                            #plt.xlabel(r"Time (MJD-%d)"%np.trunc(mjd_single[0]),fontsize=8)
                                                            plt.ylabel(r"Differential r magnitude",fontsize=8)
                                                if np.size(mjd_single[i])!=0:
                                                   plt.plot(mjd_single[i]-np.trunc(mjd_single[i][0]), dm_single[i], c=c, marker='o', markersize=2, linestyle='None')
                                                   plt.errorbar(mjd_single[i]-np.trunc(mjd_single[i][0]), dm_single[i], sigma_dm_single[i], c=c, marker='o', markersize=2, linestyle='None')

            #ifilter comparison light curves

            if np.shape(icomps)[0] > 1:
                        color = cm.rainbow(np.linspace(0, 1, len(icomps)))
                        plt.subplot(3, 1, 3)
                        plt.title(r''+FoVname+' light curves of the '+str(len(icomps))+' comparison stars in the i band\nNight of '+nights[i], fontsize=12)
                        plt.gca().invert_yaxis()
                        for j, c in zip(range(len(icomps)), color):
                                    if os.path.exists(path+'/ifilter/outputcats/doerPhot_V'+str(j+1)+'.csv'):
                                                mjd, dm, sigma_dm = np.genfromtxt(path+'/ifilter/outputcats/doerPhot_V'+str(j+1)+'.csv', delimiter=',',usecols=(8,12,13),unpack=True)
                                                nightdim = np.zeros(N_nights,int)
                                                for k in range(0,(np.size(mjd))):
                                                            if k!=0:
                                                                        if (mjd[k]-mjd[k-1]) > 0.8:
                                                                                    n_night = np.argmin(np.fabs(mjd[k-1]-imjd_ref))
                                                                                    nightdim[n_night] = k-nightdim[0]
                                                                                    for l in range(1,n_night):
                                                                                                nightdim[n_night] = nightdim[n_night]-nightdim[l]
                                                                        if k==np.size(mjd)-1:
                                                                                    n_night = np.argmin(np.fabs(mjd[k]-imjd_ref))
                                                                                    nightdim[n_night] = k+1-nightdim[0]
                                                                                    for l in range(1,n_night):
                                                                                                nightdim[n_night] = nightdim[n_night]-nightdim[l]
                                                mjd_single = []
                                                dm_single = []
                                                sigma_dm_single = []
                                                for k in range(0,N_nights):
                                                            mjd_single.append(mjd[np.sum(nightdim[0:k]):np.sum(nightdim[0:k])+nightdim[k]])
                                                            dm_single.append(dm[np.sum(nightdim[0:k]):np.sum(nightdim[0:k])+nightdim[k]])
                                                            sigma_dm_single.append(sigma_dm[np.sum(nightdim[0:k]):np.sum(nightdim[0:k])+nightdim[k]])
                                                if j==0:
                                                            plt.xlabel(r"Time (MJD-%d)"%np.trunc(mjd_single[i][0]),fontsize=8)
                                                            plt.ylabel(r"Differential i magnitude",fontsize=8)
                                                if np.size(mjd_single[i])!=0:
                                                   plt.plot(mjd_single[i]-np.trunc(mjd_single[i][0]), dm_single[i], c=c, marker='o', markersize=2, linestyle='None')
                                                   plt.errorbar(mjd_single[i]-np.trunc(mjd_single[i][0]), dm_single[i], sigma_dm_single[i], c=c, marker='o', markersize=2, linestyle='None')

            plt.savefig(path+'/'+((FoVname.replace(" ","")).replace("-","m")).replace("+","p")+'_complcs_night'+str(i+1)+'.png')
            plt.close(fig)

