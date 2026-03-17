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

search = input("Enter 'hunt' or 'target' for the type of search (default value 'hunt'):\n") or 'hunt'
while search!='hunt' and search!='target':
   search = input("Invalid value, please insert 'hunt' or 'target':\n") or 'hunt'

analysis = input("Enter 'STELLA', 'INT' or 'LCO' for the type of data (default value 'LCO'):\n") or 'LCO'
while analysis!='STELLA' and analysis!='INT' and analysis!='LCO':
   analysis = input("Invalid value, please insert 'STELLA', 'INT' or 'LCO':\n") or 'LCO'

if analysis == 'INT':
            cam = path.split('/')[-1][-1]

#Defining scale factor to use as squared region around the target, with side=2*bound*semi_axis
bound = 2
#Defining radius size (in arcmin) to use for plotting 4FGL sources nearby to our candidate target
radius = 11

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

if analysis == 'STELLA':
            dates = [s.split("-")[0].split("science")[1].split("A")[0] for s in files]
elif analysis == 'INT':
            dates = []
            for i in range(0,len(files)):
                        hdul = fits.open(path+'/rfilter/'+files[i])
                        hdr = hdul[0].header
                        if hdr['UTOBS'][0]=='0':
                                    dates.append(hdr['DATE-OBS'][:-1]+str(int(hdr['DATE-OBS'][-1])-1))
                                    if dates[i][-1]=='0':
                                                dates[i] = dates[i][:-2]+'31'
                        else:
                                    dates.append(hdr['DATE-OBS'])
                        hdul.close()
elif analysis=='LCO':
            dates = [s.split("-")[2] for s in files]
            
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

#Checking nights mjd_ref in g, r and i filters
mjd = [float(s.split('_')[2]) for s in gimages if 'median' not in s]
gmjd_ref = np.zeros(N_nights)
counter = 0
for i in range(0,np.size(mjd)):
            if i!=0:
                        if (mjd[i]-mjd[i-1]) > 0.7:
                                    gmjd_ref[counter] = mjd[i-1]
                                    counter+=1
            if i==np.size(mjd)-1:
                        gmjd_ref[counter] = mjd[i]
mjd = [float(s.split('_')[2]) for s in rimages if 'median' not in s]
rmjd_ref = np.zeros(N_nights)
counter = 0
for i in range(0,np.size(mjd)):
            if i!=0:
                        if (mjd[i]-mjd[i-1]) > 0.7:
                                    rmjd_ref[counter] = mjd[i-1]
                                    counter+=1
            if i==np.size(mjd)-1:
                        rmjd_ref[counter] = mjd[i]
mjd = [float(s.split('_')[2]) for s in iimages if 'median' not in s]
imjd_ref = np.zeros(N_nights)
counter = 0
for i in range(0,np.size(mjd)):
            if i!=0:
                        if (mjd[i]-mjd[i-1]) > 0.7:
                                    imjd_ref[counter] = mjd[i-1]
                                    counter+=1
            if i==np.size(mjd)-1:
                        imjd_ref[counter] = mjd[i]

#gfilter field plot

image = fits.open(path+'/gfilter/median_g.fits')
hdr = image[0].header
if analysis=='STELLA':
   FoVname = hdr['OBJNAME']
   pixscale = float(hdr['PIXSCALE'])
elif analysis=='INT':
   FoVname = hdr['OBJECT']
   pixscale = float(hdr['SECPPIX'])
elif analysis=='LCO':
   FoVname = hdr['OBJECT'].replace("L_","L")
   pixscale = float(hdr['PIXSCALE'])
   
gw = wcs.WCS(hdr)
gdata = image[0].data

if search == 'hunt':
   #Load the corresponding Fermi coordinates and ellipses (in degrees) of the target
   hdul = fits.open('/export/work/marcotu/gll_psc_v28.fit')
   [[ra, dec, ell_a, ell_b, ell_theta]] = [[hdul[1].data['RAJ2000'][i],hdul[1].data['DEJ2000'][i],3600./pixscale*hdul[1].data['Conf_95_SemiMajor'][i],3600./pixscale*hdul[1].data['Conf_95_SemiMinor'][i],hdul[1].data['Conf_95_PosAng'][i]] for i in range(0,len(hdul[1].data)) if FoVname.replace("L","L ")==hdul[1].data['Source_Name'][i]]
   hdul.close()

   #Load the 3FGL Fermi updated coordinates for our target and eventual closeby sources (within radius arcmin from our target), also extended sources
   hdul = fits.open('/export/work/marcotu/gll_psc_v16.fit')
   sources_3FGL = [hdul[1].data['Source_Name'][i] for i in range(0,len(hdul[1].data))]
   coord_3FGL = [[hdul[1].data['RAJ2000'][i],hdul[1].data['DEJ2000'][i],hdul[1].data['Conf_95_SemiMajor'][i],hdul[1].data['Conf_95_SemiMinor'][i],hdul[1].data['Conf_95_PosAng'][i]] for i in range(0,len(hdul[1].data))]
   sources_3FGL = np.asarray(sources_3FGL)
   coord_3FGL = np.asarray(coord_3FGL)

   distance = SkyCoord(ra, dec, frame='fk5', unit='deg').separation(SkyCoord(coord_3FGL[:,0],coord_3FGL[:,1],frame='fk5',unit='deg')).arcminute
   sources_3FGL = [sources_3FGL[i] for i in range(0,len(sources_3FGL)) if distance[i]<=radius]
   sources_3FGL = np.asarray(sources_3FGL)
   coord_3FGL = [coord_3FGL[i,:] for i in range(0,len(coord_3FGL)) if distance[i]<=radius]
   coord_3FGL = np.asarray(coord_3FGL)
   distance = [distance[i] for i in range(0,len(distance)) if distance[i]<=radius]

   #Load eventual 4FGL extended closeby sources (within radius_gamma arcmin from our target)
   #sources_4FGLe = [[hdul[2].data['Source_Name'][i],hdul[2].data['RAJ2000'][i],hdul[2].data['DEJ2000'][i],hdul[2].data['Model_SemiMajor'][i],hdul[2].data['Model_SemiMinor'][i],hdul[2].data['Model_PosAng'][i]] for i in range(0,len(hdul[2].data))]
   #distance = SkyCoord(ra, dec, frame='fk5', unit='deg').separation(SkyCoord(np.asarray(sources_4FGLe)[:,1],np.asarray(sources_4FGLe)[:,2],frame='fk5',unit='deg')).arcminute
   #sources_4FGLe = [sources_4FGLe[i] for i in range(0,len(sources_4FGLe)) if distance[i]<=radius]
   #distance = [distance[i] for i in range(0,len(distance)) if distance[i]<=radius]
   hdul.close()

else:
   #Load the target coordinates from the header of the fits
   if analysis == 'STELLA':
      ra = hdr['OBJRA']
      dec = hdr['OBJDEC']
   elif analysis == 'INT':
      ra = (float(hdr['RA'].split(':')[0])+float(hdr['RA'].split(':')[1])/60.+float(hdr['RA'].split(':')[2])/3600.)*360./24.
      if float(hdr['DEC'].split(':')[0])>0.:
         dec = (float(hdr['DEC'].split(':')[0])+float(hdr['DEC'].split(':')[1])/60.+float(hdr['DEC'].split(':')[2])/3600.)
      else:
         dec = (float(hdr['DEC'].split(':')[0])-float(hdr['DEC'].split(':')[1])/60.-float(hdr['DEC'].split(':')[2])/3600.)
   else:
      ra = (float(hdr['CAT-RA'].split(':')[0])+float(hdr['CAT-RA'].split(':')[1])/60.+float(hdr['CAT-RA'].split(':')[2])/3600.)*360./24.
      if float(hdr['CAT-DEC'].split(':')[0])>0.:
         dec = (float(hdr['CAT-DEC'].split(':')[0])+float(hdr['CAT-DEC'].split(':')[1])/60.+float(hdr['CAT-DEC'].split(':')[2])/3600.)
      else:
         dec = (float(hdr['CAT-DEC'].split(':')[0])-float(hdr['CAT-DEC'].split(':')[1])/60.-float(hdr['CAT-DEC'].split(':')[2])/3600.)

gxtar, gytar = gw.all_world2pix(ra, dec, 1)

try:
    bkg = sep.Background(gdata, bw=16, bh=16, fw=3, fh=3)
except ValueError:
    gdata = gdata.byteswap(True).newbyteorder()
    bkg = sep.Background(gdata, bw=16, bh=16, fw=3, fh=3)
bkg.subfrom(gdata)
hdu1 = fits.PrimaryHDU(data=gdata,header=hdr)
hdu2 = fits.BinTableHDU(data=image[1].data,header=image[1].header)
hdul = fits.HDUList([hdu1, hdu2])
hdul.writeto(path+'/median_g_after.fits',overwrite=True)
if search == 'hunt':
   #Cropping the median image in a squared region around the target
   gdata_crop = gdata[max(0,round(gxtar-bound*max(ell_a,ell_b))):min(np.shape(gdata)[1]-1,round(gxtar+bound*max(ell_a,ell_b))),max(0,round(gytar-bound*max(ell_a,ell_b))):min(np.shape(gdata)[0]-1,round(gytar+bound*max(ell_a,ell_b)))]
   gx = [image[1].data['x'][i] for i in range(0,len(image[1].data)) if image[1].data['x'][i]>=(gxtar-bound*max(ell_a,ell_b)) and image[1].data['x'][i]<=(gxtar+bound*max(ell_a,ell_b)) and image[1].data['y'][i]>=(gytar-bound*max(ell_a,ell_b)) and image[1].data['y'][i]<=(gytar+bound*max(ell_a,ell_b))]
   gy = [image[1].data['y'][i] for i in range(0,len(image[1].data)) if image[1].data['x'][i]>=(gxtar-bound*max(ell_a,ell_b)) and image[1].data['x'][i]<=(gxtar+bound*max(ell_a,ell_b)) and image[1].data['y'][i]>=(gytar-bound*max(ell_a,ell_b)) and image[1].data['y'][i]<=(gytar+bound*max(ell_a,ell_b))]
else:
   gdata_crop = gdata
   gx = image[1].data['x']
   gy = image[1].data['y']
ngsources = len(gx)
mg, sg = np.mean(gdata_crop), np.std(gdata_crop)

rcParams['figure.figsize'] = [9., 11.]
plt.rcParams['figure.dpi'] = 400
plt.rcParams['savefig.dpi'] = 400
fig, (ax1, ax2, ax3) = plt.subplots(3,1)
if search == 'hunt':
   ax1.set_title(r''+FoVname+' FoV bkg subtracted, g band, '+str(ng)+' images and '+str(ngsources)+' sources in the '+str(2*bound)+'$\sigma$ square', fontsize=12)
   ax1.set_xlim(max(0,gxtar-bound*max(ell_a,ell_b)), min(np.shape(gdata)[1]-1,gxtar+bound*max(ell_a,ell_b)))
   ax1.set_ylim(max(0,gytar-bound*max(ell_a,ell_b)), min(np.shape(gdata)[0]-1,gytar+bound*max(ell_a,ell_b)))
   target = Ellipse(xy=(gxtar,gytar), width=2*ell_a, height=2*ell_b, angle=ell_theta)
else:
   ax1.set_title(r''+FoVname+' FoV bkg subtracted, g band, '+str(ng)+' images and '+str(ngsources)+' sources', fontsize=12)
   target = Circle(xy=(gxtar,gytar), radius=20)
   target.set_linewidth(1.0)
ax1.imshow(gdata, interpolation='nearest', cmap='gray', vmin=mg-0.15*sg, vmax=mg+0.15*sg, origin='lower')
ax1.text(gxtar,gytar,FoVname,horizontalalignment='center',verticalalignment='center',color='cyan',fontsize=7)
target.set_facecolor('none')
target.set_edgecolor('cyan')
ax1.add_artist(target)

if search == 'hunt':
   if np.size(sources_3FGL)!=0:
      x_3FGL, y_3FGL = gw.all_world2pix(coord_3FGL[:,0], coord_3FGL[:,1], 1)
      for i in range(0,np.shape(sources_3FGL)[0]):
         if x_3FGL[i]>=max(0,gxtar-bound*max(ell_a,ell_b)) and x_3FGL[i]<=min(np.shape(gdata)[1]-1,gxtar+bound*max(ell_a,ell_b)) and y_3FGL[i]>=max(0,gytar-bound*max(ell_a,ell_b)) and y_3FGL[i]<=min(np.shape(gdata)[0]-1,gytar+bound*max(ell_a,ell_b)):
            ax1.text(x_3FGL[i],y_3FGL[i],sources_3FGL[i].replace("L ","L"),horizontalalignment='center',verticalalignment='center',color='yellow',fontsize=7)
            a = Ellipse(xy=(x_3FGL[i], y_3FGL[i]), width=2*3600./pixscale*coord_3FGL[i,2], height=2*3600./pixscale*coord_3FGL[i,3], angle=coord_3FGL[i,4])
            a.set_facecolor('none')
            a.set_edgecolor('yellow')
            ax1.add_artist(a)
   #if np.size(sources_4FGLe)!=0:
   #   x_4FGLe, y_4FGLe = gw.all_world2pix(np.asarray(sources_4FGLe)[:,1], np.asarray(sources_4FGLe)[:,2], 1)
   #   for i in range(0,np.shape(sources_4FGLe)[0]):
   #      if x_4FGLe[i]>=max(0,gxtar-bound*max(ell_a,ell_b)) and x_4FGLe[i]<=min(np.shape(gdata)[1]-1,gxtar+bound*max(ell_a,ell_b)) and y_4FGLe[i]>=max(0,gytar-bound*max(ell_a,ell_b)) and y_4FGLe[i]<=min(np.shape(gdata)[0]-1,gytar+bound*max(ell_a,ell_b)):
   #         ax1.text(x_4FGLe[i],y_4FGLe[i],sources_4FGLe[i][0].replace("L ","L"),horizontalalignment='center',verticalalignment='center',color='cyan',fontsize=7)
   #         a = Ellipse(xy=(x_4FGLe[i], y_4FGLe[i]), width=2*3600./pixscale*np.asarray(sources_4FGLe)[i,3], height=2*3600./pixscale*np.asarray(sources_4FGLe)[i,4], angle=np.asarray(sources_4FGLe)[i,5])
   #         a.set_facecolor('none')
   #         a.set_edgecolor('cyan')
   #         ax1.add_artist(a)
               
for i in range(0,ngsources):
    a = Circle(xy=(gx[i], gy[i]), radius=2.)
    a.set_linewidth(1.0)
    a.set_facecolor('none')
    a.set_edgecolor('red')
    ax1.add_artist(a)
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
    bkg = sep.Background(rdata, bw=16, bh=16, fw=3, fh=3)
except ValueError:
    rdata = rdata.byteswap(True).newbyteorder()
    bkg = sep.Background(rdata, bw=16, bh=16, fw=3, fh=3)
bkg.subfrom(rdata)
hdu1 = fits.PrimaryHDU(data=rdata,header=hdr)
hdu2 = fits.BinTableHDU(data=image[1].data,header=image[1].header)
hdul = fits.HDUList([hdu1, hdu2])
hdul.writeto(path+'/median_r_after.fits',overwrite=True)
if search == 'hunt':
   rdata_crop = rdata[max(0,round(rxtar-bound*max(ell_a,ell_b))):min(np.shape(rdata)[1]-1,round(rxtar+bound*max(ell_a,ell_b))),max(0,round(rytar-bound*max(ell_a,ell_b))):min(np.shape(rdata)[0]-1,round(rytar+bound*max(ell_a,ell_b)))]
   rx = [image[1].data['x'][i] for i in range(0,len(image[1].data)) if image[1].data['x'][i]>=(rxtar-bound*max(ell_a,ell_b)) and image[1].data['x'][i]<=(rxtar+bound*max(ell_a,ell_b)) and image[1].data['y'][i]>=(rytar-bound*max(ell_a,ell_b)) and image[1].data['y'][i]<=(rytar+bound*max(ell_a,ell_b))]
   ry = [image[1].data['y'][i] for i in range(0,len(image[1].data)) if image[1].data['x'][i]>=(rxtar-bound*max(ell_a,ell_b)) and image[1].data['x'][i]<=(rxtar+bound*max(ell_a,ell_b)) and image[1].data['y'][i]>=(rytar-bound*max(ell_a,ell_b)) and image[1].data['y'][i]<=(rytar+bound*max(ell_a,ell_b))]
else:
   rdata_crop = rdata
   rx = image[1].data['x']
   ry = image[1].data['y']
nrsources = len(rx)
mr, sr = np.mean(rdata_crop), np.std(rdata_crop)
if search == 'hunt':
   ax2.set_title(r''+FoVname+' FoV bkg subtracted, r band, '+str(nr)+' images and '+str(nrsources)+' sources in the '+str(2*bound)+'$\sigma$ square', fontsize=12)
   ax2.set_xlim(max(0,rxtar-bound*max(ell_a,ell_b)), min(np.shape(rdata)[1]-1,rxtar+bound*max(ell_a,ell_b)))
   ax2.set_ylim(max(0,rytar-bound*max(ell_a,ell_b)), min(np.shape(rdata)[0]-1,rytar+bound*max(ell_a,ell_b)))
   target = Ellipse(xy=(rxtar,rytar), width=2*ell_a, height=2*ell_b, angle=ell_theta)
else:
   ax2.set_title(r''+FoVname+' FoV bkg subtracted, r band, '+str(nr)+' images and '+str(nrsources)+' sources', fontsize=12)
   target = Circle(xy=(rxtar,rytar), radius=20)
   target.set_linewidth(1.0)
ax2.imshow(rdata, interpolation='nearest', cmap='gray', vmin=mr-0.15*sr, vmax=mr+0.15*sr, origin='lower')
ax2.text(rxtar,rytar,FoVname,horizontalalignment='center',verticalalignment='center',color='cyan',fontsize=7)
target.set_facecolor('none')
target.set_edgecolor('cyan')
ax2.add_artist(target)

if search == 'hunt':
   if np.size(sources_3FGL)!=0:
      x_3FGL, y_3FGL = rw.all_world2pix(coord_3FGL[:,0], coord_3FGL[:,1], 1)
      for i in range(0,np.shape(sources_3FGL)[0]):
         if x_3FGL[i]>=max(0,rxtar-bound*max(ell_a,ell_b)) and x_3FGL[i]<=min(np.shape(rdata)[1]-1,rxtar+bound*max(ell_a,ell_b)) and y_3FGL[i]>=max(0,rytar-bound*max(ell_a,ell_b)) and y_3FGL[i]<=min(np.shape(rdata)[0]-1,rytar+bound*max(ell_a,ell_b)):
            ax2.text(x_3FGL[i],y_3FGL[i],sources_3FGL[i].replace("L ","L"),horizontalalignment='center',verticalalignment='center',color='yellow',fontsize=7)
            a = Ellipse(xy=(x_3FGL[i], y_3FGL[i]), width=2*3600./pixscale*coord_3FGL[i,2], height=2*3600./pixscale*coord_3FGL[i,3], angle=coord_3FGL[i,4])
            a.set_facecolor('none')
            a.set_edgecolor('yellow')
            ax2.add_artist(a)
   #if np.size(sources_4FGLe)!=0:
   #   x_4FGLe, y_4FGLe = rw.all_world2pix(np.asarray(sources_4FGLe)[:,1], np.asarray(sources_4FGLe)[:,2], 1)
   #   for i in range(0,np.shape(sources_4FGLe)[0]):
   #      if x_4FGLe[i]>=max(0,rxtar-bound*max(ell_a,ell_b)) and x_4FGLe[i]<=min(np.shape(rdata)[1]-1,rxtar+bound*max(ell_a,ell_b)) and y_4FGLe[i]>=max(0,rytar-bound*max(ell_a,ell_b)) and y_4FGLe[i]<=min(np.shape(rdata)[0]-1,rytar+bound*max(ell_a,ell_b)):
   #         ax2.text(x_4FGLe[i],y_4FGLe[i],sources_4FGLe[i][0].replace("L ","L"),horizontalalignment='center',verticalalignment='center',color='cyan',fontsize=7)
   #         a = Ellipse(xy=(x_4FGLe[i], y_4FGLe[i]), width=2*3600./pixscale*np.asarray(sources_4FGLe)[i,3], height=2*3600./pixscale*np.asarray(sources_4FGLe)[i,4], angle=np.asarray(sources_4FGLe)[i,5])
   #         a.set_facecolor('none')
   #         a.set_edgecolor('cyan')
   #         ax2.add_artist(a)
                                    
for i in range(0,nrsources):
    a = Circle(xy=(rx[i], ry[i]), radius=2.)
    a.set_linewidth(1.0)
    a.set_facecolor('none')
    a.set_edgecolor('red')
    ax2.add_artist(a)
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
    bkg = sep.Background(idata, bw=16, bh=16, fw=3, fh=3)
except ValueError:
    idata = idata.byteswap(True).newbyteorder()
    bkg = sep.Background(idata, bw=16, bh=16, fw=3, fh=3)
bkg.subfrom(idata)
hdu1 = fits.PrimaryHDU(data=idata,header=hdr)
hdu2 = fits.BinTableHDU(data=image[1].data,header=image[1].header)
hdul = fits.HDUList([hdu1, hdu2])
hdul.writeto(path+'/median_i_after.fits',overwrite=True)
if search == 'hunt':
   idata_crop = idata[max(0,round(ixtar-bound*max(ell_a,ell_b))):min(np.shape(idata)[1]-1,round(ixtar+bound*max(ell_a,ell_b))),max(0,round(iytar-bound*max(ell_a,ell_b))):min(np.shape(idata)[0]-1,round(iytar+bound*max(ell_a,ell_b)))]
   ix = [image[1].data['x'][i] for i in range(0,len(image[1].data)) if image[1].data['x'][i]>=(ixtar-bound*max(ell_a,ell_b)) and image[1].data['x'][i]<=(ixtar+bound*max(ell_a,ell_b)) and image[1].data['y'][i]>=(iytar-bound*max(ell_a,ell_b)) and image[1].data['y'][i]<=(iytar+bound*max(ell_a,ell_b))]
   iy = [image[1].data['y'][i] for i in range(0,len(image[1].data)) if image[1].data['x'][i]>=(ixtar-bound*max(ell_a,ell_b)) and image[1].data['x'][i]<=(ixtar+bound*max(ell_a,ell_b)) and image[1].data['y'][i]>=(iytar-bound*max(ell_a,ell_b)) and image[1].data['y'][i]<=(iytar+bound*max(ell_a,ell_b))]
else:
   idata_crop = idata
   ix = image[1].data['x']
   iy = image[1].data['y']
nisources = len(ix)
mi, si = np.mean(idata_crop), np.std(idata_crop)
if search == 'hunt':
   ax3.set_title(r''+FoVname+' FoV bkg subtracted, i band, '+str(ni)+' images and '+str(nisources)+' sources in the '+str(2*bound)+'$\sigma$ square', fontsize=12)
   ax3.set_xlim(max(0,ixtar-bound*max(ell_a,ell_b)), min(np.shape(idata)[1]-1,ixtar+bound*max(ell_a,ell_b)))
   ax3.set_ylim(max(0,iytar-bound*max(ell_a,ell_b)), min(np.shape(idata)[0]-1,iytar+bound*max(ell_a,ell_b)))
   target = Ellipse(xy=(ixtar,iytar), width=2*ell_a, height=2*ell_b, angle=ell_theta)
else:
   ax3.set_title(r''+FoVname+' FoV bkg subtracted, i band, '+str(ni)+' images and '+str(nisources)+' sources', fontsize=12)
   target = Circle(xy=(ixtar,iytar), radius=20)
   target.set_linewidth(1.0)
ax3.imshow(idata, interpolation='nearest', cmap='gray', vmin=mi-0.15*si, vmax=mi+0.15*si, origin='lower')
ax3.text(ixtar,iytar,FoVname,horizontalalignment='center',verticalalignment='center',color='cyan',fontsize=7)
target.set_facecolor('none')
target.set_edgecolor('cyan')
ax3.add_artist(target)

if search == 'hunt':
   if np.size(sources_3FGL)!=0:
      x_3FGL, y_3FGL = iw.all_world2pix(coord_3FGL[:,0], coord_3FGL[:,1], 1)
      for i in range(0,np.shape(sources_3FGL)[0]):
         if x_3FGL[i]>=max(0,ixtar-bound*max(ell_a,ell_b)) and x_3FGL[i]<=min(np.shape(idata)[1]-1,ixtar+bound*max(ell_a,ell_b)) and y_3FGL[i]>=max(0,iytar-bound*max(ell_a,ell_b)) and y_3FGL[i]<=min(np.shape(idata)[0]-1,iytar+bound*max(ell_a,ell_b)):
            ax3.text(x_3FGL[i],y_3FGL[i],sources_3FGL[i].replace("L ","L"),horizontalalignment='center',verticalalignment='center',color='yellow',fontsize=7)
            a = Ellipse(xy=(x_3FGL[i], y_3FGL[i]), width=2*3600./pixscale*coord_3FGL[i,2], height=2*3600./pixscale*coord_3FGL[i,3], angle=coord_3FGL[i,4])
            a.set_facecolor('none')
            a.set_edgecolor('yellow')
            ax3.add_artist(a)
   #if np.size(sources_4FGLe)!=0:
   #   x_4FGLe, y_4FGLe = iw.all_world2pix(np.asarray(sources_4FGLe)[:,1], np.asarray(sources_4FGLe)[:,2], 1)
   #   for i in range(0,np.shape(sources_4FGLe)[0]):
   #      if x_4FGLe[i]>=max(0,ixtar-bound*max(ell_a,ell_b)) and x_4FGLe[i]<=min(np.shape(idata)[1]-1,ixtar+bound*max(ell_a,ell_b)) and y_4FGLe[i]>=max(0,iytar-bound*max(ell_a,ell_b)) and y_4FGLe[i]<=min(np.shape(idata)[0]-1,iytar+bound*max(ell_a,ell_b)):
   #         ax3.text(x_4FGLe[i],y_4FGLe[i],sources_4FGLe[i][0].replace("L ","L"),horizontalalignment='center',verticalalignment='center',color='cyan',fontsize=7)
   #         a = Ellipse(xy=(x_4FGLe[i], y_4FGLe[i]), width=2*3600./pixscale*np.asarray(sources_4FGLe)[i,3], height=2*3600./pixscale*np.asarray(sources_4FGLe)[i,4], angle=np.asarray(sources_4FGLe)[i,5])
   #         a.set_facecolor('none')
   #         a.set_edgecolor('cyan')
   #         ax3.add_artist(a)
                                    
for i in range(0,nisources):
    a = Circle(xy=(ix[i], iy[i]), radius=2.)
    a.set_linewidth(1.0)
    a.set_facecolor('none')
    a.set_edgecolor('red')
    ax3.add_artist(a)
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

if analysis == 'INT':
   if search == 'hunt':
      plt.savefig(path+'/'+((FoVname.replace(" ","")).replace("-","m")).replace("+","p")+'_FoV_'+str(2*bound)+'sigma_cam'+cam+'.png')
   else:
      plt.savefig(path+'/'+((FoVname.replace(" ","")).replace("-","m")).replace("+","p")+'_FoV_cam'+cam+'.png')
else:
   if search == 'hunt':
      plt.savefig(path+'/'+((FoVname.replace(" ","")).replace("-","m")).replace("+","p")+'_FoV_'+str(2*bound)+'sigma.png')
   else:
      plt.savefig(path+'/'+((FoVname.replace(" ","")).replace("-","m")).replace("+","p")+'_FoV.png')
plt.close(fig)

#gfilter comparison plot

rcParams['figure.figsize'] = [9., 11.]
plt.rcParams['figure.dpi'] = 400
plt.rcParams['savefig.dpi'] = 400
fig, (ax1, ax2, ax3) = plt.subplots(3,1)
if np.shape(gcomps)[0] > 1:
            ax1.set_title(r''+FoVname+' FoV bkg subtracted, g band, '+str(len(gcomps))+' comparison stars with $<\sigma>$=%.5f'%gsigma+' mag', fontsize=12)
            ax1.imshow(gdata, interpolation='nearest', cmap='gray', vmin=mg-0.15*sg, vmax=mg+0.15*sg, origin='lower')
            ax1.text(gxtar,gytar,FoVname,horizontalalignment='center',verticalalignment='center',color='cyan',fontsize=7)
            if search == 'hunt':
               target = Ellipse(xy=(gxtar, gytar), width=2*ell_a, height=2*ell_b, angle=ell_theta)
            else:
               target = Circle(xy=(gxtar,gytar), radius=20)
               target.set_linewidth(1.0)
            target.set_facecolor('none')
            target.set_edgecolor('cyan')
            ax1.add_artist(target)
            if search == 'hunt':
               if np.size(sources_3FGL)!=0:
                  x_3FGL, y_3FGL = gw.all_world2pix(coord_3FGL[:,0], coord_3FGL[:,1], 1)
                  for i in range(0,np.shape(sources_3FGL)[0]):
                     if x_3FGL[i]>=max(0,gxtar-bound*max(ell_a,ell_b)) and x_3FGL[i]<=min(np.shape(gdata)[1]-1,gxtar+bound*max(ell_a,ell_b)) and y_3FGL[i]>=max(0,gytar-bound*max(ell_a,ell_b)) and y_3FGL[i]<=min(np.shape(gdata)[0]-1,gytar+bound*max(ell_a,ell_b)):
                        ax1.text(x_3FGL[i],y_3FGL[i],sources_3FGL[i].replace("L ","L"),horizontalalignment='center',verticalalignment='center',color='yellow',fontsize=7)
                        a = Ellipse(xy=(x_3FGL[i], y_3FGL[i]), width=2*3600./pixscale*coord_3FGL[i,2], height=2*3600./pixscale*coord_3FGL[i,3], angle=coord_3FGL[i,4])
                        a.set_facecolor('none')
                        a.set_edgecolor('yellow')
                        ax1.add_artist(a)
               #if np.size(sources_4FGLe)!=0:
               #   x_4FGLe, y_4FGLe = gw.all_world2pix(np.asarray(sources_4FGLe)[:,1], np.asarray(sources_4FGLe)[:,2], 1)
               #   for i in range(0,np.shape(sources_4FGLe)[0]):
               #      if x_4FGLe[i]>=max(0,gxtar-bound*max(ell_a,ell_b)) and x_4FGLe[i]<=min(np.shape(gdata)[1]-1,gxtar+bound*max(ell_a,ell_b)) and y_4FGLe[i]>=max(0,gytar-bound*max(ell_a,ell_b)) and y_4FGLe[i]<=min(np.shape(gdata)[0]-1,gytar+bound*max(ell_a,ell_b)):
               #         ax1.text(x_4FGLe[i],y_4FGLe[i],sources_4FGLe[i][0].replace("L ","L"),horizontalalignment='center',verticalalignment='center',color='cyan',fontsize=7)
               #         a = Ellipse(xy=(x_4FGLe[i], y_4FGLe[i]), width=2*3600./pixscale*np.asarray(sources_4FGLe)[i,3], height=2*3600./pixscale*np.asarray(sources_4FGLe)[i,4], angle=np.asarray(sources_4FGLe)[i,5])
               #         a.set_facecolor('none')
               #         a.set_edgecolor('cyan')
               #         ax1.add_artist(a)
            for i in range(0,len(gcomps)):
                        #[xcomp, ycomp] = np.dot(Rinv,[gcomps[i,0]-raref, gcomps[i,1]-decref])
                        #[xcomp, ycomp] = [xcomp+xref, ycomp+yref]
                        xcomp, ycomp = gw.all_world2pix(gcomps[i,0], gcomps[i,1], 1)
                        a = Circle(xy=(xcomp, ycomp), radius=20.)
                        a.set_linewidth(0.8)
                        a.set_facecolor('none')
                        a.set_edgecolor('red')
                        ax1.add_artist(a)

#rfilter comparison plot
if np.shape(rcomps)[0] > 1:
            ax2.set_title(r''+FoVname+' FoV bkg subtracted, r band, '+str(len(rcomps))+' comparison stars with $<\sigma>$=%.5f'%rsigma+' mag', fontsize=12)
            ax2.imshow(rdata, interpolation='nearest', cmap='gray', vmin=mr-0.15*sr, vmax=mr+0.15*sr, origin='lower')
            ax2.text(rxtar,rytar,FoVname,horizontalalignment='center',verticalalignment='center',color='y',fontsize=7)
            if search == 'hunt':
               target = Ellipse(xy=(rxtar, rytar), width=2*ell_a, height=2*ell_b, angle=ell_theta)
            else:
               target = Circle(xy=(rxtar,rytar), radius=20)
               target.set_linewidth(1.0)
            target.set_facecolor('none')
            target.set_edgecolor('cyan')
            ax2.add_artist(target)
            if search == 'hunt':
               if np.size(sources_3FGL)!=0:
                  x_3FGL, y_3FGL = rw.all_world2pix(coord_3FGL[:,0], coord_3FGL[:,1], 1)
                  for i in range(0,np.shape(sources_3FGL)[0]):
                     if x_3FGL[i]>=max(0,rxtar-bound*max(ell_a,ell_b)) and x_3FGL[i]<=min(np.shape(rdata)[1]-1,rxtar+bound*max(ell_a,ell_b)) and y_3FGL[i]>=max(0,gytar-bound*max(ell_a,ell_b)) and y_3FGL[i]<=min(np.shape(gdata)[0]-1,gytar+bound*max(ell_a,ell_b)):
                        ax2.text(x_3FGL[i],y_3FGL[i],sources_3FGL[i].replace("L ","L"),horizontalalignment='center',verticalalignment='center',color='yellow',fontsize=7)
                        a = Ellipse(xy=(x_3FGL[i], y_3FGL[i]), width=2*3600./pixscale*coord_3FGL[i,2], height=2*3600./pixscale*coord_3FGL[i,3], angle=coord_3FGL[i,4])
                        a.set_facecolor('none')
                        a.set_edgecolor('yellow')
                        ax2.add_artist(a)
               #if np.size(sources_4FGLe)!=0:
               #   x_4FGLe, y_4FGLe = rw.all_world2pix(np.asarray(sources_4FGLe)[:,1], np.asarray(sources_4FGLe)[:,2], 1)
               #   for i in range(0,np.shape(sources_4FGLe)[0]):
               #      if x_4FGLe[i]>=max(0,rxtar-bound*max(ell_a,ell_b)) and x_4FGLe[i]<=min(np.shape(rdata)[1]-1,rxtar+bound*max(ell_a,ell_b)) and y_4FGLe[i]>=max(0,rytar-bound*max(ell_a,ell_b)) and y_4FGLe[i]<=min(np.shape(rdata)[0]-1,rytar+bound*max(ell_a,ell_b)):
               #         ax2.text(x_4FGLe[i],y_4FGLe[i],sources_4FGLe[i][0].replace("L ","L"),horizontalalignment='center',verticalalignment='center',color='cyan',fontsize=7)
               #         a = Ellipse(xy=(x_4FGLe[i], y_4FGLe[i]), width=2*3600./pixscale*np.asarray(sources_4FGLe)[i,3], height=2*3600./pixscale*np.asarray(sources_4FGLe)[i,4], angle=np.asarray(sources_4FGLe)[i,5])
               #         a.set_facecolor('none')
               #         a.set_edgecolor('cyan')
               #         ax2.add_artist(a)
            for i in range(0,len(rcomps)):
                        xcomp, ycomp = rw.all_world2pix(rcomps[i,0], rcomps[i,1], 1)
                        a = Circle(xy=(xcomp, ycomp), radius=20.)
                        a.set_linewidth(0.8)
                        a.set_facecolor('none')
                        a.set_edgecolor('red')
                        ax2.add_artist(a)
                        a = Circle(xy=(xcomp, ycomp), radius=20.)
                        a.set_linewidth(0.8)
                        a.set_facecolor('none')
                        a.set_edgecolor('red')
                        ax2.add_artist(a)

#ifilter comparison plot
if np.shape(icomps)[0] > 1:
            ax3.set_title(r''+FoVname+' FoV bkg subtracted, i band, '+str(len(icomps))+' comparison stars with $<\sigma>$=%.5f'%isigma+' mag', fontsize=12)
            ax3.imshow(idata, interpolation='nearest', cmap='gray', vmin=mi-0.15*si, vmax=mi+0.15*si, origin='lower')
            ax3.text(ixtar,iytar,FoVname,horizontalalignment='center',verticalalignment='center',color='cyan',fontsize=7)
            if search == 'hunt':
               target = Ellipse(xy=(ixtar, iytar), width=2*ell_a, height=2*ell_b, angle=ell_theta)
            else:
               target = Circle(xy=(ixtar,iytar), radius=20)
               target.set_linewidth(1.0)
            target.set_facecolor('none')
            target.set_edgecolor('cyan')
            ax3.add_artist(target)
            if search == 'hunt':
               if np.size(sources_3FGL)!=0:
                  x_3FGL, y_3FGL = iw.all_world2pix(coord_3FGL[:,0], coord_3FGL[:,1], 1)
                  for i in range(0,np.shape(sources_3FGL)[0]):
                     if x_3FGL[i]>=max(0,ixtar-bound*max(ell_a,ell_b)) and x_3FGL[i]<=min(np.shape(idata)[1]-1,ixtar+bound*max(ell_a,ell_b)) and y_3FGL[i]>=max(0,iytar-bound*max(ell_a,ell_b)) and y_3FGL[i]<=min(np.shape(idata)[0]-1,iytar+bound*max(ell_a,ell_b)):
                        ax3.text(x_3FGL[i],y_3FGL[i],sources_3FGL[i].replace("L ","L"),horizontalalignment='center',verticalalignment='center',color='yellow',fontsize=7)
                        a = Ellipse(xy=(x_3FGL[i], y_3FGL[i]), width=2*3600./pixscale*coord_3FGL[i,2], height=2*3600./pixscale*coord_3FGL[i,3], angle=coord_3FGL[i,4])
                        a.set_facecolor('none')
                        a.set_edgecolor('yellow')
                        ax3.add_artist(a)
               #if np.size(sources_4FGLe)!=0:
               #   x_4FGLe, y_4FGLe = iw.all_world2pix(np.asarray(sources_4FGLe)[:,1], np.asarray(sources_4FGLe)[:,2], 1)
               #   for i in range(0,np.shape(sources_4FGLe)[0]):
               #      if x_4FGLe[i]>=max(0,ixtar-bound*max(ell_a,ell_b)) and x_4FGLe[i]<=min(np.shape(idata)[1]-1,ixtar+bound*max(ell_a,ell_b)) and y_4FGLe[i]>=max(0,iytar-bound*max(ell_a,ell_b)) and y_4FGLe[i]<=min(np.shape(idata)[0]-1,iytar+bound*max(ell_a,ell_b)):
               #         ax3.text(x_4FGLe[i],y_4FGLe[i],sources_4FGLe[i][0].replace("L ","L"),horizontalalignment='center',verticalalignment='center',color='cyan',fontsize=7)
               #         a = Ellipse(xy=(x_4FGLe[i], y_4FGLe[i]), width=2*3600./pixscale*np.asarray(sources_4FGLe)[i,3], height=2*3600./pixscale*np.asarray(sources_4FGLe)[i,4], angle=np.asarray(sources_4FGLe)[i,5])
               #         a.set_facecolor('none')
               #         a.set_edgecolor('cyan')
               #         ax3.add_artist(a)
            for i in range(0,len(icomps)):
                        xcomp, ycomp = iw.all_world2pix(icomps[i,0], icomps[i,1], 1)
                        a = Circle(xy=(xcomp, ycomp), radius=20.)
                        a.set_linewidth(0.8)
                        a.set_facecolor('none')
                        a.set_edgecolor('red')
                        ax3.add_artist(a)

if analysis == 'INT':
            plt.savefig(path+'/'+((FoVname.replace(" ","")).replace("-","m")).replace("+","p")+'_compFoV_cam'+cam+'.png')
else:
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
                                                                        if (mjd[k]-mjd[k-1]) > 0.7:
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
                                                                        if (mjd[k]-mjd[k-1]) > 0.7:
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
                                                                        if (mjd[k]-mjd[k-1]) > 0.7:
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
                                                            plt.ylabel(r"Differential i magnitude",fontsize=8)
                                                if np.size(mjd_single[i])!=0:
                                                   plt.xlabel(r"Time (MJD-%d)"%np.trunc(mjd_single[i][0]),fontsize=8)
                                                   plt.plot(mjd_single[i]-np.trunc(mjd_single[i][0]), dm_single[i], c=c, marker='o', markersize=2, linestyle='None')
                                                   plt.errorbar(mjd_single[i]-np.trunc(mjd_single[i][0]), dm_single[i], sigma_dm_single[i], c=c, marker='o', markersize=2, linestyle='None')

            if analysis == 'INT':
                        plt.savefig(path+'/'+((FoVname.replace(" ","")).replace("-","m")).replace("+","p")+'_complcs_night'+str(i+1)+'_cam'+cam+'.png')
            else:
                        plt.savefig(path+'/'+((FoVname.replace(" ","")).replace("-","m")).replace("+","p")+'_complcs_night'+str(i+1)+'.png')
            plt.close(fig)

