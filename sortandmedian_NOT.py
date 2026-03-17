#!/bin/python

import numpy as np
from astropy.io import fits
import os
import ccdproc
from astropy.nddata import CCDData
from astropy.time import Time
from ccdproc import Combiner, wcs_project

path = os.getcwd()

os.system('osclean_gri.sh')

if os.path.exists(path+'/median_g.fits'):
    os.remove(path+'/median_g.fits')
if os.path.exists(path+'/median_r.fits'):
    os.remove(path+'/median_r.fits')
if os.path.exists(path+'/median_i.fits'):
    os.remove(path+'/median_i.fits')

#gxmin = 0
#rxmin = 0
#ixmin = 0
#zxmin = 0
#xmax = 4058
#gymin = 0
#rymin = 0
#iymin = 0
#zymin = 0
#ymax = 4063
    
images = [f for f in os.listdir(path) if f.endswith('.fits')]
gimages = []
gtimes = []
gdates = []
rimages = []
rtimes = []
rdates = []
iimages = []
itimes = []
idates = []

for i in range(0,(np.size(images))):
    image = fits.open(images[i])
    hdr = image[0].header
    filter = hdr['FAFLTNM'][0]
    if hdr['DATE-OBS'].split('T')[1][0]=='0':
        date = hdr['DATE-OBS'].split('T')[0][:-1]+str(int(hdr['DATE-OBS'].split('T')[0][-1])-1)
        if date[-1]=='0':
            date = date[:-2]+'31'
    else:
        date = hdr['DATE-OBS'].split('T')[0]
    time = Time(hdr['DATE-OBS'], format='isot', scale='utc').mjd
    if filter=='g':
        gimages.append(images[i])
        gtimes.append(time)
        gdates.append(date)
    elif filter=='r':
        rimages.append(images[i])
        rtimes.append(time)
        rdates.append(date)
    elif filter=='i':
        iimages.append(images[i])
        itimes.append(time)
        idates.append(date)
    image.close()

gimages = [gimages[i] for i in np.argsort(gtimes)]
gdates = [gdates[i] for i in np.argsort(gtimes)]
rimages = [rimages[i] for i in np.argsort(rtimes)]
rdates = [rdates[i] for i in np.argsort(rtimes)]
iimages = [iimages[i] for i in np.argsort(itimes)]
idates = [idates[i] for i in np.argsort(itimes)]
np.savetxt(path+'/gfilter.txt',gimages,fmt="%s")
np.savetxt(path+'/rfilter.txt',rimages,fmt="%s")
np.savetxt(path+'/ifilter.txt',iimages,fmt="%s")

#Save the images taken during the night with the best sampling
nights = []
for i in range(0,(np.size(rdates))):
    if i==0:
        nights.append(rdates[i])
    else:
        if rdates[i]!=rdates[i-1]:
            nights.append(rdates[i])
counter = np.zeros(np.size(nights))
for i in range(0,(np.size(nights))):
    for j in range(0,(np.size(rdates))):
        if rdates[j]==nights[i]:
            counter[i]+=1
best_night = nights[np.argmax(counter)]
gimages_best = [gimages[i] for i in range(0,np.size(gimages)) if gdates[i]==best_night]
rimages_best = [rimages[i] for i in range(0,np.size(rimages)) if rdates[i]==best_night]
iimages_best = [iimages[i] for i in range(0,np.size(iimages)) if idates[i]==best_night]
np.savetxt(path+'/gfilter_best.txt',gimages_best,fmt="%s")
np.savetxt(path+'/rfilter_best.txt',rimages_best,fmt="%s")
np.savetxt(path+'/ifilter_best.txt',iimages_best,fmt="%s")

#Separate images in different lists associated to different nights

gimages_single = []
rimages_single = []
iimages_single = []
for i in range(0,(np.size(nights))):
    gimages_single.append([gimages[j] for j in range(0,np.size(gimages)) if gdates[j]==nights[i]])
    rimages_single.append([rimages[j] for j in range(0,np.size(rimages)) if rdates[j]==nights[i]])
    iimages_single.append([iimages[j] for j in range(0,np.size(iimages)) if idates[j]==nights[i]])

#Eventually perform the median excluding a night with bad quality data
gimages_median = gimages_single[0]
rimages_median = rimages_single[0]
iimages_median = iimages_single[0]

#Upper limit on number of images to combine, to avoid memory problems
n_max = 35
            
#median in gfilter

print("##### G FILTER #######")
reprojected = []
for i in range(0,min((np.size(gimages_median)),n_max)):
    ccd = CCDData.read(gimages_median[i], unit="adu")
    print(ccd)
    try:
        print(ccd.size)
    except AttributeError:
        pass
    #ccd = ccdproc.trim_image(ccd[gxmin:gxmax, gymin:gymax])
    if i==0:
        with fits.open(gimages_median[i]) as image:
            # image = fits.open(gimages_median[i])
            hdr = image[0].header
            image.close()
        ccd0 = ccd
    repro_ccd = wcs_project(ccd,ccd0)
    reprojected.append(repro_ccd)
print("Reprojection done")
combiner = Combiner(reprojected)
del reprojected
print("Combiner initialised")
median = np.array(combiner.median_combine())
del combiner
print("Combiner finished")
hdr['INSTRUME'] = 'median'
hdu = fits.PrimaryHDU(data=median,header=hdr)
hdul = fits.HDUList([hdu])
hdul.writeto(path+'/median_g.fits',overwrite=True)
del hdr, hdu, hdul

#median in rfilter

print("##### R FILTER #######")
reprojected = []
for i in range(0,min((np.size(rimages_median)),n_max)):
    ccd = CCDData.read(rimages_median[i], unit="adu")
    print(ccd)
    try:
        print(ccd.size)
    except AttributeError:
        pass
    #ccd = ccdproc.trim_image(ccd[rxmin:rxmax, rymin:rymax])
    if i==0:
        with fits.open(rimages_median[i]) as image:
            #image = fits.open(rimages_median[i])
            hdr = image[0].header
            image.close()
        ccd0 = ccd
    #print(np.shape(ccd))
    repro_ccd = wcs_project(ccd,ccd0)
    reprojected.append(repro_ccd)
print("Reprojection done")
combiner = Combiner(reprojected)
del reprojected
print("Combiner initialised")
median = np.array(combiner.median_combine())
del combiner
print("Combiner finished")
hdr['INSTRUME'] = 'median'
hdu = fits.PrimaryHDU(data=median,header=hdr)
hdul = fits.HDUList([hdu])
hdul.writeto(path+'/median_r.fits',overwrite=True)
del hdr, hdu, hdul

#median in ifilter

print("##### I FILTER #######")
reprojected = []
for i in range(0,min((np.size(iimages_median)),n_max)):
    ccd = CCDData.read(iimages_median[i], unit="adu")
    print(ccd)
    try:
        print(ccd.size)
    except AttributeError:
        pass
    #ccd = ccdproc.trim_image(ccd[ixmin:ixmax, iymin:iymax])
    if i==0:
        with fits.open(iimages_median[i]) as image:
            #image = fits.open(iimages_median[i])
            hdr = image[0].header
            image.close()
        ccd0 = ccd
    repro_ccd = wcs_project(ccd,ccd0)
    reprojected.append(repro_ccd)
print("Reprojection done")
combiner = Combiner(reprojected)
del reprojected
print("Combiner initialised")
median = np.array(combiner.median_combine())
del combiner
print("Combiner finished")
hdr['INSTRUME'] = 'median'
hdu = fits.PrimaryHDU(data=median,header=hdr)
hdul = fits.HDUList([hdu])
hdul.writeto(path+'/median_i.fits',overwrite=True)
del hdr, hdu, hdul
        
        
    



   
    

