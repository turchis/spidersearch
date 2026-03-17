#!/bin/python

import numpy as np
from astropy.io import fits
import os
import ccdproc
from astropy.nddata import CCDData
from ccdproc import Combiner, wcs_project

path = os.getcwd()

analysis = input("Enter 'STELLA', 'INT' or 'LCO' for the type of data (default value 'LCO'):\n") or 'LCO'
while analysis!='STELLA' and analysis!='INT' and analysis!='LCO':
   analysis = input("Invalid value, please insert 'STELLA', 'INT' or 'LCO':\n") or 'LCO'

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
#xmax = 4058
#gymin = 0
#rymin = 0
#iymin = 0
#ymax = 4063
    
images = [f for f in os.listdir(path) if f.endswith('.fits')]
gimages = []
gtimes = []
gdates = []
#gdim = []
rimages = []
rtimes = []
rdates = []
#rdim = []
iimages = []
itimes = []
idates = []
#idim = []

if analysis == 'STELLA':
    for i in range(0,(np.size(images))):
        image = fits.open(images[i])
        hdr = image[0].header
        filter = hdr['FILTER']
        date = images[i].split("-")[0].split("science")[1].split("A")[0]
        time = hdr['JD-OBS']
        if filter=='gp':
            gimages.append(images[i])
            gtimes.append(time)
            gdates.append(date)
            #gdim.append(np.shape(data))
        elif filter=='rp':
            rimages.append(images[i])
            rtimes.append(time)
            rdates.append(date)
            #rdim.append(np.shape(data))
        elif filter=='ip':
            iimages.append(images[i])
            itimes.append(time)
            idates.append(date)
            #idim.append(np.shape(data))
        image.close()

elif analysis == 'LCO':
    for i in range(0,(np.size(images))):
        image = fits.open(images[i])
        hdr = image[0].header
        filter = hdr['FILTER']
        date = images[i].split("-")[2]
        time = hdr['MJD-OBS']
        if filter=='gp':
            gimages.append(images[i])
            gtimes.append(time)
            gdates.append(date)
            #gdim.append(np.shape(data))
        elif filter=='rp':
            rimages.append(images[i])
            rtimes.append(time)
            rdates.append(date)
            #rdim.append(np.shape(data))
        elif filter=='ip':
            iimages.append(images[i])
            itimes.append(time)
            idates.append(date)
            #idim.append(np.shape(data))
        image.close()
        
elif analysis == 'INT':
    for i in range(0,(np.size(images))):
        image = fits.open(images[i])
        hdr = image[0].header
        filter = hdr['WFFBAND']
        if hdr['UTOBS'][0]=='0':
            date = hdr['DATE-OBS'][:-1]+str(int(hdr['DATE-OBS'][-1])-1)
            if date[-1]=='0':
                date = date[:-2]+'31'
        else:
            date = hdr['DATE-OBS']
        time = hdr['MJD-OBS']
        if filter=='g':
            gimages.append(images[i])
            gtimes.append(time)
            gdates.append(date)
            #gdim.append(np.shape(data))
        elif filter=='r':
            rimages.append(images[i])
            rtimes.append(time)
            rdates.append(date)
            #rdim.append(np.shape(data))
        elif filter=='i':
            iimages.append(images[i])
            itimes.append(time)
            idates.append(date)
            #idim.append(np.shape(data))
        image.close()
#gymax, gxmax = gdim[0][0], gdim[0][1]
#for i in range(1,len(gdim)):
#    if gdim[i][0] < gymax or gdim[i][1] < gxmax:
#        gymax, gxmax = gdim[i][0], gdim[i][1]
gimages = [gimages[i] for i in np.argsort(gtimes)]
gdates = [gdates[i] for i in np.argsort(gtimes)]
#rymax, rxmax = rdim[0][0], rdim[0][1]
#for i in range(1,len(rdim)):
#    if rdim[i][0] < rymax or rdim[i][1] < rxmax:
#        rymax, rxmax = rdim[i][0], rdim[i][1]
rimages = [rimages[i] for i in np.argsort(rtimes)]
rdates = [rdates[i] for i in np.argsort(rtimes)]
#iymax, ixmax = idim[0][0], idim[0][1]
#for i in range(1,len(idim)):
#    if idim[i][0] < iymax or idim[i][1] < ixmax:
#        iymax, ixmax = idim[i][0], idim[i][1]
iimages = [iimages[i] for i in np.argsort(itimes)]
idates = [idates[i] for i in np.argsort(itimes)]
np.savetxt(path+'/gfilter.txt',gimages,fmt="%s")
np.savetxt(path+'/rfilter.txt',rimages,fmt="%s")
np.savetxt(path+'/ifilter.txt',iimages,fmt="%s")
#os.system("conda activate iraf27; cl | log | conda deactivate ")
#os.system("conda activate iraf27; cl | imcombine @gfilter.txt median_g combine=median offset='wcs' | imcombine @rfilter.txt median_r combine=median offset='wcs' | imcombine @ifilter.txt median_i combine=median offset='wcs' | log")

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
n_max = 30

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

#if(os.path.isdir(path+'/gfilter'))
#            shutil.rmtree(path+'/gfilter')
#os.system('mkdir '+path+'/gfilter')
#os.system('mv median_g.fits '+path+'/gfilter')
#for i in range(0,(len(gimages))):
#            os.system('mv '+gimages[i]+' '+path+'/gfilter')

#if(os.path.isdir(path+'/rfilter')):
#            shutil.rmtree(path+'/rfilter')
#os.system('mkdir '+path+'/rfilter')
#rimages = np.genfromtxt(path+'/rfilter.txt',dtype='str')
#os.system('mv median_r.fits '+path+'/rfilter')
#for i in range(0,(len(rimages))):
#            os.system('mv '+rimages[i]+' '+path+'/rfilter')

#if(os.path.isdir(path+'/ifilter')):
#            shutil.rmtree(path+'/ifilter')
#os.system('mkdir '+path+'/ifilter')
#iimages = np.genfromtxt(path+'/ifilter.txt',dtype='str')
#os.system('mv median_i.fits '+path+'/ifilter')
#for i in range(0,(len(iimages))):
#            os.system('mv '+iimages[i]+' '+path+'/ifilter')
