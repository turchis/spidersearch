#!/bin/python

import numpy as np
import os

#Edit fits files with the correct wcs solution via astrometry.net
path = os.getcwd()

if not os.path.isdir(path+'/wcsbrokenimages'):
   os.mkdir(path+'/wcsbrokenimages')

images = [f for f in os.listdir(path) if f.endswith('.fits')]

for i in range(0,(np.size(images))):
   print('Processing wcs solution of '+images[i]+' via Astrometry.net')
   #Standard wcs solution fixing
   #os.system('solve-field -O -p '+images[i]+' --guess-scale -N '+images[i].split('.fits')[0]+'astr.fits --fits-image')
   #os.system('solve-field -O -p '+images[i]+' --guess-scale -N '+images[i].split('.fits')[0]+'astr.fits --fits-image')
   os.system('solve-field -O -p '+images[i]+' --guess-scale -N '+images[i].split('.fits')[0]+'astr.fits --fits-image')
   #os.system('solve-field -O -p '+images[i]+' --scale-low 4.0 --scale-high 5.0 -u arcsecperpix --downsample 4 -w 4096 -e 4096 -N '+images[i].split('.fits')[0]+'astr.fits --fits-image')
   #For NOT 2x2 binned fits files with bad wcs solution, try this
   #os.system('solve-field -O -p '+images[i]+' -L 0.4274 -H 0.4276 -u arcsecperpix --downsample 4 -w 1074 -e 1051 -N '+images[i].split('.fits')[0]+'astr.fits --fits-image')
   #os.system('solve-field -O -p '+images[i]+' -L 0.4274 -H 0.4276 -N '+images[i].split('.fits')[0]+'astr.fits --fits-image')
   #os.system('solve-field -O -p '+images[i]+' -N '+images[i].split('.fits')[0]+'astr.fits --fits-image')
   os.remove(images[i].split('.fits')[0]+'.axy')
   os.remove(images[i].split('.fits')[0]+'.corr')
   os.remove(images[i].split('.fits')[0]+'-indx.xyls')
   os.remove(images[i].split('.fits')[0]+'.match')
   os.remove(images[i].split('.fits')[0]+'.rdls')
   os.remove(images[i].split('.fits')[0]+'.solved')
   os.remove(images[i].split('.fits')[0]+'.wcs')

   if os.path.exists(images[i].split('.fits')[0]+'astr.fits'):
      os.remove(images[i])
      os.rename(images[i].split('.fits')[0]+'astr.fits',images[i])
   else:
      os.rename(path+'/'+images[i],path+'/wcsbrokenimages/'+images[i])

   
