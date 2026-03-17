#!/bin/python

import numpy as np
import os
import shutil
from astrometry_net_client import Session
from astrometry_net_client import FileUpload

#Edit fits files with the correct wcs solution via astrometry.net
path = os.getcwd()

images = [f for f in os.listdir(path) if f.endswith('.fits')]
s = Session(api_key = 'teztkcltmttgqlal')

if not os.path.isdir(path+'/wcsbrokenimages'):
   os.mkdir(path+'/wcsbrokenimages')

for i in range(0,(np.size(images))):
   print('Processing wcs solution of '+images[i]+' via Astrometry.net')
   upl = FileUpload(images[i], session=s)
   submission = upl.submit()
   submission.until_done()
   job = submission.jobs[0]
   job.until_done()

   if job.success():
       os.remove(images[i])
       print(images[i]+' processed successfully')
       os.system('wget https://nova.astrometry.net/new_fits_file/'+str(job.id))
       os.rename(str(job.id),images[i])
   else:
       print(images[i]+' processing failed')
       os.rename(path+'/'+images[i],path+'/wcsbrokenimages/'+images[i])

   print(job.info())
