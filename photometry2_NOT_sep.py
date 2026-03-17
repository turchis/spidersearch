#!/home/gudrun/marcotu/python-envs/astrosource/bin/python

import logging
from urllib.parse import urljoin

import numpy as np
import matplotlib.pyplot as plt
import os

os.system('source ~/python-envs/astrosource/bin/activate')

import shutil
from astropy.table import Table
from astropy.io import fits
from astropy import wcs
from astropy.utils.data import get_pkg_data_filename
import sep
from requests import HTTPError

#from banzai.utils import stats, array_utils
#from banzai.utils.photometry_utils import get_reference_sources, match_catalogs, to_magnitude, fit_photometry
#from banzai.data import DataTable

#from skimage import measure

sep.set_sub_object_limit(int(1e5))

#def radius_of_contour(contour, source):
#    x = contour[:, 1]
#    y = contour[:, 0]
#    x_center = (source['xmax'] - source['xmin'] + 1) / 2.0 - 0.5
#    y_center = (source['ymax'] - source['ymin'] + 1) / 2.0 - 0.5

#    return np.percentile(np.sqrt((x - x_center)**2.0 + (y - y_center)** 2.0), 90)

def photometry(image,image_file,threshold,min_area,star_scale,cfwhm,min_fwhm,max_fwhm):
            # Increase the internal buffer size in sep. This is most necessary for crowded fields.
            ny, nx = image[0].data.shape
            #print(nx,ny)
            
            sep.set_extract_pixstack(int(nx * ny - 1))
            data = image[0].data
            hdr = image[0].header

            # Fits can be backwards byte order, so fix that if need be and subtract
            # the background
            try:
                bkg = sep.Background(data, bw=4, bh=4, fw=3, fh=3)
            except ValueError:
                data = data.byteswap(True).newbyteorder()
                bkg = sep.Background(data, bw=4, bh=4, fw=3, fh=3)
            bkg.subfrom(data)
            bkg_rms = bkg.rms()
            #hdu = fits.PrimaryHDU(data)
            #hdul = fits.HDUList([hdu])
            #hdul.writeto('median_r_after.fits',overwrite=True)
            #plt.imsave(path+'/bkg_iband.png',bkg.back(), cmap='gray', origin='lower')
            #print(bkg.globalrms)
            #bkg_globalrms = bkg.globalrms

            # Do an initial source detection
            sources = sep.extract(data, threshold, minarea=min_area, deblend_nthresh=32, deblend_cont=0.005, err=bkg_rms)
           # print(len(sources))

            # Convert the detections into a table
            sources = Table(sources)

            # We remove anything with a detection flag >= 8
            # This includes memory overflows and objects that are too close the edge
            sources = sources[sources['flag'] < 8]

            sources = prune_nans_from_table(sources)
            #print(len(sources))

            # Calculate the ellipticity
            sources['ellipticity'] = 1.0 - (sources['b'] / sources['a'])

            # Fix any value of theta that are invalid due to floating point rounding
            # -pi / 2 < theta < pi / 2
            sources['theta'][sources['theta'] > (np.pi / 2.0)] -= np.pi
            sources['theta'][sources['theta'] < (-np.pi / 2.0)] += np.pi

            # Calculate the kron radius
            #kronrad, krflag = sep.kron_radius(data, sources['x'], sources['y'],
                                              #sources['a'], sources['b'],
                                              #sources['theta'], 6.0)
            #sources['flag'] |= krflag
            #sources['kronrad'] = kronrad

            # Calcuate the equivilent of flux_auto
            #flux, fluxerr, flag = sep.sum_ellipse(data, sources['x'], sources['y'],
                                                  #sources['a'], sources['b'],
                                                  #np.pi / 2.0, 2.5 * kronrad,
                                                  #subpix=1, err=bkg_rms)

             # Calculate flux profiles at 50% and FWHMs of the stars and exclude the stars with fwhm > max_fwhm
            fwhm, flag = sep.flux_radius(data, sources['x'], sources['y'],
                                               6.0 * sources['a'], 0.5,
                                               normflux=sources['flux'], subpix=5)
            fwhm = 2 * fwhm
            sources['flag'] |= flag
            sources['targfwhm'] = fwhm
            #sources = sources[sources['targfwhm'] >= min_fwhm]
            sources = sources[sources['targfwhm'] <= max_fwhm]
            #print(len(sources))
            #sources['fwhm'] = 2 * np.sqrt(np.log(2) * (sources['a']**2 + sources['b']**2))
            #sources = sources[sources['fwhm'] <= max_fwhm]
            #print(len(sources))
            #sources['fwhm'] = np.nan
            #sources['fwtm'] = np.nan
            # Here we estimate contours
            #for source in sources:
            #    if source['flag'] == 0:
            #        for ratio, keyword in zip([0.5, 0.1], ['fwhm', 'fwtm']):
            #            contours = measure.find_contours(data[source['ymin']: source['ymax'] + 1,
            #                                             source['xmin']: source['xmax'] + 1],
            #                                             ratio * source['peak'])
            #            if contours:
                            # If there are multiple contours like a donut might have take the outer
            #                contour_radii = [radius_of_contour(contour, source) for contour in contours]
            #                source[keyword] = 2.0 * np.nanmax(contour_radii)
            #sources = sources[sources['fwhm'] > 0.]
            #sources = sources[sources['fwhm'] <= max_fwhm]
            #print(len(sources))
            # Put in the catalogue a more accurate estimate of the FWHM obtained by averaging FWHMs of the comparison stars
            sources['fwhm'] = cfwhm*np.ones(len(sources))

            # Set apertures radii r=star_scale*fwhm
            aprad = star_scale*cfwhm*np.ones(len(sources))
            sources['aprad'] = aprad
                                                  
            # Do circular aperture photometry
            #(from ALFOSC documentation)
            pixel_scale = 0.2138*2
            gain = float(image[0].header['GAIN'])
            flux, fluxerr, flag = sep.sum_circle(data, sources['x'], sources['y'],
                                                     aprad, gain=gain, err=bkg_rms)
            sources['flux'] = flux
            sources['fluxerr'] = fluxerr
            sources['flag'] |= flag
            
            

            # Do circular aperture photometry for diameters of 1" to 6"
            #pixel_scale = float(image[0].header['PIXSCALE'])
            #gain = float(image[0].header['AVGAIN'])
            #for diameter in [1, 2, 3, 4, 5, 6]:
            #    flux, fluxerr, flag = sep.sum_circle(data, sources['x'], sources['y'],
            #                                         diameter / 2.0 / pixel_scale, gain=gain, err=bkg_rms)
            #    sources['fluxaper{0}'.format(diameter)] = flux
            #    sources['fluxerr{0}'.format(diameter)] = fluxerr
            #    sources['flag'] |= flag

            # Measure the flux profile at 50%
            fluxrad50, flag = sep.flux_radius(data, sources['x'], sources['y'],
                                               6.0 * sources['a'], 0.5,
                                               normflux=sources['flux'], subpix=5)
            sources['flag'] |= flag
            sources['fluxrad50'] = fluxrad50
            #flux_radii, flag = sep.flux_radius(data, sources['x'], sources['y'],
            #                                   6.0 * sources['a'], [0.25, 0.5, 0.75],
            #                                   normflux=sources['flux'], subpix=5)
            #sources['flag'] |= flag
            #sources['fluxrad25'] = flux_radii[:, 0]
            #sources['fluxrad50'] = flux_radii[:, 1]
            #sources['fluxrad75'] = flux_radii[:, 2]

            # Cut individual bright pixels. Often cosmic rays
            sources = sources[sources['fluxrad50'] > 0.5]
            #print(len(sources))


            # Calculate the windowed positions
            sig = 2.0 / 2.35 * sources['fwhm']
            xwin, ywin, flag = sep.winpos(data, sources['x'], sources['y'], sig)
            #sources['flag'] |= flag
            sources['xwin'] = xwin
            sources['ywin'] = ywin

            # Calculate the average background at each source
            bkgflux, fluxerr, flag = sep.sum_circle(bkg.back(), sources['x'], sources['y'],
                                                    sources['aprad'], subpix=1)
            # masksum, fluxerr, flag = sep.sum_ellipse(mask, sources['x'], sources['y'],
            #                                         sources['a'], sources['b'], np.pi / 2.0,
            #                                         2.5 * kronrad, subpix=1)

            background_area = (sources['aprad']**2) * np.pi # - masksum
            sources['background'] = bkgflux
            sources['background'][background_area > 0] /= background_area[background_area > 0]
            # Update the catalog to match fits convention instead of python array convention
            sources['x'] += 1.0
            sources['y'] += 1.0

            sources['xpeak'] += 1
            sources['ypeak'] += 1

            sources['xwin'] += 1.0
            sources['ywin'] += 1.0

            sources['theta'] = np.degrees(sources['theta'])

            catalog = sources['x', 'y', 'xwin', 'ywin', 'xpeak', 'ypeak',
                              'flux', 'fluxerr', 'peak', 'background', 'targfwhm', 'fwhm',
                              'a', 'b', 'theta', 'aprad', 'ellipticity',
                              'fluxrad50', 'x2', 'y2', 'xy', 'flag']

            # Add the units and description to the catalogs
            catalog['x'].unit = 'pixel'
            catalog['x'].description = 'X coordinate of the object'
            catalog['y'].unit = 'pixel'
            catalog['y'].description = 'Y coordinate of the object'
            catalog['xwin'].unit = 'pixel'
            catalog['xwin'].description = 'Windowed X coordinate of the object'
            catalog['ywin'].unit = 'pixel'
            catalog['ywin'].description = 'Windowed Y coordinate of the object'
            catalog['xpeak'].unit = 'pixel'
            catalog['xpeak'].description = 'X coordinate of the peak'
            catalog['ypeak'].unit = 'pixel'
            catalog['ypeak'].description = 'Windowed Y coordinate of the peak'
            catalog['flux'].unit = 'count'
            catalog['flux'].description = 'Flux within the circular aperture'
            catalog['fluxerr'].unit = 'count'
            catalog['fluxerr'].description = 'Error on the flux within circular aperture'
            catalog['peak'].unit = 'count'
            catalog['peak'].description = 'Peak flux (flux at xpeak, ypeak)'
            #for diameter in [1, 2, 3, 4, 5, 6]:
            #    catalog['fluxaper{0}'.format(diameter)].unit = 'count'
            #    catalog['fluxaper{0}'.format(diameter)].description = 'Flux from fixed circular aperture: {0}" diameter'.format(diameter)
            #    catalog['fluxerr{0}'.format(diameter)].unit = 'count'
            #    catalog['fluxerr{0}'.format(diameter)].description = 'Error on Flux from circular aperture: {0}"'.format(diameter)

            catalog['background'].unit = 'count'
            catalog['background'].description = 'Average background value in the aperture'
            catalog['targfwhm'].unit = 'pixel'
            catalog['targfwhm'].description = 'FWHM of the object'
            catalog['fwhm'].unit = 'pixel'
            catalog['fwhm'].description = 'Average FWHM of the comparison stars'
            #catalog['fwtm'].unit = 'pixel'
            #catalog['fwtm'].description = 'Full-Width Tenth Maximum'
            catalog['a'].unit = 'pixel'
            catalog['a'].description = 'Semi-major axis of the object'
            catalog['b'].unit = 'pixel'
            catalog['b'].description = 'Semi-minor axis of the object'
            catalog['theta'].unit = 'degree'
            catalog['theta'].description = 'Position angle of the object'
            catalog['aprad'].unit = 'pixel'
            catalog['aprad'].description = 'Aperture radius used for extraction'
            #catalog['kronrad'].unit = 'pixel'
            #catalog['kronrad'].description = 'Kron radius used for extraction'
            catalog['ellipticity'].description = 'Ellipticity'
            #catalog['fluxrad25'].unit = 'pixel'
            #catalog['fluxrad25'].description = 'Radius containing 25% of the flux'
            catalog['fluxrad50'].unit = 'pixel'
            catalog['fluxrad50'].description = 'Radius containing 50% of the flux'
            #catalog['fluxrad75'].unit = 'pixel'
            #catalog['fluxrad75'].description = 'Radius containing 75% of the flux'
            catalog['x2'].unit = 'pixel^2'
            catalog['x2'].description = 'Variance on X coordinate of the object'
            catalog['y2'].unit = 'pixel^2'
            catalog['y2'].description = 'Variance on Y coordinate of the object'
            catalog['xy'].unit = 'pixel^2'
            catalog['xy'].description = 'XY covariance of the object'
            catalog['flag'].description = 'Bit mask of extraction/photometry flags'

            catalog.sort('flux')
            catalog.reverse()

            # Save some background statistics in the header
            #mean_background = stats.sigma_clipped_mean(bkg.back(), 5.0)
            #image[0].header['L1MEAN'] = (mean_background,
            #                        '[counts] Sigma clipped mean of frame background')

            median_background = np.median(bkg.back())
            image[0].header['L1MEDIAN'] = (median_background,
                                      '[counts] Median of frame background')

            #std_background = stats.robust_standard_deviation(bkg.back())
            #image[0].header['L1sigma'] = (std_background,
            #                         '[counts] Robust std dev of frame background')

            # Save some image statistics to the header
            good_objects = catalog['flag'] == 0
            for quantity in ['fwhm', 'ellipticity', 'theta']:
                good_objects = np.logical_and(good_objects, np.logical_not(np.isnan(catalog[quantity])))
            if good_objects.sum() == 0:
                image[0].header['L1FWHM'] = ('NaN', '[arcsec] Frame FWHM in arcsec')
                image[0].header['L1FWTM'] = ('NaN', 'Ratio of FWHM to Full-Width Tenth Max')

                image[0].header['L1ELLIP'] = ('NaN', 'Mean image ellipticity (1-B/A)')
                image[0].header['L1ELLIPA'] = ('NaN', '[deg] PA of mean image ellipticity')
            else:
                seeing = np.nanmedian(catalog['fwhm'][good_objects]) * pixel_scale
                image[0].header['L1FWHM'] = (seeing, '[arcsec] Frame FWHM in arcsec')
                #image[0].header['L1FWTM'] = (np.nanmedian(catalog['fwtm'][good_objects] / catalog['fwhm'][good_objects]),
                #                        'Ratio of FWHM to Full-Width Tenth Max')

                #mean_ellipticity = stats.sigma_clipped_mean(catalog['ellipticity'][good_objects], 3.0)
                #image[0].header['L1ELLIP'] = (mean_ellipticity, 'Mean image ellipticity (1-B/A)')

                #mean_position_angle = stats.sigma_clipped_mean(catalog['theta'][good_objects], 3.0)
                #image[0].header['L1ELLIPA'] = (mean_position_angle,'[deg] PA of mean image ellipticity')

            if len(image)>1:
                fits.update(image_file,np.array(catalog),1)
            else:
                fits.append(image_file,np.array(catalog))
            with fits.open(image_file, mode='update') as filehandle:
                        filehandle[0].header['EXTNAME'] = 'SCI'
                        filehandle[1].header['EXTNAME'] = 'CAT'

            w = wcs.WCS(hdr)
            ra, dec = w.all_pix2world(sources['x'], sources['y'], 1)
            regr = open(image_file.split('fits')[0]+'reg','w')
            print("# Region file format: DS9 version 4.1",file=regr)
            print("global dashlist=8 3 width=1 font=\"helvetica 10 normal roman\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1",file=regr)
            print("fk5",file=regr)
            for i in range(0,len(sources)):
                print("circle(%f,%f,%f\") #  color=green width=1" % (ra[i],dec[i],aprad[i]*pixel_scale),file=regr)
            regr.close()
                
            return image

def prune_nans_from_table(table):
    nan_in_row = np.zeros(len(table), dtype=bool)
    for col in table.colnames:
        nan_in_row |= np.isnan(table[col])
    return table[~nan_in_row]

path = os.getcwd()

gimages = np.genfromtxt(path+'/gfilter.txt',dtype='str')
gimages = gimages.tolist()
ng = len(np.genfromtxt(path+'/gfilter/usedImages.txt',dtype='str'))-1
gcomps = np.atleast_2d(np.genfromtxt(path+'/gfilter/compsUsed.csv', delimiter=','))
rimages = np.genfromtxt(path+'/rfilter.txt',dtype='str')
rimages = rimages.tolist()
nr = len(np.genfromtxt(path+'/rfilter/usedImages.txt',dtype='str'))-1
rcomps = np.atleast_2d(np.genfromtxt(path+'/rfilter/compsUsed.csv', delimiter=','))
iimages = np.genfromtxt(path+'/ifilter.txt',dtype='str')
iimages = iimages.tolist()
ni = len(np.genfromtxt(path+'/ifilter/usedImages.txt',dtype='str'))-1
icomps = np.atleast_2d(np.genfromtxt(path+'/ifilter/compsUsed.csv', delimiter=','))

#getting fwhm of the comparison stars for each frame and average bewtween them

if np.shape(gcomps)[0]>1:
            fwhm = np.zeros((len(gcomps),ng))
            for i in range(0,len(gcomps)):
                        cfwhm = np.genfromtxt(path+'/gfilter/outputcats/doerPhot_V'+str(i+1)+'.csv', delimiter=',',usecols=(6),unpack=True)
                        for j in range(0,len(cfwhm)):
                                    fwhm[i,j] = cfwhm[j]
            gfwhm_avg = np.zeros(ng)
            for i in range(0,ng):
                        gfwhm_avg[i] = np.mean(fwhm[:,i])

if np.shape(rcomps)[0]>1:
            fwhm = np.zeros((len(rcomps),nr))
            for i in range(0,len(rcomps)):
                        cfwhm = np.genfromtxt(path+'/rfilter/outputcats/doerPhot_V'+str(i+1)+'.csv', delimiter=',',usecols=(6),unpack=True)
                        for j in range(0,len(cfwhm)):
                                    fwhm[i,j] = cfwhm[j]
            rfwhm_avg = np.zeros(nr)
            for i in range(0,nr):
                        rfwhm_avg[i] = np.mean(fwhm[:,i])
if np.shape(icomps)[0]>1:
            fwhm = np.zeros((len(icomps),ni))
            for i in range(0,len(icomps)):
                        cfwhm = np.genfromtxt(path+'/ifilter/outputcats/doerPhot_V'+str(i+1)+'.csv', delimiter=',',usecols=(6),unpack=True)
                        for j in range(0,len(cfwhm)):
                                    fwhm[i,j] = cfwhm[j]
            ifwhm_avg = np.zeros(ni)
            for i in range(0,ni):
                        ifwhm_avg[i] = np.mean(fwhm[:,i])
   
#Gets in input the factors for the extraction radius and the detection threshold 
fext = float(input("Enter the scale factor for the FWHM used for the extraction radius in the 1st source detection (default value 1.2):\n") or 1.2)
while fext<=0.:
            fext = float(input("Invalid value, please insert a positive extraction factor:\n") or 1.2)
s = float(input("Enter the factor s for the detection threshold s*bkg_rms (default value 2.0):\n") or 2.0)
while s<=0.:
            s = float(input("Invalid value, please insert a positive extraction factor:\n"))

if np.shape(gcomps)[0]>1:
            image_file = path+'/gfilter/median_g.fits'
            image = fits.open(path+'/gfilter/median_g.fits')
            gfwhm = np.mean(gfwhm_avg)
            photometry(image,image_file=image_file,threshold=s,min_area=2,star_scale=fext,cfwhm=gfwhm,min_fwhm=2,max_fwhm=5)
            print('median_g.fits')
            image.close()
            for i in range(0,(len(gimages))):
                        image_file = path+'/gfilter/'+gimages[i]  
                        image = fits.open(path+'/gfilter/'+gimages[i])
                        if gfwhm_avg[min([i,ng-1])]==0.:
                                    gfwhm_avg[min([i,ng-1])] = gfwhm
                        photometry(image,image_file=image_file,threshold=s,min_area=2,star_scale=fext,cfwhm=gfwhm_avg[min([i,ng-1])],min_fwhm=2,max_fwhm=5)
                        print(gimages[i])
                        image.close()
                        print("image should now be closed")

if np.shape(rcomps)[0]>1:
            image_file = path+'/rfilter/median_r.fits'
            image = fits.open(path+'/rfilter/median_r.fits')
            rfwhm = np.mean(rfwhm_avg)
            photometry(image,image_file=image_file,threshold=s,min_area=2,star_scale=fext,cfwhm=rfwhm,min_fwhm=2,max_fwhm=5)
            print('median_r.fits')
            image.close()
            for i in range(0,(len(rimages))):
                        image_file = path+'/rfilter/'+rimages[i]  
                        image = fits.open(path+'/rfilter/'+rimages[i])
                        if rfwhm_avg[min([i,nr-1])]==0.:
                                    rfwhm_avg[min([i,nr-1])] = rfwhm
                        photometry(image,image_file=image_file,threshold=s,min_area=2,star_scale=fext,cfwhm=rfwhm_avg[min([i,nr-1])],min_fwhm=2,max_fwhm=5)
                        print(rimages[i])
                        image.close()
                        print("image should now be closed")

if np.shape(icomps)[0]>1:
            image_file = path+'/ifilter/median_i.fits'  
            image = fits.open(path+'/ifilter/median_i.fits')
            ifwhm = np.mean(ifwhm_avg)
            photometry(image,image_file=image_file,threshold=1.5,min_area=0,star_scale=fext,cfwhm=ifwhm,min_fwhm=0,max_fwhm=30)
            print('median_i.fits')
            image.close()
            for i in range(0,(len(iimages))):
                        image_file = path+'/ifilter/'+iimages[i]  
                        image = fits.open(path+'/ifilter/'+iimages[i])
                        if ifwhm_avg[min([i,ni-1])]==0.:
                                    ifwhm_avg[min([i,ni-1])] = ifwhm
                        photometry(image,image_file=image_file,threshold=1.5,min_area=0,star_scale=fext,cfwhm=ifwhm_avg[min([i,ni-1])],min_fwhm=0,max_fwhm=30)
                        print(iimages[i])
                        image.close()
                        print("image should now be closed")

#image_file = path+'/ifilter/'+iimages[6]  
#image = fits.open(path+'/ifilter/'+iimages[6])
#if ifwhm_avg[min([6,ni-1])]==0.:
#            ifwhm_avg[min([6,ni-1])] = ifwhm
#photometry(image,image_file=image_file,threshold=1.5,min_area=0,star_scale=fext,cfwhm=ifwhm_avg[min([6,ni-1])],min_fwhm=0,max_fwhm=30)
#print(iimages[6])
#image.close()
#print("image should now be closed")

