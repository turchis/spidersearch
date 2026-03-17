#!/home/gudrun/marcotu/python-envs/astrosource/bin/python

import time
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
import sep
from requests import HTTPError

#from banzai.utils import stats, array_utils
#from banzai.utils.photometry_utils import get_reference_sources, match_catalogs, to_magnitude, fit_photometry
#from banzai.data import DataTable

#from skimage import measure

sep.set_sub_object_limit(int(1e4))

#def radius_of_contour(contour, source):
#    x = contour[:, 1]
#    y = contour[:, 0]
#    x_center = (source['xmax'] - source['xmin'] + 1) / 2.0 - 0.5
#    y_center = (source['ymax'] - source['ymin'] + 1) / 2.0 - 0.5

#    return np.percentile(np.sqrt((x - x_center)**2.0 + (y - y_center)** 2.0), 90)

def photometry(image,image_file,threshold,min_area,star_scale,min_fwhm,max_fwhm,analysis):
            # Increase the internal buffer size in sep. This is most necessary for crowded fields.
            ny, nx = image[0].data.shape
            #print(nx,ny)
            
            sep.set_extract_pixstack(int(nx * ny - 1))
            data = image[0].data
            hdr = image[0].header

            # Fits can be backwards byte order, so fix that if need be and subtract
            # the background
            try:
                bkg = sep.Background(data, bw=8, bh=8, fw=3, fh=3)
            except ValueError:
                data = data.byteswap(True).newbyteorder()
                bkg = sep.Background(data, bw=8, bh=8, fw=3, fh=3)
            bkg.subfrom(data)
            bkg_rms = bkg.rms()
            #hdu = fits.PrimaryHDU(data)
            #hdul = fits.HDUList([hdu])
            #hdul.writeto('median_r_after.fits',overwrite=True)
            #plt.imsave(path+'/bkg_iband.png',bkg.back(), cmap='gray', origin='lower')
            #print(bkg.globalrms)
            #bkg_globalrms = bkg.globalrms

            # Do an initial source detection
            sources = sep.extract(data, threshold, minarea=min_area, deblend_nthresh=32, deblend_cont=0.00001, err=bkg_rms)
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
            sources['fwhm'] = fwhm
            #sources = sources[sources['fwhm'] <= max_fwhm]
            #sources = sources[sources['fwhm'] >= min_fwhm]
            #print(len(sources))
            #sources['fwhm'] = 2 * np.sqrt(np.log(2) * (sources['x2'] + sources['y2']))
            sources = sources[sources['fwhm'] <= max_fwhm]
            sources = sources[sources['fwhm'] >= min_fwhm]
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

            # Set apertures radii r=star_scale*fwhm
            aprad = star_scale*sources['fwhm']
            sources['aprad'] = aprad
                                                  
            # Do circular aperture photometry
            if analysis == 'STELLA':
                        pixel_scale = float(image[0].header['PIXSCALE'])
                        gain = float(image[0].header['AVGAIN'])
            elif analysis == 'LCO':
                        pixel_scale = float(image[0].header['PIXSCALE'])
                        gain = float(image[0].header['GAIN'])
            elif analysis == 'INT':
                        pixel_scale = float(image[0].header['SECPPIX'])
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
                              'flux', 'fluxerr', 'peak', 'background', 'fwhm',
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
            catalog['fwhm'].unit = 'pixel'
            catalog['fwhm'].description = 'FWHM of the object'
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

analysis = input("Enter 'STELLA', 'INT' or 'LCO' for the type of data (default value 'LCO'):\n") or 'LCO'
while analysis!='STELLA' and analysis!='INT' and analysis!='LCO':
   analysis = input("Invalid value, please insert 'STELLA', 'INT' or 'LCO':\n") or 'INT'

#Gets in input the factor for the extraction radius and the detection threshold 
star_scale = float(input("Enter the factor for the extraction radius in the 1st source detection (default value 1.2):\n") or 1.2)
while star_scale<=0.:
            star_scale = float(input("Invalid value, please insert a positive extraction factor:\n"))
s = float(input("Enter the factor s for the detection threshold s*bkg_rms (default value 2.0):\n") or 2.0)
while s<=0.:
            s = float(input("Invalid value, please insert a positive threshold:\n"))

#Gets in input respectively lowcounts, hicounts, lowestcounts, thresholdcounts, closerejectd, targetradius and ignoreedgefraction (check astrosource*.py script for more details)
#lowc = float(input("Enter an integer number for lowcounts to be used in run 1 and 2 of Astrosource:\n"))
#while lowc<0. or lowc%1!=0.:
#            lowc = float(input("Invalid value, please insert a non negative integer number:\n"))
#lowc = int(lowc)

#highc = float(input("Enter an integer number for highcounts to be used in run 1 and 2 of Astrosource:\n"))
#while highc<=0. or highc%1!=0.:
#            highc = float(input("Invalid value, please insert a positive integer number:\n"))
#highc = int(highc)

#lowestc = float(input("Enter an integer number for lowestcounts to be used in run 1 and 2 of Astrosource:\n"))
#while lowestc<0. or lowestc%1!=0.:
#            lowestc = float(input("Invalid value, please insert a non negative integer number:\n"))
#lowestc = int(lowestc)

#threshc = float(input("Enter an integer number for thresholdcounts to be used in run 1 and 2 of Astrosource:\n"))
#while threshc<=0. or threshc%1!=0.:
#            threshc = float(input("Invalid value, please insert a positive integer number:\n"))
#threshc = int(threshc)

#varsthresh = float(input("Enter an integer number for varsearchthresh to be used in run 1 and 2 of Astrosource:\n"))
#while varsthresh<=0. or varsthresh%1!=0.:
#            varsthresh = float(input("Invalid value, please insert a positive integer number:\n"))
#varsthresh = int(varsthresh)

#clrejd = float(input("Enter a positive number for closerejectd (in arcsec) to be used in run 1 and 2 of Astrosource:\n"))
#while clrejd<=0.:
#            clrejd = float(input("Invalid value, please insert a positive number:\n"))

#trad = float(input("Enter a positive number for targetradius (in arcsec) to be used in run 1 and 2 of Astrosource:\n"))
#while trad<=0.:
#            trad = float(input("Invalid value, please insert a positive number:\n"))

#ignedg = float(input("Enter a number between 0 and 1 (extremes included) for ignoreedgefraction to be used in run 1 and 2 of Astrosource:\n"))
#while ignedg<0. or ignedg>1.:
#            ignedg = float(input("Invalid value, please insert a number between 0 and 1 (extremes included):\n"))

start_time = time.time()

image_file = path+'/median_g.fits' 
image = fits.open(path+'/median_g.fits')
photometry(image,image_file=image_file,threshold=2,min_area=2,star_scale=star_scale,min_fwhm=2,max_fwhm=10,analysis=analysis)
print('median_g.fits')
image.close()
print("image should now be closed")

image_file = path+'/median_r.fits' 
image = fits.open(path+'/median_r.fits')
photometry(image,image_file=image_file,threshold=2,min_area=2,star_scale=star_scale,min_fwhm=2,max_fwhm=10,analysis=analysis)
print('median_r.fits')
image.close()
print("image should now be closed")

image_file = path+'/median_i.fits' 
image = fits.open(path+'/median_i.fits')
photometry(image,image_file=image_file,threshold=2,min_area=2,star_scale=star_scale,min_fwhm=2,max_fwhm=10,analysis=analysis)
print('median_i.fits')
image.close()
print("image should now be closed")

images = [f for f in os.listdir(path) if f.endswith('.fits') and not f.startswith('median')]
for i in range(0,(len(images))):
            image_file = path+'/'+images[i] 
            image = fits.open(path+'/'+images[i])
            photometry(image,image_file=image_file,threshold=s,min_area=2,star_scale=star_scale,min_fwhm=2,max_fwhm=10,analysis=analysis)
            print(images[i])
            image.close()
            print("image should now be closed")

end_time = time.time()
execution_time = end_time - start_time
print(f"Execution time: {execution_time:.4f} seconds")

#print(path+'/rfilter')
#os.makedirs('./gfilter', exist_ok=True)
#os.makedirs('./rfilter', exist_ok=True)
#os.makedirs('./ifilter', exist_ok=True)
#if(os.path.isdir(path+'/gfilter')):
#            shutil.rmtree(path+'/gfilter')
#os.system('mkdir '+path+'/gfilter')
#gimages = np.genfromtxt('./gfilter.txt',dtype='str')
#os.system('mv median_g.fits ./gfilter/')
#for i in range(0,(len(gimages))):
#            os.system('mv '+gimages[i]+' ./gfilter/')

#if(os.path.isdir(path+'/rfilter')):
#            shutil.rmtree(path+'/rfilter')
#os.system('mkdir '+path+'/rfilter')
#rimages = np.genfromtxt('./rfilter.txt',dtype='str')
#shutil.copyfile('median_r.fits', './rfilter/median_r.fits')
#for i in range(0,(len(rimages))):
            #shutil.copyfile(rimages[i], ' ./rfilter/'+rimages[i])

#if(os.path.isdir(path+'/ifilter')):
#            shutil.rmtree(path+'/ifilter')
#os.system('mkdir '+path+'/ifilter')
#iimages = np.genfromtxt('./ifilter.txt',dtype='str')
#os.system('mv median_i.fits ifilter/')
#for i in range(0,(len(iimages))):
#            os.system('mv '+iimages[i]+' ./ifilter')

#1st run of Astrosource in the filters g, r and i, to find a set of about 10 stable comparison stars

#Read ra and dec of the target source from one of the fits headers
#image = fits.open(path+'/gfilter/median_g.fits')
#hdr = image[0].header
#ra_hms = hdr['CAT-RA']
#ra_hms = ra_hms.split(':')
#ra = (float(ra_hms[0])+float(ra_hms[1])/60.+float(ra_hms[2])/3600.)*360./24.
#dec_dms = hdr['CAT-DEC']
#dec_dms = dec_dms.split(':')
#if float(dec_dms[0])<0.:
#            dec = float(dec_dms[0])-float(dec_dms[1])/60.-float(dec_dms[2])/3600.
#else:
#            dec = float(dec_dms[0])+float(dec_dms[1])/60.+float(dec_dms[2])/3600.
#image.close()
#ra = 260.232875
#dec = -5.561111

#Perform the 1st Astrosource run to find the comparison stars
#os.chdir(path+'/gfilter')
#os.system('astrosource_LCO_2016 --ra '+str(ra)+' --dec '+str(dec)+' --indir '+path+'/gfilter --format fits --full --lowcounts '+str(lowc)+' --hicounts '+str(highc)+' --lowestcounts '+str(lowestc)+' --thresholdcounts '+str(threshc)+' --closerejectd '+str(clrejd)+' --targetradius '+str(trad)+' --varsearch --varsearchthresh '+str(varsthresh)+' --ignoreedgefraction '+str(ignedg)+' --verbose')
#os.chdir(path+'/rfilter')
#os.system('astrosource_LCO_2016 --ra '+str(ra)+' --dec '+str(dec)+' --indir '+path+'/rfilter --format fits --full --lowcounts '+str(lowc)+' --hicounts '+str(highc)+' --lowestcounts '+str(lowestc)+' --thresholdcounts '+str(threshc)+' --closerejectd '+str(clrejd)+' --targetradius '+str(trad)+' --varsearch --varsearchthresh '+str(varsthresh)+' --ignoreedgefraction '+str(ignedg)+' --verbose')
#os.chdir(path+'/ifilter')
#os.system('astrosource_LCO_2016 --ra '+str(ra)+' --dec '+str(dec)+' --indir '+path+'/ifilter --format fits --full --lowcounts '+str(lowc)+' --hicounts '+str(highc)+' --lowestcounts '+str(lowestc)+' --thresholdcounts '+str(threshc)+' --closerejectd '+str(clrejd)+' --targetradius '+str(trad)+' --varsearch --varsearchthresh '+str(varsthresh)+' --ignoreedgefraction '+str(ignedg)+' --verbose')

#Extract file for coordinates of the comparison star
#gcomps = np.genfromtxt(path+'/gfilter/compsUsed.csv', delimiter=',')
#np.savetxt(path+'/gfilter/compsCoord.csv',[gcomps[:,0],gcomps[:,1]])
#rcomps = np.genfromtxt(path+'/rfilter/compsUsed.csv', delimiter=',')
#np.savetxt(path+'/rfilter/compsCoord.csv',[rcomps[:,0],rcomps[:,1]])
#icomps = np.genfromtxt(path+'/ifilter/compsUsed.csv', delimiter=',')
#np.savetxt(path+'/ifilter/compsCoord.csv',[icomps[:,0],icomps[:,1]])

#Perform the 2nd Astrosource run to extract lightcurves of the comparison stars
#os.chdir(path+'/gfilter')
#os.system('astrosource_LCO_2016 --target-file '+path+'/gfilter/compsCoord.csv --indir '+path+'/gfilter --comparison --phot --plot --targetradius '+str(trad)+' --usescreenedcomps --usecompsused --usecompletedcalib --verbose')
#os.chdir(path+'/rfilter')
#os.system('astrosource_LCO_2016 --target-file '+path+'/rfilter/compsCoord.csv --indir '+path+'/rfilter --comparison --phot --plot --targetradius '+str(trad)+' --usescreenedcomps --usecompsused --usecompletedcalib --verbose')
#os.chdir(path+'/ifilter')
#os.system('astrosource_LCO_2016 --target-file '+path+'/ifilter/compsCoord.csv --indir '+path+'/ifilter --comparison --phot --plot --targetradius '+str(trad)+' --usescreenedcomps --usecompsused --usecompletedcalib --verbose')



#class PhotometricCalibrator(Stage):
#    color_to_fit = {'gp': 'g-r', 'rp': 'r-i', 'ip': 'r-i', 'zs': 'i-z'}
#
#    def __init__(self, runtime_context):
#        super(PhotometricCalibrator, self).__init__(runtime_context)
#
#    def do_stage(self, image):
#        if image.filter not in ['gp', 'rp', 'ip', 'zs']:
#            return image
#
#        if image['CAT'] is None:
#            logger.warning("Not photometrically calibrating image because no catalog exists", image=image)
#            return image
#
#        if image.meta.get('WCSERR', 1) > 0:
#            logger.warning("Not photometrically calibrating image because WCS solution failed", image=image)
#            return image
#
#        try:
#            # Get the sources in the frame
#           reference_catalog = get_reference_sources(image.meta,
#                                                      urljoin(self.runtime_context.REFERENCE_CATALOG_URL, '/image'),
#                                                      nx=image.shape[1], ny=image.shape[0])
#        except HTTPError as e:
#            logger.error(f'Error retrieving photometric reference catalog: {e}', image=image)
#            return image
#
#        # Match the catalog to the detected sources
#        good_sources = np.logical_and(image['CAT'].data['flag'] == 0, image['CAT'].data['flux'] > 0.0)
#        matched_catalog = match_catalogs(image['CAT'].data[good_sources], reference_catalog)
#
#        if len(matched_catalog) == 0:
#            logger.error('No matching sources found. Skipping zeropoint determination', image=image)
#            return image
#        # catalog_mag = instrumental_mag + zeropoint + color_coefficient * color
#        # Fit the zeropoint and color_coefficient rejecting outliers
#        # Note the zero index here in the filter name is because we only store teh first letter of the filter name
#        try:
#            zeropoint, zeropoint_error, color_coefficient, color_coefficient_error = \
#                fit_photometry(matched_catalog, image.filter[0], self.color_to_fit[image.filter], image.exptime)
#        except:
#            logger.error(f"Error fitting photometry: {logs.format_exception()}", image=image)
#            return image
#
#        # Save the zeropoint, color coefficient and errors to header
#        image.meta['L1ZP'] = zeropoint, "Instrumental zeropoint [mag]"
#        image.meta['L1ZPERR'] = zeropoint_error, "Error on Instrumental zeropoint [mag]"
#        image.meta['L1COLORU'] = self.color_to_fit[image.filter], "Color used for calibration"
#        image.meta['L1COLOR'] = color_coefficient, "Color coefficient [mag]"
#        image.meta['L1COLERR'] = color_coefficient_error, "Error on color coefficient [mag]"
#        # Calculate the mag of each of the items in the catalog (without the color term) saving them
#        image['CAT'].data['mag'], image['CAT'].data['magerr'] = to_magnitude(image['CAT'].data['flux'], image['CAT'].data['fluxerr'],
#                                                                             zeropoint, image.exptime)
#        return image
