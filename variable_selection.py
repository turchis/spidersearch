#!/bin/python

import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits
from astropy import wcs
from astropy.coordinates import SkyCoord
#import statistics as st
import os

path=os.getcwd()
#Lower limit of for getting rid of saturate/too bright sources to be spiders
lowerdm=-2.
#Lower limit of stdev for selecting candidates
#lower_stdev=0.1
#Number of decimal digits to which round the coordinates for identifying candidates showing variability simultaneously in all three bands
ndig=3
#Width of the step in dm
dm_step = 0.1
#Minimum number of stars detected in the magnitude bin to perform a statistical variability analysis
N_min = 5
#First cut-off of sigma_target/sigma_median for selecting candidates
cut_1 = 1.
#Second cut-off of sigma_target/sigma_median for selecting candidates
cut_2 = 0.0
#Third cut-off of sigma_target/sigma_median for selecting candidates
#cut_3 = 1.0

#Defining scale factor to use as squared region around the target, with side=2*bound*semimaj_axis
bound = 1.0

analysis = input("Enter 'STELLA', 'INT' or 'LCO' for the type of data (default value 'LCO'):\n") or 'LCO'
while analysis!='STELLA' and analysis!='INT' and analysis!='LCO':
   analysis = input("Invalid value, please insert 'STELLA', 'INT' or 'LCO':\n") or 'LCO'

#gfilter

#Read the target name from one of the fits headers and load the corresponding Fermi coordinates and ellipses (in degrees)
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
w = wcs.WCS(hdr)
image.close()
#Load the corresponding Fermi coordinates and ellipses (in degrees) of the target

#print(ra,dec)

hdul = fits.open('/export/work/marcotu/gll_psc_v28.fit')
[[ra, dec, ell_a, ell_b]] = [[hdul[1].data['RAJ2000'][i],hdul[1].data['DEJ2000'][i],3600./pixscale*hdul[1].data['Conf_95_SemiMajor'][i],3600./pixscale*hdul[1].data['Conf_95_SemiMinor'][i]] for i in range(0,len(hdul[1].data)) if FoVname.replace("L","L ")==hdul[1].data['Source_Name'][i] for i in range(0,len(hdul[1].data)) if FoVname.replace("L","L ")==hdul[1].data['Source_Name'][i]]
hdul.close()
                               
sources_g = np.genfromtxt(path+'/gfilter/starVariability.csv', delimiter=',')
sources_g = np.asarray(sources_g)
dmbins_g = np.arange(np.min(sources_g[:,2]), np.max(sources_g[:,2]), dm_step)

xsources, ysources = w.all_world2pix(sources_g[:,0], sources_g[:,1], 1)
xtar, ytar = w.all_world2pix(ra, dec, 1)
 
#Producing csv files with parameters of the variable sources in the gband with sigmam >= cut_1*median (list1), sigma >= cut_2*median sigma<cut_1*median (list2) and sigma >= cut_3*median sigma<cut_2*median (list3), with also the corresponding region files 
#if os.path.exists(path+'/gfilter/varSources_1.csv'):
#    os.remove(path+'/gfilter/varSources_1.csv')
var1 = open(path+'/gfilter/varSources_1_COBIPULSEmeth.csv','w')
#if os.path.exists(path+'/gfilter/varSources_2.csv'):
#    os.remove(path+'/gfilter/varSources_2.csv')
#var2 = open(path+'/gfilter/varSources_2.csv','w')
#if os.path.exists(path+'/gfilter/varSourcescoord_1.csv'):
#    os.remove(path+'/gfilter/varSourcescoord_1.csv')
varc1 = open(path+'/gfilter/varSourcescoord_1_COBIPULSEmeth.csv','w')
#if os.path.exists(path+'/gfilter/varSourcescoord_2.csv'):
#    os.remove(path+'/gfilter/varSourcescoord_2.csv')
#varc2 = open(path+'/gfilter/varSourcescoord_2.csv','w')
#if os.path.exists(path+'/Sourcesg_1.reg'):
#    os.remove(path+'/Sourcesg_1.reg')
reg1 = open(path+'/Sourcesg_1_COBIPULSEmeth.reg','w')
#if os.path.exists(path+'/Sourcesg_2.reg'):
#    os.remove(path+'/Sourcesg_2.reg')
#reg2 = open(path+'/Sourcesg_2.reg','w')

print("# Region file format: DS9 version 4.1 (via variable_selection.py)",file=reg1)
print("global dashlist=8 3 width=1 font=\"helvetica 10 normal roman\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1",file=reg1)
print("fk5",file=reg1)
#print("# Region file format: DS9 version 4.1 (via variable_selection.py)",file=reg2)
#print("global dashlist=8 3 width=1 font=\"helvetica 10 normal roman\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1",file=reg2)
#print("fk5",file=reg2)

#Array for the variable sources in the gband with sigmam >= cut_1*median
varSourcesg_1 = []
#Array for the variable sources in the gband with cut_2*median <= sigmam < cut_1*median 
#varSourcesg_2 = []
counter1 = 0
#counter2 = 0
for i in range(0,len(dmbins_g)):
    tempsources = []
    for j in range(0,len(sources_g)):
        if (sources_g[j,2] >= dmbins_g[i] and sources_g[j,2] < dmbins_g[i]+dm_step):
            tempsources.append(sources_g[j,3])
    #Performing a variability analysis only when the statistics is not too poor (number of stars in the mag bin >= N_min)
    if len(tempsources) >= N_min:
                median = np.median(tempsources)
                for j in range(0,len(sources_g)):
                            #Retaining only the sources inside the squared region with side=2*bound*semimaj_axis and with dm>=lowerdm
                            if (sources_g[j,2] >= lowerdm and sources_g[j,2] >= dmbins_g[i] and sources_g[j,2] < dmbins_g[i]+dm_step and xsources[j]>=xtar-bound*max(ell_a,ell_b) and xsources[j]<=xtar+bound*max(ell_a,ell_b) and ysources[j]>=ytar-bound*max(ell_a,ell_b) and ysources[j]<=ytar+bound*max(ell_a,ell_b)):
                                        if (sources_g[j,3]>=cut_1*median):
                                                    counter1+=1
                                                    varSourcesg_1.append([sources_g[j,0],sources_g[j,1],sources_g[j,2],sources_g[j,3]])
                                                    print(sources_g[j,0], sources_g[j,1], sources_g[j,2], sources_g[j,3], sep=',', file=var1)
                                                    print(sources_g[j,0], sources_g[j,1], sep=',', file=varc1)
                                                    print("circle(%f,%f,3.0\") # text={%d} color=red width=1" % (sources_g[j,0],sources_g[j,1],int(counter1)),file=reg1)
                                        #elif (sources_g[j,3]>=cut_2*median and sources_g[j,3]<cut_1*median):
                                        #            counter2+=1
                                        #            varSourcesg_2.append([sources_g[j,0],sources_g[j,1],sources_g[j,2],sources_g[j,3]])
                                        #            print(sources_g[j,0], sources_g[j,1], sources_g[j,2], sources_g[j,3], sep=',', file=var2)
                                        #            print(sources_g[j,0], sources_g[j,1], sep=',', file=varc2)
                                        #            print("circle(%f,%f,3.0\") # text={%d} color=green width=1" % (sources_g[j,0],sources_g[j,1],int(counter2)),file=reg2)
    else:
                for j in range(0,len(sources_g)):
                            #Retaining only the sources inside the squared region with side=2*bound*semimaj_axis and with dm>=lowerdm
                            if (sources_g[j,2] >= lowerdm and sources_g[j,2] >= dmbins_g[i] and sources_g[j,2] < dmbins_g[i]+dm_step and xsources[j]>=xtar-bound*max(ell_a,ell_b) and xsources[j]<=xtar+bound*max(ell_a,ell_b) and ysources[j]>=ytar-bound*max(ell_a,ell_b) and ysources[j]<=ytar+bound*max(ell_a,ell_b)):
                                        counter1+=1
                                        varSourcesg_1.append([sources_g[j,0],sources_g[j,1],sources_g[j,2],sources_g[j,3]])
                                        print(sources_g[j,0], sources_g[j,1], sources_g[j,2], sources_g[j,3], sep=',', file=var1)
                                        print(sources_g[j,0], sources_g[j,1], sep=',', file=varc1)
                                        print("circle(%f,%f,3.0\") # text={%d} color=red width=1" % (sources_g[j,0],sources_g[j,1],int(counter1)),file=reg1)

var1.close()
#var2.close()
varc1.close()
#varc2.close()
reg1.close()
#reg2.close()

varSourcesg_1 = np.asarray(varSourcesg_1)
#varSourcesg_2 = np.asarray(varSourcesg_2)

print('Number of candidate variable sources found in the g band above cut_1:',len(varSourcesg_1))

#rfilter

image = fits.open(path+'/rfilter/median_r.fits')
hdr = image[0].header
w = wcs.WCS(hdr)
image.close()

sources_r = np.genfromtxt(path+'/rfilter/starVariability.csv', delimiter=',')
sources_r = np.asarray(sources_r)
dmbins_r = np.arange(np.min(sources_r[:,2]), np.max(sources_r[:,2]), dm_step)

xsources, ysources = w.all_world2pix(sources_r[:,0], sources_r[:,1], 1)
xtar, ytar = w.all_world2pix(ra, dec, 1)

#Producing csv files with parameters of the variable sources in the rband with sigmam >= cut_1*median (list1) and sigma >= cut_2*median sigma<cut_1*median (list2), with also the region files
#if os.path.exists(path+'/rfilter/varSources_1.csv'):
#    os.remove(path+'/rfilter/varSources_1.csv')
var1 = open(path+'/rfilter/varSources_1_COBIPULSEmeth.csv','w')
#if os.path.exists(path+'/rfilter/varSources_2.csv'):
#    os.remove(path+'/rfilter/varSources_2.csv')
#var2 = open(path+'/rfilter/varSources_2.csv','w')
#if os.path.exists(path+'/rfilter/varSourcescoord_1.csv'):
#    os.remove(path+'/rfilter/varSourcescoord_1.csv')
varc1 = open(path+'/rfilter/varSourcescoord_1_COBIPULSEmeth.csv','w')
#if os.path.exists(path+'/rfilter/varSourcescoord_2.csv'):
#    os.remove(path+'/rfilter/varSourcescoord_2.csv')
#varc2 = open(path+'/rfilter/varSourcescoord_2.csv','w')
#if os.path.exists(path+'/Sourcesr_1.reg'):
#    os.remove(path+'/Sourcesr_1.reg')
reg1 = open(path+'/Sourcesr_1_COBIPULSEmeth.reg','w')
#if os.path.exists(path+'/Sourcesr_2.reg'):
#    os.remove(path+'/Sourcesr_2.reg')
#reg2 = open(path+'/Sourcesr_2.reg','w')

print("# Region file format: DS9 version 4.1 (via variable_selection.py)",file=reg1)
print("global dashlist=8 3 width=1 font=\"helvetica 10 normal roman\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1",file=reg1)
print("fk5",file=reg1)
#print("# Region file format: DS9 version 4.1 (via variable_selection.py)",file=reg2)
#print("global dashlist=8 3 width=1 font=\"helvetica 10 normal roman\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1",file=reg2)
#print("fk5",file=reg2)

varSourcesr_1 = []
#varSourcesr_2 = []
counter1 = 0
#counter2 = 0
for i in range(0,len(dmbins_r)):
    tempsources = []
    for j in range(0,len(sources_r)):
        if (sources_r[j,2] >= dmbins_r[i] and sources_r[j,2] < dmbins_r[i]+dm_step):
            tempsources.append(sources_r[j,3])
    #Performing a variability analysis only when the statistics is not too poor (number of stars in the mag bin >= N_min)
    if len(tempsources) >= N_min:
                median = np.median(tempsources)
                for j in range(0,len(sources_r)):
                            #Retaining only the sources inside the squared region with side=2*bound*semimaj_axis and with dm>=lowerdm
                            if (sources_r[j,2] >= lowerdm and sources_r[j,2] >= dmbins_r[i] and sources_r[j,2] < dmbins_r[i]+dm_step and xsources[j]>=xtar-bound*max(ell_a,ell_b) and xsources[j]<=xtar+bound*max(ell_a,ell_b) and ysources[j]>=ytar-bound*max(ell_a,ell_b) and ysources[j]<=ytar+bound*max(ell_a,ell_b)):
                                        if (sources_r[j,3]>=cut_1*median):
                                                    counter1+=1
                                                    varSourcesr_1.append([sources_r[j,0],sources_r[j,1],sources_r[j,2],sources_r[j,3]])
                                                    print(sources_r[j,0], sources_r[j,1], sources_r[j,2], sources_r[j,3], sep=',', file=var1)
                                                    print(sources_r[j,0], sources_r[j,1], sep=',', file=varc1)
                                                    print("circle(%f,%f,3.0\") # text={%d} color=red width=1" % (sources_r[j,0],sources_r[j,1],int(counter1)),file=reg1)
                                        #elif (sources_r[j,3]>=cut_2*median and sources_r[j,3]<cut_1*median):
                                        #            counter2+=1
                                        #            varSourcesr_2.append([sources_r[j,0],sources_r[j,1],sources_r[j,2],sources_r[j,3]])
                                        #            print(sources_r[j,0], sources_r[j,1], sources_r[j,2], sources_r[j,3], sep=',', file=var2)
                                        #            print(sources_r[j,0], sources_r[j,1], sep=',', file=varc2)
                                        #            print("circle(%f,%f,3.0\") # text={%d} color=green width=1" % (sources_r[j,0],sources_r[j,1],int(counter2)),file=reg2)
    else:
                for j in range(0,len(sources_r)):
                            #Retaining only the sources inside the squared region with side=2*bound*semimaj_axis and with dm>=lowerdm
                            if (sources_r[j,2] >= lowerdm and sources_r[j,2] >= dmbins_r[i] and sources_r[j,2] < dmbins_r[i]+dm_step and xsources[j]>=xtar-bound*max(ell_a,ell_b) and xsources[j]<=xtar+bound*max(ell_a,ell_b) and ysources[j]>=ytar-bound*max(ell_a,ell_b) and ysources[j]<=ytar+bound*max(ell_a,ell_b)):
                                        counter1+=1
                                        varSourcesr_1.append([sources_r[j,0],sources_r[j,1],sources_r[j,2],sources_r[j,3]])
                                        print(sources_r[j,0], sources_r[j,1], sources_r[j,2], sources_r[j,3], sep=',', file=var1)
                                        print(sources_r[j,0], sources_r[j,1], sep=',', file=varc1)
                                        print("circle(%f,%f,3.0\") # text={%d} color=red width=1" % (sources_r[j,0],sources_r[j,1],int(counter1)),file=reg1)

var1.close()
#var2.close()
varc1.close()
#varc2.close()
reg1.close()
#reg2.close()

varSourcesr_1 = np.asarray(varSourcesr_1)
#varSourcesr_2 = np.asarray(varSourcesr_2)

print('Number of candidate variable sources found in the r band above cut_1:',len(varSourcesr_1))

#ifilter

image = fits.open(path+'/ifilter/median_i.fits')
hdr = image[0].header
w = wcs.WCS(hdr)
image.close()
            
sources_i = np.genfromtxt(path+'/ifilter/starVariability.csv', delimiter=',')
sources_i = np.asarray(sources_i)
dmbins_i = np.arange(np.min(sources_i[:,2]), np.max(sources_i[:,2]), dm_step)

xsources, ysources = w.all_world2pix(sources_i[:,0], sources_i[:,1], 1)
xtar, ytar = w.all_world2pix(ra, dec, 1)

#Producing csv files with parameters of the variable sources in the iband with sigmam >= cut_1*median (list1) and sigma >= cut_2*median sigma<cut_1*median (list2), with also the region files
#if os.path.exists(path+'/ifilter/varSources_1.csv'):
#    os.remove(path+'/ifilter/varSources_1.csv')
var1 = open(path+'/ifilter/varSources_1_COBIPULSEmeth.csv','w')
#if os.path.exists(path+'/ifilter/varSources_2.csv'):
#    os.remove(path+'/ifilter/varSources_2.csv')
#var2 = open(path+'/ifilter/varSources_2.csv','w')
#if os.path.exists(path+'/ifilter/varSourcescoord_1.csv'):
#    os.remove(path+'/ifilter/varSourcescoord_1.csv')
varc1 = open(path+'/ifilter/varSourcescoord_1_COBIPULSEmeth.csv','w')
#if os.path.exists(path+'/ifilter/varSourcescoord_2.csv'):
#    os.remove(path+'/ifilter/varSourcescoord_2.csv')
#varc2 = open(path+'/ifilter/varSourcescoord_2.csv','w')
#if os.path.exists(path+'/Sourcesi_1.reg'):
#    os.remove(path+'/Sourcesi_1.reg')
reg1 = open(path+'/Sourcesi_1_COBIPULSEmeth.reg','w')
#if os.path.exists(path+'/Sourcesi_2.reg'):
#    os.remove(path+'/Sourcesi_2.reg')
#reg2 = open(path+'/Sourcesi_2.reg','w')

print("# Region file format: DS9 version 4.1 (via variable_selection.py)",file=reg1)
print("global dashlist=8 3 width=1 font=\"helvetica 10 normal roman\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1",file=reg1)
print("fk5",file=reg1)
#print("# Region file format: DS9 version 4.1 (via variable_selection.py)",file=reg2)
#print("global dashlist=8 3 width=1 font=\"helvetica 10 normal roman\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1",file=reg2)
#print("fk5",file=reg2)

varSourcesi_1 = []
#varSourcesi_2 = []
counter1 = 0
#counter2 = 0
for i in range(0,len(dmbins_i)):
    tempsources = []
    for j in range(0,len(sources_i)):
        if (sources_i[j,2] >= dmbins_i[i] and sources_i[j,2] < dmbins_i[i]+dm_step):
            tempsources.append(sources_i[j,3])
    #Performing a variability analysis only when the statistics is not too poor (number of stars in the mag bin >= N_min)
    if len(tempsources) >= N_min:
                median = np.median(tempsources)
                for j in range(0,len(sources_i)):
                            #Retaining only the sources inside the squared region with side=2*bound*semimaj_axis and with dm>=lowerdm
                            if (sources_i[j,2] >= lowerdm and sources_i[j,2] >= dmbins_i[i] and sources_i[j,2] < dmbins_i[i]+dm_step and xsources[j]>=xtar-bound*max(ell_a,ell_b) and xsources[j]<=xtar+bound*max(ell_a,ell_b) and ysources[j]>=ytar-bound*max(ell_a,ell_b) and ysources[j]<=ytar+bound*max(ell_a,ell_b)):
                                        if (sources_i[j,3]>=cut_1*median):
                                                    counter1+=1
                                                    varSourcesi_1.append([sources_i[j,0],sources_i[j,1],sources_i[j,2],sources_i[j,3]])
                                                    print(sources_i[j,0], sources_i[j,1], sources_i[j,2], sources_i[j,3], sep=',', file=var1)
                                                    print(sources_i[j,0], sources_i[j,1], sep=',', file=varc1)
                                                    print("circle(%f,%f,3.0\") # text={%d} color=red width=1" % (sources_i[j,0],sources_i[j,1],int(counter1)),file=reg1)
                                        #elif (sources_i[j,3]>=cut_2*median and sources_i[j,3]<cut_1*median):
                                        #            counter2+=1
                                        #            varSourcesi_2.append([sources_i[j,0],sources_i[j,1],sources_i[j,2],sources_i[j,3]])
                                        #            print(sources_i[j,0], sources_i[j,1], sources_i[j,2], sources_i[j,3], sep=',', file=var2)
                                        #            print(sources_i[j,0], sources_i[j,1], sep=',', file=varc2)
                                        #            print("circle(%f,%f,3.0\") # text={%d} color=green width=1" % (sources_i[j,0],sources_i[j,1],int(counter2)),file=reg2)
    else:
                for j in range(0,len(sources_i)):
                            #Retaining only the sources inside the squared region with side=2*bound*semimaj_axis and with dm>=lowerdm
                            if (sources_i[j,2] >= lowerdm and sources_i[j,2] >= dmbins_i[i] and sources_i[j,2] < dmbins_i[i]+dm_step and xsources[j]>=xtar-bound*max(ell_a,ell_b) and xsources[j]<=xtar+bound*max(ell_a,ell_b) and ysources[j]>=ytar-bound*max(ell_a,ell_b) and ysources[j]<=ytar+bound*max(ell_a,ell_b)):
                                        counter1+=1
                                        varSourcesi_1.append([sources_i[j,0],sources_i[j,1],sources_i[j,2],sources_i[j,3]])
                                        print(sources_i[j,0], sources_i[j,1], sources_i[j,2], sources_i[j,3], sep=',', file=var1)
                                        print(sources_i[j,0], sources_i[j,1], sep=',', file=varc1)
                                        print("circle(%f,%f,3.0\") # text={%d} color=red width=1" % (sources_i[j,0],sources_i[j,1],int(counter1)),file=reg1)

var1.close()
#var2.close()
varc1.close()
#varc2.close()
reg1.close()
#reg2.close()

varSourcesi_1 = np.asarray(varSourcesi_1)
#varSourcesi_2 = np.asarray(varSourcesi_2)

print('Number of candidate variable sources found in the i band above cut_1:',len(varSourcesi_1))

#Plotting the three sigma plots for band g, r and i with variable sources in separated bands (grey all sources, red variable candidates list1, green variable candidates list2)
markersize=8
figure=plt.figure(figsize=(12,12), dpi=150)
plt.subplot(3, 1, 1)
plt.xlabel(r'$\Delta$m (mag in the g band)')
plt.ylabel(r'$\sigma$ (mag)')
plt.yscale("log")
plt.plot(sources_g[:,2],sources_g[:,3],'ob',markerfacecolor='0.8',markeredgecolor='black',markersize=markersize)
plt.plot(varSourcesg_1[:,2],varSourcesg_1[:,3],'ob',markerfacecolor='red',markeredgecolor='black',markersize=markersize)
#plt.plot(varSourcesg_2[:,2],varSourcesg_2[:,3],'ob',markerfacecolor='green',markeredgecolor='black',markersize=markersize)

plt.subplot(3, 1, 2)
plt.xlabel(r'$\Delta$m (mag in the r band)')
plt.ylabel(r'$\sigma$ (mag)')
plt.yscale("log")
plt.plot(sources_r[:,2],sources_r[:,3],'ob',markerfacecolor='0.8',markeredgecolor='black',markersize=markersize)
plt.plot(varSourcesr_1[:,2],varSourcesr_1[:,3],'ob',markerfacecolor='red',markeredgecolor='black',markersize=markersize)
#plt.plot(varSourcesr_2[:,2],varSourcesr_2[:,3],'ob',markerfacecolor='green',markeredgecolor='black',markersize=markersize)
                                
plt.subplot(3, 1, 3)
plt.xlabel(r'$\Delta$m (mag in the i band)') 
plt.ylabel(r'$\sigma$ (mag)')
plt.yscale("log")
plt.plot(sources_i[:,2],sources_i[:,3],'ob',markerfacecolor='0.8',markeredgecolor='black',markersize=markersize)
plt.plot(varSourcesi_1[:,2],varSourcesi_1[:,3],'ob',markerfacecolor='red',markeredgecolor='black',markersize=markersize)
#plt.plot(varSourcesi_2[:,2],varSourcesi_2[:,3],'ob',markerfacecolor='green',markeredgecolor='black',markersize=markersize)

#plt.show()
#figure = plt.gcf()
#figure.set_size_inches(15, 15)
plt.savefig(path+"/sigmavsdelta_gri_"+str(2*bound)+"sigma_COBIPULSEmeth.png")
