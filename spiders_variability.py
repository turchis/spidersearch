#!/home/gudrun/marcotu/python-envs/astrosource/bin/python

import os
import numpy as np
from astropy.io import fits

path = os.getcwd()

analysis = input("Enter 'STELLA', 'INT' or 'LCO' for the type of data (default value 'LCO'):\n") or 'LCO'
while analysis!='STELLA' and analysis!='INT' and analysis!='LCO':
   analysis = input("Invalid value, please insert 'STELLA', 'INT' or 'LCO':\n") or 'LCO'

#Gets in input respectively lowcounts, hicounts, lowestcounts, thresholdcounts, closerejectd, targetradius and ignoreedgefraction (check astrosource*.py script for more details)
lowc = float(input("Enter an integer number for lowcounts to be used in run 1 and 2 of Astrosource (default value 200):\n") or "200")
while lowc<0. or lowc%1!=0.:
    lowc = float(input("Invalid value, please insert a non negative integer number (default value 200):\n") or "200")
lowc = int(lowc)

highc_g = float(input("Enter an integer number for highcounts in g-band to be used in run 1 and 2 of Astrosource (default value 50000000):\n") or "50000000")
while highc_g<=0. or highc_g%1!=0.:
    highc_g = float(input("Invalid value, please insert a positive integer number:\n") or "50000000")
highc_g = int(highc_g)

highc_r = float(input("Enter an integer number for highcounts in r-band to be used in run 1 and 2 of Astrosource (default value 50000000):\n") or "50000000")
while highc_r<=0. or highc_r%1!=0.:
    highc_r = float(input("Invalid value, please insert a positive integer number:\n") or "50000000")
highc_r = int(highc_r)

highc_i = float(input("Enter an integer number for highcounts in i-band to be used in run 1 and 2 of Astrosource (default value 50000000):\n") or "50000000")
while highc_i<=0. or highc_i%1!=0.:
    highc_i = float(input("Invalid value, please insert a positive integer number:\n") or "50000000")
highc_i = int(highc_i)

lowestc = float(input("Enter an integer number for lowestcounts to be used in run 1 and 2 of Astrosource (default value 0):\n") or "0")
while lowestc<0. or lowestc%1!=0.:
    lowestc = float(input("Invalid value, please insert a non negative integer number:\n") or "0")
lowestc = int(lowestc)

threshc_g = float(input("Enter an integer number for thresholdcounts in g-band to be used in run 1 and 2 of Astrosource (default value 2000000):\n") or "2000000")
while threshc_g<=0. or threshc_g%1!=0.:
    threshc_g = float(input("Invalid value, please insert a positive integer number:\n") or "2000000")
threshc_g = int(threshc_g)

threshc_r = float(input("Enter an integer number for thresholdcounts in r-band to be used in run 1 and 2 of Astrosource (default value 1500000):\n") or "1500000")
while threshc_r<=0. or threshc_r%1!=0.:
    threshc_r = float(input("Invalid value, please insert a positive integer number:\n") or "1500000")
threshc_r = int(threshc_r)

threshc_i = float(input("Enter an integer number for thresholdcounts in i-band to be used in run 1 and 2 of Astrosource (default value 1500000):\n") or "1500000")
while threshc_i<=0. or threshc_i%1!=0.:
    threshc_i = float(input("Invalid value, please insert a positive integer number:\n") or "1500000")
threshc_i = int(threshc_i)

varsthresh = float(input("Enter an integer number for varsearchthresh to be used in run 1 and 2 of Astrosource (default value 0):\n") or "0")
while varsthresh<0. or varsthresh%1!=0.:
    varsthresh = float(input("Invalid value, please insert a positive integer number:\n") or "0")
varsthresh = int(varsthresh)

clrejd = float(input("Enter a positive number for closerejectd (in arcsec) to be used in run 1 and 2 of Astrosource (default value 3):\n") or "3")
while clrejd<=0.:
    clrejd = float(input("Invalid value, please insert a positive number:\n") or "3")

trad = float(input("Enter a positive number for targetradius (in arcsec) to be used in run 1 and 2 of Astrosource (default value 5):\n") or "5")
while trad<=0.:
    trad = float(input("Invalid value, please insert a positive number:\n") or "5")

ignedg = float(input("Enter a number between 0 and 1 (extremes included) for ignoreedgefraction to be used in run 1 and 2 of Astrosource (default value 0):\n") or "0")
while ignedg<0. or ignedg>1.:
    ignedg = float(input("Invalid value, please insert a number between 0 and 1 (extremes included):\n") or "0")

#Read the target name from one of the fits headers and load the corresponding Fermi coordinates and ellipses (in degrees)
if(os.path.isdir(path+'/gfilter')):
    image = fits.open(path+'/gfilter/median_g.fits')
else:
    image = fits.open(path+'/median_g.fits')
hdr = image[0].header

if analysis == 'STELLA':
   FoVname = hdr['OBJNAME']
   ra = hdr['OBJRA']
   dec = hdr['OBJDEC']
elif analysis == 'INT':
   FoVname = hdr['OBJECT']
   ra = (float(hdr['RA'].split(':')[0])+float(hdr['RA'].split(':')[1])/60.+float(hdr['RA'].split(':')[2])/3600.)*360./24.
   if float(hdr['DEC'].split(':')[0])>0.:
      dec = (float(hdr['DEC'].split(':')[0])+float(hdr['DEC'].split(':')[1])/60.+float(hdr['DEC'].split(':')[2])/3600.)
   else:
      dec = (float(hdr['DEC'].split(':')[0])-float(hdr['DEC'].split(':')[1])/60.-float(hdr['DEC'].split(':')[2])/3600.)
else:
   FoVname = hdr['OBJECT']
   ra = (float(hdr['CAT-RA'].split(':')[0])+float(hdr['CAT-RA'].split(':')[1])/60.+float(hdr['CAT-RA'].split(':')[2])/3600.)*360./24.
   if float(hdr['CAT-DEC'].split(':')[0])>0.:
      dec = (float(hdr['CAT-DEC'].split(':')[0])+float(hdr['CAT-DEC'].split(':')[1])/60.+float(hdr['CAT-DEC'].split(':')[2])/3600.)
   else:
      dec = (float(hdr['CAT-DEC'].split(':')[0])-float(hdr['CAT-DEC'].split(':')[1])/60.-float(hdr['CAT-DEC'].split(':')[2])/3600.)
            
image.close()

#Load the corresponding Fermi coordinates and ellipses (in degrees) of the target
#with open('/export/work/marcotu/gll_psc_v16.reg', 'r') as infile:
#            for line in infile:
#                        if FoVname.replace("L","L ") in line:
#                                    ra = float(line.replace(" ","").split("(")[1].split(",")[0])
#                                    dec = float(line.replace(" ","").split("(")[1].split(",")[1].split(")")[0])
#print(ra,dec)

OLD_PATH=os.environ["PATH"]
ENV_PATH="/home/gudrun/marcotu/python-envs/astrosource/bin"
os.environ["PATH"]=ENV_PATH + ":" + OLD_PATH

#os.system('source /home/gudrun/marcotu/python-envs/astrosource/bin/activate')
os.system('osfits_gri.sh')
#os.system('osfits_griz.sh')
#Perform the 1st Astrosource run to find the comparison stars
os.chdir(path+'/gfilter')
os.system('astrosource_'+analysis+' --ra '+str(ra)+' --dec '+str(dec)+' --indir $(pwd) --format fits --full --lowcounts '+str(lowc)+' --hicounts '+str(highc_g)+' --lowestcounts '+str(lowestc)+' --thresholdcounts '+str(threshc_g)+' --closerejectd '+str(clrejd)+' --targetradius '+str(trad)+' --no2mass --varsearch --varsearchthresh '+str(varsthresh)+' --ignoreedgefraction '+str(ignedg)+' --verbose --bjd')
os.chdir(path+'/rfilter')
os.system('astrosource_'+analysis+' --ra '+str(ra)+' --dec '+str(dec)+' --indir $(pwd) --format fits --full --lowcounts '+str(lowc)+' --hicounts '+str(highc_r)+' --lowestcounts '+str(lowestc)+' --thresholdcounts '+str(threshc_r)+' --closerejectd '+str(clrejd)+' --targetradius '+str(trad)+' --no2mass --varsearch --varsearchthresh '+str(varsthresh)+' --ignoreedgefraction '+str(ignedg)+' --verbose --bjd')
os.chdir(path+'/ifilter')
os.system('astrosource_'+analysis+' --ra '+str(ra)+' --dec '+str(dec)+' --indir $(pwd) --format fits --full --lowcounts '+str(lowc)+' --hicounts '+str(highc_i)+' --lowestcounts '+str(lowestc)+' --thresholdcounts '+str(threshc_i)+' --closerejectd '+str(clrejd)+' --targetradius '+str(trad)+' --no2mass --varsearch --varsearchthresh '+str(varsthresh)+' --ignoreedgefraction '+str(ignedg)+' --verbose --bjd')
#os.chdir(path+'/zfilter')
#os.system('astrosource_'+analysis+' --ra '+str(ra)+' --dec '+str(dec)+' --indir $(pwd) --format fits --full --lowcounts '+str(lowc)+' --hicounts '+str(highc)+' --lowestcounts '+str(lowestc)+' --thresholdcounts '+str(threshc)+' --closerejectd '+str(clrejd)+' --targetradius '+str(trad)+' --no2mass --varsearch --varsearchthresh '+str(varsthresh)+' --varsearchminimages 0. --ignoreedgefraction '+str(ignedg)+' --verbose --bjd')
#os.system('astrosource_'+analysis+' --ra '+str(ra)+' --dec '+str(dec)+' --indir $(pwd) --format fits --full --lowcounts '+str(lowc)+' --hicounts '+str(highc)+' --lowestcounts '+str(lowestc)+' --thresholdcounts 200000 --closerejectd '+str(clrejd)+' --targetradius '+str(trad)+' --no2mass --varsearch --varsearchthresh '+str(varsthresh)+' --ignoreedgefraction '+str(ignedg)+' --verbose --bjd')

#Getting the coords files for comparison stars
ra_c, dec_c = np.genfromtxt(path+'/gfilter/compsUsed.csv', delimiter=',',usecols=(0,1),unpack=True)
np.savetxt(path+"/gfilter/compsCoord.csv", np.column_stack((ra_c, dec_c)), delimiter=",", fmt='%0.8f')
reg = open(path+'/gfilter/compsCoord.reg','w')
print("# Region file format: DS9 version 4.1",file=reg)
print("global dashlist=8 3 width=1 font=\"helvetica 10 normal roman\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1",file=reg)
print("fk5",file=reg)
for i in range(0,np.size(ra_c)):
   if np.size(ra_c) > 1:
      print("circle(%f,%f,3.0\") # color=red width=1" % (ra_c[i],dec_c[i]),file=reg)
   else:
      print("circle(%f,%f,3.0\") # color=red width=1" % (ra_c,dec_c),file=reg)
reg.close()
ra_c, dec_c = np.genfromtxt(path+'/rfilter/compsUsed.csv', delimiter=',',usecols=(0,1),unpack=True)
np.savetxt(path+"/rfilter/compsCoord.csv", np.column_stack((ra_c, dec_c)), delimiter=",", fmt='%0.8f')
reg = open(path+'/rfilter/compsCoord.reg','w')
print("# Region file format: DS9 version 4.1",file=reg)
print("global dashlist=8 3 width=1 font=\"helvetica 10 normal roman\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1",file=reg)
print("fk5",file=reg)
for i in range(0,np.size(ra_c)):
   if np.size(ra_c) > 1:
      print("circle(%f,%f,3.0\") # color=red width=1" % (ra_c[i],dec_c[i]),file=reg)
   else:
      print("circle(%f,%f,3.0\") # color=red width=1" % (ra_c,dec_c),file=reg)
reg.close()
ra_c, dec_c = np.genfromtxt(path+'/ifilter/compsUsed.csv', delimiter=',',usecols=(0,1),unpack=True)
reg = open(path+'/ifilter/compsCoord.reg','w')
print("# Region file format: DS9 version 4.1",file=reg)
print("global dashlist=8 3 width=1 font=\"helvetica 10 normal roman\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1",file=reg)
print("fk5",file=reg)
for i in range(0,np.size(ra_c)):
   if np.size(ra_c) > 1:
      print("circle(%f,%f,3.0\") # color=red width=1" % (ra_c[i],dec_c[i]),file=reg)
   else:
      print("circle(%f,%f,3.0\") # color=red width=1" % (ra_c,dec_c),file=reg)
reg.close()
np.savetxt(path+"/ifilter/compsCoord.csv", np.column_stack((ra_c, dec_c)), delimiter=",", fmt='%0.8f')
#ra_c, dec_c = np.genfromtxt(path+'/zfilter/compsUsed.csv', delimiter=',',usecols=(0,1),unpack=True)
#np.savetxt(path+"/zfilter/compsCoord.csv", np.column_stack((ra_c, dec_c)), delimiter=",", fmt='%0.8f')


#Fold comparison light curves to be shown in the FoV plots
os.chdir(path+'/gfilter')
os.system('astrosource_'+analysis+' --target-file $(pwd)/compsCoord.csv --indir $(pwd) --comparison --phot --plot --targetradius '+str(trad)+' --usescreenedcomps --usecompsused --usecompletedcalib --verbose')
os.chdir(path+'/rfilter')
os.system('astrosource_'+analysis+' --target-file $(pwd)/compsCoord.csv --indir $(pwd) --comparison --phot --plot --targetradius '+str(trad)+' --usescreenedcomps --usecompsused --usecompletedcalib --verbose')
os.chdir(path+'/ifilter')
os.system('astrosource_'+analysis+' --target-file $(pwd)/compsCoord.csv --indir $(pwd) --comparison --phot --plot --targetradius '+str(trad)+' --usescreenedcomps --usecompsused --usecompletedcalib --verbose')
#os.chdir(path+'/zfilter')
#os.system('astrosource_'+analysis+' --target-file $(pwd)/compsCoord.csv --indir $(pwd) --comparison --phot --plot --targetradius '+str(trad)+' --usescreenedcomps --usecompsused --usecompletedcalib --verbose')

os.environ["PATH"]=OLD_PATH
