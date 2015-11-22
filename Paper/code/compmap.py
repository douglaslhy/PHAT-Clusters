"""compmap
make completeness map plots for artificial star tests
inputs:  apXXX/apXXX.dst.fake.fits
outputs:  apXXX/apXXX_compmap.png, apXXX_475map.txt, apXXX_814map.txt
parameter:  name of cluster
example call:  >>> %run compmap.py ap1"""

import sys
import pylab as plt
import numpy as np
plt.ion()
import pdb
import pyfits
from scipy.interpolate import interp1d
import matplotlib.image as mpimg
from matplotlib import rc
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)



def detected(stars):
	"""compute percentage of total stars that are detected
	stars is array of magout-magin values"""

	n=len(stars)
	ud=np.where(abs(stars) > 0.75)				#finds undetected stars
	udmag=stars[ud]
	countud=len(udmag)					#counts number of undetected stars
	if n > 0:
		return float(n-countud) / float(n)		#gives percent detected
	else:
		return np.nan


def meandiff(stars):
	"""compute avg of mag difference of stars detected
	stars is array of magout-magin values"""

	n=len(stars)
	det=np.where(abs(stars) <= 0.75)				#finds detected stars
	detmag=stars[det]
	if len(detmag) > 0:
		return np.percentile(detmag, 50)			#gives 50% difference
	else:
		return np.nan

def stddiff(stars):
	"""compute std deviation of mag difference of stars detected
	stars is array of magout-magin values"""

	n=len(stars)
	det=np.where(abs(stars) <= 0.75)	       		#finds detected stars
	detmag=stars[det]
	if len(detmag) > 0:
		return detmag.std()				#gives std deviation of difference			
	else:
		return np.nan


def percentlow(stars):
	"""compute low percentage of mag difference of stars detected
	stars is array of magout-magin values"""

	n=len(stars)
	det=np.where(abs(stars) <= 0.75)	       		#finds detected stars
	detmag=stars[det]
	if len(detmag) > 0:
		return np.percentile(detmag, 16)		#gives 16th percentile of difference
	else:
		return np.nan

def percentup(stars):
	"""compute upper percentage of mag difference of stars detected
	stars is array of magout-magin values"""

	n=len(stars)
	det=np.where(abs(stars) <= 0.75)	       		#finds detected stars
	detmag=stars[det]
	if len(detmag) > 0:
		return np.percentile(detmag, 84)		#gives 84th percentile of difference	
	else:
		return np.nan


def calcbins(number, start, inc):

	"""calculate bins"""

	binstart=[start]*(number+1)				
	bins=[0]*(number+1)
	for i in range(len(binstart)):
		bins[i]=binstart[i]+(i*inc)
	return bins

def findbins(stars, bin1, bin2):
	
	"""finds objects within specified bin
	stars is array to bin over, bin1 is lower limit of bin, bin2 is upper limit of bin"""

	binstars = np.where((stars >= bin1) & (stars < bin2))
	
	return binstars

def findbins2d(dec, ra, decbin1, decbin2, rabin1, rabin2):
	
	"""finds objects within specified ra and dec bin
	ra and dec are arrays to bin over, bin1 is lower limit of bin, bin2 is upper limit of bin
	returns indices where ra, dec are within bin"""

	binstars = np.where((dec >= decbin1) & (dec < decbin2) & (ra >= rabin1) & (ra < rabin2))
	
	return binstars

def findcompmag(comp):

	"""calculates 50% completeness magnitude given array of completeness values"""

	magbins=calcbins(23, 19.25, 0.5)	
	f = interp1d(magbins,comp)			#interpolate over bins
	xnew = np.linspace(19.5,30,100)
	#solve for where f(xnew) == 0.5
	compdiff = abs(f(xnew) - 0.5)			#calculate 50% completeness magnitude
	a = np.where(compdiff == np.nanmin(compdiff))
	compmag = xnew[a]

	return np.median(compmag)

def findcomppercent(comp):

	"""calculates completeness percent at 24.47 mag given array of completeness values"""

	magbins=calcbins(23, 19.25, 0.5)	
	f = interp1d(magbins,comp)			#interpolate over bins
	#find comp percent at f(24.47)
	comppercent=f(24.47)

	return np.median(comppercent)

def calcdist(ra, dec):

	"calculate arcsec distance from ra and dec in degrees"
	radist=(max(ra)-min(ra))*np.cos(np.radians(np.mean(dec)))*3600.
	decdist= (max(dec)-min(dec))*3600.

	return radist, decdist


def plotcomp(name, completeness1, completeness2, magbins):
	
	"plot completeness vs magnitude"

	plt.close('all')
	ax=plt.subplot(111)
	ax.plot(magbins, completeness1)
	ax.plot(magbins, completeness2)
	plt.xlabel('Magnitude In')
	plt.ylabel('Completeness')
	plt.legend(("F475W", "F814W"))
	plt.savefig(name+'/'+name+'_completeness.png', dpi=150, bbox_inches='tight')

def plotdiff(name, magdiffavg, maglow, magup, magbins, band):
	
	"plot difference in input vs output magnitude vs magnitude"

	plt.close('all')
	ax=plt.subplot(111)
	yerrlow = [a-b for a,b in zip(magdiffavg, maglow)]
	yerrup = [a-b for a,b in zip(magup, magdiffavg)]
	plt.errorbar(magbins, magdiffavg, yerr=[yerrlow, yerrup])
	plt.xlabel('Magnitude In')
	plt.ylabel('Magnitude (Out - In)')
	plt.title(name+' '+band)
	plt.ylim(np.nanmin(maglow)-0.1, np.nanmax(magup)+0.1)
	plt.savefig(name+'/'+name+'_diff'+band+'.png', dpi=150, bbox_inches='tight')
	

def plotbias(name, colordiffavg, colorlow, colorup, midcbin):
	
	"plot difference in input vs output magnitude vs magnitude"

	plt.close('all')
	ax=plt.subplot(111)
	yerrlow = [a-b for a,b in zip(colordiffavg, colorlow)]
	yerrup = [a-b for a,b in zip(colorup, colordiffavg)]
	plt.errorbar(midcbin, colordiffavg, yerr=[yerrlow, yerrup])
	plt.xlabel('Color In')
	plt.ylabel('Color (Out - In)')
	plt.title(name+' Color Bias')
	plt.ylim(np.nanmin(colorlow)-0.1, np.nanmax(colorup)+0.1)
	plt.savefig(name+'/'+name+'_colorbias.png', dpi=150, bbox_inches='tight')


def plotmap(name, compmag, compperc, decbins, rabins, decdist, radist, nast):

	"plot completness mag as a function of ra and dec for each cluster"

	#plot comp mag
	plt.close('all')
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
        fig = plt.figure(figsize=(30,7))

	ax0 = plt.subplot(131)
	plt.imshow(compmag, interpolation='nearest', cmap=plt.cm.Blues)
	plt.ylim(0,40)
	plt.xlim(len(rabins),0)
	plt.xticks([50,25,0],['-'+str(round(radist/2, 2))+'"', 'center', '+'+str(round(radist/2, 2))+'"'], fontsize=18)
	plt.yticks([50*0.758,25*0.758,1],['-'+str(round(decdist/2, 2))+'"', 'center', '+'+str(round(decdist/2, 2))+'"'], fontsize=18)
	plt.title('50\% Completeness Mag', fontsize=24)
	cbar=plt.colorbar(shrink=0.9)


	#plot log(number of asts)
	ax1 = plt.subplot(132)
	plt.imshow(np.log10(nast),  interpolation='nearest', cmap=plt.cm.Blues)
	plt.ylim(0,40)
	plt.xlim(len(rabins),0)
	plt.xticks([50,25,0],['-'+str(round(radist/2, 2))+'"', 'center', '+'+str(round(radist/2, 2))+'"'], fontsize=18)
	plt.yticks([50*0.758,25*0.758,1],['-'+str(round(decdist/2, 2))+'"', 'center', '+'+str(round(decdist/2, 2))+'"'], fontsize=18)
	plt.title('log(number of ASTs)', fontsize=24)
	cbar=plt.colorbar(shrink=0.9)


	
	#plot completeness value for magnitude of 24.47
	ax2 = plt.subplot(133)
	plt.imshow(compperc, interpolation='nearest', cmap=plt.cm.Blues)
	plt.ylim(0,40)
	plt.xlim(len(rabins),0)
	plt.xticks([50,25,0],['-'+str(round(radist/2, 2))+'"', 'center', '+'+str(round(radist/2, 2))+'"'], fontsize=18)
	plt.yticks([50*0.758,25*0.758,1],['-'+str(round(decdist/2, 2))+'"', 'center', '+'+str(round(decdist/2, 2))+'"'], fontsize=18)
	plt.title('Completeness \% at 24.47 Mag', fontsize=24)
	cbar=plt.colorbar(shrink=0.9)



	#plot image of cluster
	#ax3 = plt.subplot(224)
	#plt.imshow(img[100:200,100:200], cmap=plt.cm.Blues)
	#plt.ylim(0,100)
	#plt.xticks(fontsize=10)
	#plt.yticks(fontsize=10)
	#plt.title(name, fontsize=12)

	#save map as .png
        #plt.subplots_adjust(left=0.2, right=0.8)
	plt.savefig(name+'/'+name+'_compmap.png', dpi=150, bbox_inches='tight')
	

name = sys.argv[1]

#read in fits data
hdulist = pyfits.open(name+'/'+name+'.dst.fake.fits')

#read in image of cluster
#img=mpimg.imread(name+'/'+name+'_img.jpg')

tbdata = hdulist[1].data
cols = hdulist[1].columns
mag1in = tbdata.field('F475W_IN')
mag1out = tbdata.field('F475W_VEGA')
mag2in = tbdata.field('F814W_IN')
mag2out = tbdata.field('F814W_VEGA')

mag1rec = mag1out - mag1in			#magout - magin
mag2rec = mag2out - mag2in			#magout - magin
rad = tbdata.field('RADIUSIN')
rad = 0.05*rad					#convert from pixels to arcsec
radeg = tbdata.field('RA_J2000')
decdeg = tbdata.field('DEC_J2000')
ra=(radeg*3600.)/0.05				#convert from deg to pixels
dec=(decdeg*3600.)/0.05				#convert from deg to pixels
radiff = max(ra)-min(ra)
decdiff = max(dec)-min(dec)

#make position bins
nbins=50
inc = radiff/nbins
rabins = calcbins(nbins, min(ra), inc)
decbins = calcbins(nbins, min(dec), inc)
nbinstot = len(rabins)*len(decbins)

compmag1=np.zeros(shape=(nbins,nbins))
compmag2=np.zeros(shape=(nbins,nbins))
compperc1=np.zeros(shape=(nbins,nbins))
compperc2=np.zeros(shape=(nbins,nbins))
nast=np.zeros(shape=(nbins,nbins))
detstars=np.zeros(shape=(nbins,nbins))
#compmag1=[]	
#compmag2=[]

#loop over all regions
for j in range(len(decbins)-1):
	for k in range(len(rabins)-1):

		#make region bins by calling findbins2d function
		regbin = findbins2d(dec, ra, decbins[j], decbins[j+1], rabins[k], rabins[k+1])		
		mag1inreg = mag1in[regbin]
		mag2inreg = mag2in[regbin]
		mag1outreg = mag1rec[regbin]
		mag2outreg = mag2rec[regbin]
		nast[j,k] = len(mag1inreg)		
	
		bins = calcbins(24, 19, 0.5)			#calculate mag bins
		midbin=calcbins(23, 19.25, 0.5)
		completeness1=[]
		completeness2=[]
		mag1outavg=[]
		mag2outavg=[]
		mag1outstd=[]
		mag2outstd=[]
		mag1low=[]
		mag1up=[]
		mag2low=[]
		mag2up=[]

		#loop over all magnitude bins for each region
		for i in range(len(bins)-1):
			magbin1 = findbins(mag1inreg, bins[i], bins[i+1])	#make mag bins by calling findbins function
			magbin2 = findbins(mag2inreg, bins[i], bins[i+1])			
			mag1inbin = mag1inreg[magbin1]
			mag2inbin = mag2inreg[magbin2]
			mag1outbin = mag1outreg[magbin1]
			mag2outbin = mag2outreg[magbin2]
			completeness1.append(detected(mag1outbin))		#find completeness for each mag bin
			completeness2.append(detected(mag2outbin))		

			mag1outavg.append(meandiff(mag1outbin))			#calculate average diff for bin	
			mag2outavg.append(meandiff(mag2outbin))
			mag1outstd.append(stddiff(mag1outbin))			#calculate std dev of diff for bin
			mag2outstd.append(stddiff(mag2outbin))

			mag1low.append(percentlow(mag1outbin))	
			mag1up.append(percentup(mag1outbin))
			mag2low.append(percentlow(mag2outbin))
			mag2up.append(percentup(mag2outbin))

		#find mag at 50% completeness for each region
		compmag1[j,k]=findcompmag(completeness1)
		compmag2[j,k]=findcompmag(completeness2)
		compperc1[j,k]=findcomppercent(completeness1)
		compperc2[j,k]=findcomppercent(completeness2)
		detstars[j,k] = detected(mag1outreg)

		#print mag1inreg, completeness1

#save these completeness mags
np.savetxt(name+'/'+name+'_475map.txt', compmag1)
np.savetxt(name+'/'+name+'_814map.txt', compmag2)
radegbin = [(0.05*(r/3600.)) for r in rabins]
decdegbin = [(0.05*(d/3600.)) for d in decbins]
incdeg = 0.05*(inc/3600.)
#calculate distances
radist, decdist=calcdist(radeg, decdeg)
#plotmap(name, compmag1, compperc1, img, decbins, rabins, decdist, radist, nast)

plotmap(name, compmag1, compperc1, decbins, rabins, decdist, radist, nast)

