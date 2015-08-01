import numpy as np
import pylab as plt
import os, sys, shutil
import pyfits
import output
from matplotlib import cm

name = sys.argv[1]

'''plot cluster CMD and background CMD'''


#read in data from fits files
def readdata(name):
	''' read in phot, fake, and sky files '''

	phot_file = pyfits.open(name+ '/' +name+ '_phot.fits')
	phot_data = phot_file[1].data
	phot_f475 = phot_data.field('F475W_VEGA')
	phot_f814 = phot_data.field('F814W_VEGA')
	phot = np.zeros([len(phot_f475), 2])
	phot[:,0] = phot_f475
	phot[:,1] = phot_f814

	sky_file = pyfits.open(name + '/' +name+ '_sky.fits')
	sky_data = sky_file[1].data
	sky_f475 = sky_data.field('F475W_VEGA')
	sky_f814 = sky_data.field('F814W_VEGA')
	bg = np.zeros([len(sky_f475), 2])
	bg[:,0] = sky_f475
	bg[:,1] = sky_f814
	
	fake_file = pyfits.open(name+ '/' +name+ '.dst.fake.fits')
	fake_data = fake_file[1].data
	out_f475 = fake_data.field('F475W_VEGA')
	out_f814 = fake_data.field('F814W_VEGA')
	in_f475 = fake_data.field('F475W_IN')
	in_f814 = fake_data.field('F814W_IN')
	diff_f475 = rec_diff(out_f475, in_f475)
	diff_f814 = rec_diff(out_f814, in_f814)

	fake = np.zeros([len(out_f475), 4])
	fake[:,0] = out_f475
	fake[:,1] = out_f814
	fake[:,2] = diff_f475			#out - in
	fake[:,3] = diff_f814			#out - in

	return phot, fake, bg

def rec_diff(out_mag, in_mag):
	'''tests for recovery for asts
	outputs out-in, if recovered, 99.999 if not recovered'''

	out_in = out_mag - in_mag
	rec = np.abs(out_in) <= 0.75		#test for recovery
	diff = np.zeros(len(out_mag))
	diff[rec] = out_in[rec]			#recovered stars = out - in
	diff[~rec] = 99.999			#unrecovered stars = 99.999
	return diff

#find best fit values
def bestfit(matchoutput,mass):
	"""print values associated with minimum fit value
	1st argument = matchoutput (from apXXX/console_apXXX.txt)
	2nd argument = mass calculated from matchoutput"""
	fit = matchoutput[:,5]
	minfit = np.where(fit == fit.min())
	age = matchoutput[:,3]
	av = matchoutput[:,0]
	z = matchoutput[:,4]
	print age[minfit]
	print av[minfit]
	print mass[minfit]
	print z[minfit]
	return age[minfit], av[minfit], mass[minfit], z[minfit]

#plot cmds
def plotdata(phot, bg, bestage, av):
	"""plot cluster cmd data"""
	m475=phot[:,0]
	m814=phot[:,1]
	color=m475-m814
	ax0 = plt.subplot(121)
	plt.scatter(color, m475, edgecolor='None', facecolor='black', s=20)
	isocol, iso475 = plotiso(bestage, av)
	print isocol, iso475
	plt.xticks([np.arange(-0.25, 2.26, 0.25)])
	#plt.ylim(29.99, 17.01, )
	plt.ylim(29.99, 16.)
	plt.xlim(-0.49, 4.49)
	plt.xlabel('F475W-F814W', fontsize=16)
	plt.ylabel('F475W', fontsize=16)
	plt.title('(a) Cluster ' + name.upper(), fontsize=16)
	ax1 = plt.subplot(122, sharex=ax0, sharey=ax0)
	plt.scatter(bg[:,0]-bg[:,1], bg[:,0], edgecolor='None', facecolor='red', s=20)
	plt.setp(ax1.get_yticklabels(), visible=False)
	plt.ylim(27.99, 14.01)	
	plt.xlabel('F475W-F814W', fontsize=16)
	plt.title('(b) Background')

#overplot isochrones
def plotiso(bestage, av):
	"""find isochrone for the cluster's age from Match, and overplot it"""
	age = np.genfromtxt('iso_ages_solar.dat', skip_header=12, usecols=0, dtype=None)
	f475 = np.genfromtxt('iso_ages_solar.dat', skip_header=12, usecols=9, dtype=None)
	f814 = np.genfromtxt('iso_ages_solar.dat', skip_header=12, usecols=10, dtype=None)
	f475 = f475+24.47
	f814 = f814+24.47
	iso = np.where(age == bestage)
	#add reddening
	f475red, f814red = red(f475[iso], f814[iso], av)
	redcol = f475red - f814red
	plt.plot(redcol, f475red)
	return redcol, f475red

#calc reddening
def red(f475, f814, av):
	"""calculate reddened magnitudes"""
	#extinction values from Cardelli extinction law			
	R475 = 1.19119
	R814 = 0.60593
	#add extinction
	A475 = av * R475
	A814 = av * R814
	#calculate reddened mags
	f475red = f475 + A475
	f814red = f814 + A814
	return f475red, f814red


#read in data
phot, fake, bg = readdata(name)

#read in output from console file into matchoutput
matchoutput = np.genfromtxt(name+'/console_'+name+'.txt', skip_header=10, skip_footer=2)
#calculate weights
fit = matchoutput[:,5]	
fit = fit-min(fit)				#make minimum fit = 0
weights = np.exp(-0.5*fit)			#convert to weights

#calculate mass
mass = np.log10(matchoutput[:,6]*(10**(matchoutput[:,3]+0.1) - 10**matchoutput[:,3])+1)


#print values assoc. w/ lowest fit value
bestval=bestfit(matchoutput,mass)
fitlist=[]
for item in bestval:
	fitlist.append(float(item[0]))

bestage=fitlist[0]
#av=fitlist[1]
av=0.
print av


#plot cmd w/ reddened isochrone
plt.close('all')
plotdata(phot, bg, bestage, av)

plt.ylim(27.99, 16.)
plt.xlim(-0.49, 4.49)
plt.subplots_adjust(wspace=0.)
plt.savefig(name+ '/' +name+'_cmd_test2.png', dpi=150, bbox_inches='tight')

