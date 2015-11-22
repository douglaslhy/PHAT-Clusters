import numpy as np
import pylab as plt
import os, sys, shutil
import pyfits
from matplotlib import cm
from matplotlib import rc
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)


'''plots CMD along with best fit Match isochrone and expectation value integrated isochrone
reads phot and sky files from each ap directory
reads results from int_match_results_54.txt
'''



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



#plot cmds
def plotdata(phot, bg, bestage, av, int_age, int_av, match_width, int_width):
	"""plot cluster cmd data"""
	m475=phot[:,0]
	m814=phot[:,1]
	color=m475-m814
	ax0 = plt.subplot(121)
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')
	plt.scatter(color, m475, edgecolor='None', facecolor='black', s=20)
	isocol, iso475 = plotiso(bestage, av, match_width)
        int_isocol, int_iso475 = plotiso(int_age, int_av, int_width)
	#print isocol, iso475
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
def plotiso(bestage, av, width):
	"""find isochrone for the cluster's age from Match, and overplot it"""
	age = np.genfromtxt('iso_ages_solar.dat', skip_header=12, usecols=0, dtype=None)
	f475 = np.genfromtxt('iso_ages_solar.dat', skip_header=12, usecols=9, dtype=None)
	f814 = np.genfromtxt('iso_ages_solar.dat', skip_header=12, usecols=10, dtype=None)
	f475 = f475+24.47
	f814 = f814+24.47
        isodiff = np.abs(age - bestage)
        iso = np.where(np.min(isodiff) == isodiff)
	#iso = np.where(age == bestage)
	#add reddening
	f475red, f814red = red(f475[iso], f814[iso], av)
	redcol = f475red - f814red
	plt.plot(redcol, f475red, alpha=0.5, lw=width)
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



#read in fit results
d=np.genfromtxt('int_match_results_54.txt')

clust_id=d[:,0].astype('int')
int_age=d[:,1]
int_age_low=d[:,2]
int_age_high=d[:,3]
int_mass=d[:,4]
int_mass_low=d[:,5]
int_mass_high=d[:,6]
int_av=d[:,7]
int_av_low=d[:,8]
int_av_high=d[:,9]
bf_age=d[:,10]
age=d[:,12]
age_low=d[:,11]
age_high=d[:,13]
bf_mass=d[:,14]
mass=d[:,16]
mass_low=d[:,15]
mass_high=d[:,17]
bf_av=d[:,18]
av=d[:,20]
av_low=d[:,19]
av_high=d[:,21]
nstars=d[:,22]


name = 'ap24'

#find place in result arrays where clust_id == name
name_id = int(name.split('ap')[1])
place = np.where(clust_id == name_id)

#find default values for width of isochrone
match_width = np.zeros(len(bf_age))
int_width = np.zeros(len(bf_age))
young = bf_age <= 8.5
match_width[young] = 4
match_width[~young] = 2
int_width[young] = 2
int_width[~young] = 4


#read in photometry
phot, fake, bg = readdata(name)


#plot cmd w/ reddened isochrone
plt.close('all')
plotdata(phot, bg, bf_age[place], bf_av[place], int_age[place], int_av[place], match_width[place], int_width[place])

plt.ylim(27.99, 16.)
plt.xlim(-0.49, 4.49)
plt.subplots_adjust(wspace=0.)
plt.savefig(name+ '/' +name+'_cmd_test.png', dpi=150, bbox_inches='tight')

