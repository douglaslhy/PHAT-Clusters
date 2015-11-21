# !/usr/bin/python
""" 
run_match
Runs MATCH version 2.5 on chex

Runs AP clusters through the MATCH CMD fitting program, outputs fitting results and several analysis plots
Example Run Command: python run_match.py ap1 0 &
1st argument:  cluster name
2nd argument:  dav value (0 to 0.5)

run from one directory up from apXXX
Inputs: writeparam.py; output.py; write_fits.py; pg.py; iso_ages_solar.dat; apXXX/apXXX_phot.fits; apXXX/apXXX.dst.fake.fits; apXXX/apXXX_sky.fits 
Outputs: apXXX/apXXX_cmd.png; apXXX/param.sfh; apXXX/apXXX_out.cmd; apXXX/console_apXXX.txt; apXXX/apXXX_out; apXXX/apXXX_pdf.png; apXXX/apXXX_jointpdf_av_age.png
"""

import numpy as np
from matplotlib import use
use('Agg')
import pylab as plt
import os, sys, shutil
import pyfits
import output
from matplotlib import cm
import writeparam
import write_fits
import pg



def readdata(name):
	''' read in cluster photometry, background photometry, and artificial star tests from fits files
        input:  cluster name, e.g. 'ap1'
        outputs:  array of cluster photometry (F475, F814);
        array of background photometry (F475, F814);
        array of artificial star tests (F475 out, F814 out, F475 out-in, F814 out-in) '''

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
        inputs:  output magnitude, input magnitude
	output: out-in if star is  recovered, 99.999 if star is not recovered'''

	out_in = out_mag - in_mag
	rec = np.abs(out_in) <= 0.75		#test for recovery
	diff = np.zeros(len(out_mag))
	diff[rec] = out_in[rec]			#recovered stars = out - in
	diff[~rec] = 99.999			#unrecovered stars = 99.999
	return diff
	

def make_dat(name, phot, fake, bg):
	'''make photmetry, background, and artifical star dat files to be read by MATCH
        inputs:  cluster name; photometry array, artificial star array, background array
        outpus:  apXXX/apXXX_phot.dat; apXXX/apXXX_fake.dat; apXXX/bg.dat'''
	np.savetxt(name+ '/' +name+ '_phot.dat', phot, fmt='%f')
	np.savetxt(name+ '/' +name+ '_fake.dat', fake, fmt='%f')
	np.savetxt(name+ '/bg.dat', bg, fmt='%f')


def completeness(fake):
	"""  calculate completeness from artificial star tests
        input:  array of artificial stars  (F475 out, F814 out, F475 out-in, F814 out-in) 
        outputs:  array of completeness (completeness for F475, completeness for F814)"""
	m1, m2, dm1, dm2 = fake.T
	mag1range = np.arange(m1.min(), m1.max(), 0.25)
	mag2range = np.arange(m2.min(), m2.max(), 0.25)
	comp1 = np.zeros(len(mag1range))
	comp2 = np.zeros(len(mag2range))
	#assert()
	for i in range(1, len(mag1range)):
		w = np.where((m1 >= mag1range[i-1]) & (m1 < mag1range[i]))
		wbad = np.where(dm1[w] < 90)
		comp1[i] = np.float(len(wbad[0]))/len(w[0])
	for i in range(1, len(mag2range)):
		w = np.where((m2 >= mag2range[i-1]) & (m2 < mag2range[i]))
		wbad = np.where(dm2[w] < 90)
		comp2[i] = np.float(len(wbad[0]))/len(w[0])

	comp50m1 = np.abs(comp1 - complim).argmin()
	comp50m2 = np.abs(comp2 - complim).argmin()
	return mag1range[comp50m1], mag2range[comp50m2]


def bestfit(matchoutput, mass):
	"""print values associated with minimum fit value
	inputs:  matchoutput (from apXXX/console_apXXX.txt); mass calculated from matchoutput
        outputs:  best fit age, av, mass, z"""
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


def plotdata(phot, bg, bestage, av):
	"""plot cluster cmd data
        inputs:  array of cluster photometry; array of background photometry; best fit age; best fit av """
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
	plt.xlim(-0.49, 3.49)
	plt.xlabel('F475W-F814W', fontsize=16)
	plt.ylabel('F475W', fontsize=16)
	plt.title('(a) Cluster ' + name.upper(), fontsize=16)
	ax1 = plt.subplot(122, sharex=ax0, sharey=ax0)
	plt.scatter(bg[:,0]-bg[:,1], bg[:,0], edgecolor='None', facecolor='red', s=20)
	plt.setp(ax1.get_yticklabels(), visible=False)
	plt.ylim(27.99, 14.01)	
	plt.xlabel('F475W-F814W', fontsize=16)
	plt.title('(b) Background')


def plotiso(bestage, av):
	"""find isochrone for the cluster's age from iso_ages_solar.dat, and overplot it
        inputs:  age, av
        outputs:  reddening array; reddened f475 magnitude"""
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
	"""calculate reddened magnitudes using Cardelli extinction law
        inputs:  f475 magnitude, f814 magnitude, av
        outputs:  reddened f475 magnitude, reddened f814 magnitude"""
	R475 = 1.19119
	R814 = 0.60593
	#add extinction
	A475 = av * R475
	A814 = av * R814
	#calculate reddened mags
	f475red = f475 + A475
	f814red = f814 + A814
	return f475red, f814red


def plotpdf(name, data, weights, label, xmin, xmax, bestfit, bdelta=0.1):
	"""  make pdf plots for age, av, and z
	inputs:  cluster name; array of parameters (age, av, or z) from matchoutput file; fit values from matchoutput file; 
        string name of parameter; minimum value of parameter to plot; maximum value of parameter to plot; best fit value of parameter; bin size
        output:  apXXX/apXXX_XXXptiles.txt"""

   	hdelta=bdelta/2.
  	nbins=np.around((data.max()-data.min())/bdelta)+1
  	xbins=np.linspace(np.around(data.min()-hdelta,3),np.around(data.max()+bdelta+hdelta,3),nbins+2)	#calculate bins
  	h = np.histogram(data, bins=xbins, weights=weights)						#create histogram of data
  	hist = h[0]/np.sum(h[0])									#normalize histogram
  	ptiles = output.weighted_percentile(data, weights, [0.16, 0.5, 0.84])				#compute percentiles
	np.savetxt(name+'/'+name+'_'+label+'ptiles.txt', ptiles)					#save these weighted percentiles
	plt.axvspan(ptiles[0], ptiles[2], color='0.8', alpha=0.5)
	plt.step(np.linspace(data.min(), data.max(), len(h[0])), hist, color='black', lw=3)
	plt.axvline(ptiles[1], color='r', lw=2)
	plt.axvline(bestfit, color='b', lw=2)
	plt.ylim(np.min(hist.min(), 0), hist.max()*1.05)
	plt.xlim(xmin, xmax)
	plt.xlabel(label)


def plotmass(name, data, weights, label, bestfit, bdelta=0.1):
	"""  make pdf plot for mass -- x-axis is scaled by mass data rather than being an input
	inputs:  cluster name; mass array from matchoutput file; fit values from matchoutput file; 
        string name of parameter; best fit value of parameter; bin size
        output:  apXXX/apXXX_XXXptiles.txt""" 

   	hdelta=bdelta/2.
  	nbins=np.around((data.max()-data.min())/bdelta)+1
  	xbins=np.linspace(np.around(data.min()-hdelta,3),np.around(data.max()+bdelta+hdelta,3),nbins+2)	#calculate bins
  	h = np.histogram(data, bins=xbins, weights=weights)						#create histogram of data
  	hist = h[0]/np.sum(h[0])									#normalize histogram
  	ptiles = output.weighted_percentile(data, weights, [0.16, 0.5, 0.84])				#compute percentiles
	np.savetxt(name+'/'+name+'_'+label+'ptiles.txt', ptiles)					#save these weighted percentiles
	hist = h[0]/np.sum(h[0])
	plt.axvspan(ptiles[0], ptiles[2], color='0.8', alpha=0.5)
	plt.step(np.linspace(data.min(), data.max(), len(h[0])), hist, color='black', lw=3)
	plt.axvline(ptiles[1], color='r', lw=2)
	plt.axvline(bestfit, color='b', lw=2)
	plt.ylim(np.min(hist.min(), 0), hist.max()*1.05)
	plt.xlabel(label)
	plt.xlim(data.min()*0.95, data.max()*1.05)


def jointclusterpdfs(name, matchoutput, weights):
	""" plots joint cluster pdf showing Av vs log age, color-coded by weight
	inputs:  cluster name; matchoutput (from apXXX/console_apXXX.txt); 
        weights calculated above from matchoutput (from apXXX/console_apXXX.txt)
        output:  apXXX/apXXX_jointpdf_av_age.png"""
	cmap = cm.YlGnBu_r
	#cmap.set_gamma(0.5)

	sigma = np.exp(-0.5*(matchoutput[:,5].min() + np.arange(0, 4, 1)**2))

	plt.close('all')
	plt.figure(2, figsize=(8,8))

	w3s = np.where(weights >= np.exp(-0.5*(matchoutput[:,5].min()+9)))

	plt.scatter(matchoutput[:,0], matchoutput[:,3], color='0.7', s=15)
	plt.scatter(matchoutput[:,0][w3s], matchoutput[:,3][w3s], c=weights[w3s], cmap=cmap, edgecolor='None', s=55)
	cb  = plt.colorbar(ticks=sigma)
	#cb.set_ticklabels(('', '1-$\sigma$', '2-$\sigma$', '3-$\sigma$'), update_ticks=True)
	cb.ax.set_yticklabels(('', '68$\%$', '95$\%$', '99.7$\%$'))
	plt.xlim(matchoutput[:,0].min() - 0.1, matchoutput[:,0].max()+0.1)
	plt.ylim(matchoutput[:,3].min() - 0.1, matchoutput[:,3].max()+0.1)
	plt.xlabel('Av', fontsize=20)
	plt.ylabel('Log Age', fontsize=20)
	plt.title(name)
	plt.savefig(name+'/'+name+'_jointpdf_av_age.png', dpi=150, bbox_inches='tight')



#cluster name, e.g. 'ap1'
name = sys.argv[1]

#differential extinction, usually set to 0
dav = sys.argv[2]


#pre-processing steps

#read in cluster photometry, background photometry, and artificial star tests
phot, fake, bg = readdata(name)

#create dat files for MATCH
make_dat(name, phot, fake, bg)

#create param file
orig_stdout = sys.stdout
f = file(name+'/param.sfh','w')
sys.stdout = f
writeparam.main(name)
sys.stdout = orig_stdout
f.close()



#run match 
script = '/usr/local/bin/match2.5/bin/calcsfh {0}/param.sfh {0}/{0}_phot.dat {0}/{0}_fake.dat {0}/{0}_out -full -dav=0.{1} -davy=0 -ssp > {0}/console_{0}.txt'.format(name, dav)
os.system(script)



#post-processing analysis/plots

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
av=fitlist[1]
#save these best fit values
np.savetxt(name+'/'+name+'_bestfit.txt', fitlist)

#plot cmd w/ reddened isochrone
plt.close('all')
plotdata(phot, bg, bestage, av)

plt.ylim(27.99, 16.)
plt.xlim(-0.49, 3.49)
plt.subplots_adjust(wspace=0.)
plt.savefig(name+'/'+name+'_cmd.png', dpi=150, bbox_inches='tight')

#make pdf plots and save to file
plt.close('all')
ax0 = plt.subplot(221)
plotpdf(name, matchoutput[:,3], weights, 'Log_Age', matchoutput[:,3].min()*0.95, matchoutput[:,3].max()*1.05, fitlist[0], bdelta=0.1)
ax1 = plt.subplot(222)
plotpdf(name, matchoutput[:,0], weights, 'Av', -0.05, matchoutput[:,0].max()+0.05, fitlist[1], bdelta=0.05)
ax2 = plt.subplot(223)
plotpdf(name, matchoutput[:,4], weights, 'Log_Z', -0.2, 0.2, fitlist[3])
ax3 = plt.subplot(224)
plotmass(name, mass, weights, 'Log_M', fitlist[2], bdelta=0.05)
plt.savefig(name+'/'+name+'_pdf.png', dpi=150, bbox_inches='tight')

#make joint pdf plot
jointclusterpdfs(name, matchoutput, weights)

#make residual plot
pg.main(name)

#write fits table
write_fits.main(name)




