import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
import scipy.optimize as so
import matplotlib.cm as cm
import pyfits
from matplotlib import rc
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

'''comparison plots for syn clusters using synindex_v3.fits and synindex_6band.dat
reads fake cluster results from 'B17-F01-' +str(batch[i])+'/' +clust+ '/' +clust+ '_bestfit.txt'
and 'B06-F11-' +str(batch[i])+'/' +clust+ '/' +clust+ '_bestfit.txt' 
'''


#read in input data
input = pyfits.open('synindex_v3.fits')
input_data = input[1].data
synid1 = input_data.field('synid')
batch1 = input_data.field('syniter')
name1 = input_data.field('name')
in_age1 = input_data.field('logage')
in_mass1 = input_data.field('logmass')
in_av1 = input_data.field('av')

synid2 = np.loadtxt('synindex_6band.dat', usecols=[0]).astype(int)
batch2 = np.loadtxt('synindex_6band.dat', usecols=[1]).astype(int)
name2 = np.genfromtxt('synindex_6band.dat', usecols=(2), dtype=None)
in_age2 = np.loadtxt('synindex_6band.dat', usecols=[7]).astype(float)
in_mass2 = np.loadtxt('synindex_6band.dat', usecols=[8]).astype(float)
in_av2 = np.loadtxt('synindex_6band.dat', usecols=[9]).astype(float)

#concatenate arrays
synid = np.concatenate([synid1, synid2])
batch = np.concatenate([batch1, batch2])
name = np.concatenate([name1, name2])
in_age = np.concatenate([in_age1, in_age2])
in_mass = np.concatenate([in_mass1, in_mass2])
in_av = np.concatenate([in_av1, in_av2])

#create arrays for match results
bf_age = np.zeros(len(name))
bf_mass = np.zeros(len(name))
bf_av = np.zeros(len(name))


#loop over input data to find match results
field=np.zeros(len(name)).astype('str')

for i in range(len(name)):
    field[i] = name[i].split('-')[1]

    try:
        if field[i] == 'F01':
            clust = name[i].split('_')[0]
            best_fit = np.loadtxt('B17-F01-' +str(batch[i])+'/' +clust+ '/' +clust+ '_bestfit.txt')

            bf_age[i] = best_fit[0]
            bf_av[i] = best_fit[1]
            bf_mass[i] = best_fit[2]

        elif field[i] == 'F11':

            clust = name[i].split('_')[0]
            best_fit = np.loadtxt('B06-F11-' +str(batch[i])+'/' +clust+ '/' +clust+ '_bestfit.txt')

            bf_age[i] = best_fit[0]
            bf_av[i] = best_fit[1]
            bf_mass[i] = best_fit[2]

    except IOError:
        print str(name[i])+ '_bestfit.txt does not exist'


#def plot_one_param():

good = np.where(bf_mass > 0)
in_age=in_age[good]
bf_age=bf_age[good]
in_mass=in_mass[good]
bf_mass=bf_mass[good]
in_av=in_av[good]
bf_av=bf_av[good]



plt.close('all')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

#make mass bins for color-code
bin1 = np.where((in_mass > 0) & (in_mass < 3))	 	#purple
bin2 = np.where((in_mass >= 3) & (in_mass < 3.5)) 	#blue
bin3 = np.where((in_mass >= 3.5) & (in_mass < 4)) 	#green
bin4 = np.where((in_mass >= 4) & (in_mass < 5)) 	#red
bin5 = np.where(in_mass >= 5)				#orange

plt.scatter(in_av[bin1] + np.random.normal(0, 0.03, len(in_av[bin1])), (in_av[bin1] - bf_av[bin1]),  s=15., edgecolors='none', alpha=0.5, color='purple')
plt.scatter(in_av[bin2] + np.random.normal(0, 0.03, len(in_av[bin2])), (in_av[bin2] - bf_av[bin2]),  s=15., edgecolors='none', alpha=0.5, color='blue')
plt.scatter(in_av[bin3] + np.random.normal(0, 0.03, len(in_av[bin3])), (in_av[bin3] - bf_av[bin3]),  s=15., edgecolors='none', alpha=0.5, color='green')
plt.scatter(in_av[bin4] + np.random.normal(0, 0.03, len(in_av[bin4])), (in_av[bin4] - bf_av[bin4]),  s=15., edgecolors='none', alpha=0.5, color='red')
plt.axhline(y=0, color='black', lw=2)
plt.axhline(y=1, color='black', ls='dashed', lw=2)
plt.axhline(y=-1, color='black', ls='dashed', lw=2)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlabel('Input $A_{V}$ [mag]', fontsize=20)
plt.ylabel('Input $A_{V}$ [mag] - MATCH $A_{V}$ [mag]', fontsize=20)
plt.xlim([0.2, 2.2])
plt.ylim([-2.6,2.6])


plt.close('all')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.scatter(in_mass[bin1] + np.random.normal(0, 0.03, len(in_mass[bin1])), (in_mass[bin1] - bf_mass[bin1]),  s=15., edgecolors='none', alpha=0.5, color='purple')
plt.scatter(in_mass[bin2] + np.random.normal(0, 0.03, len(in_mass[bin2])), (in_mass[bin2] - bf_mass[bin2]),  s=15., edgecolors='none', alpha=0.5, color='blue')
plt.scatter(in_mass[bin3] + np.random.normal(0, 0.03, len(in_mass[bin3])), (in_mass[bin3] - bf_mass[bin3]),  s=15., edgecolors='none', alpha=0.5, color='green')
plt.scatter(in_mass[bin4] + np.random.normal(0, 0.03, len(in_mass[bin4])), (in_mass[bin4] - bf_mass[bin4]),  s=15., edgecolors='none', alpha=0.5, color='red')
plt.axhline(y=0, color='black', lw=2)
plt.axhline(y=1, color='black', ls='dashed', lw=2)
plt.axhline(y=-1, color='black', ls='dashed', lw=2)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlabel('Input $\log_{10}$(M [$M_{\odot}$])', fontsize=20)
plt.ylabel('Input $\log_{10}$(M [$M_{\odot}$]) - MATCH $\log_{10}$(M [$M_{\odot}$])', fontsize=20)
plt.xlim([2.2,4.2])
plt.ylim([-2.6,2.6])


plt.close('all')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.scatter(in_age[bin1] + np.random.normal(0, 0.03, len(in_age[bin1])), (in_age[bin1] - bf_age[bin1])+ np.random.normal(0, 0.02, len(in_age[bin1])),  s=15., edgecolors='none', alpha=0.5, color='purple')
plt.scatter(in_age[bin2] + np.random.normal(0, 0.03, len(in_age[bin2])), (in_age[bin2] - bf_age[bin2])+ np.random.normal(0, 0.02, len(in_age[bin2])),  s=15., edgecolors='none', alpha=0.5, color='blue')
plt.scatter(in_age[bin3] + np.random.normal(0, 0.03, len(in_age[bin3])), (in_age[bin3] - bf_age[bin3])+ np.random.normal(0, 0.02, len(in_age[bin3])),  s=15., edgecolors='none', alpha=0.5, color='green')
plt.scatter(in_age[bin4] + np.random.normal(0, 0.03, len(in_age[bin4])), (in_age[bin4] - bf_age[bin4])+ np.random.normal(0, 0.02, len(in_age[bin4])),  s=15., edgecolors='none', alpha=0.5, color='red')
plt.axhline(y=0, color='black', lw=2)
plt.axhline(y=1, color='black', ls='dashed', lw=2)
plt.axhline(y=-1, color='black', ls='dashed', lw=2)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlabel('Input $\log_{10}$(t [yr])', fontsize=20)
plt.ylabel('Input $\log_{10}$(t [yr]) - MATCH $\log_{10}$(t [yr])', fontsize=20)
plt.xlim([6.3, 10.0])
plt.ylim([-3,3])


#mass vs av
plt.close('all')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')


plt.axhline(y=0, color='black')
plt.axhline(y=1, color='black', ls='dashed')
plt.axhline(y=-1, color='black', ls='dashed')
plt.axvline(x=0, color='black')
plt.axvline(x=1, color='black', ls='dashed')
plt.axvline(x=-1, color='black', ls='dashed')
plt.scatter((in_mass - bf_mass), (in_av - bf_av) + np.random.normal(0, 0.01, len(in_age)), s=15., edgecolors='none', alpha=0.5, c=in_age - bf_age)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlabel('Input $\log_{10}$(M [$M_{\odot}$]) - MATCH $\log_{10}$(M [$M_{\odot}$])', fontsize=20)
plt.ylabel('Input $A_{V}$ [mag] - MATCH $A_{V}$ [mag]', fontsize=20)
plt.xlim([-2.6, 2.6])
plt.ylim([-2.2, 2.2])
cbar = plt.colorbar(pad=0)
cbar.set_label('Input $\log_{10}$(t [yr]) - MATCH $\log_{10}$(t [yr])', fontsize=18)
#cbar.ax.tick_params(labelsize=15)
for t in cbar.ax.get_yticklabels():
    t.set_fontsize(15)



#age vs mass

plt.close('all')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')


plt.axhline(y=0, color='black')
plt.axhline(y=1, color='black', ls='dashed')
plt.axhline(y=-1, color='black', ls='dashed')
plt.axvline(x=0, color='black')
plt.axvline(x=1, color='black', ls='dashed')
plt.axvline(x=-1, color='black', ls='dashed')
plt.scatter((in_age - bf_age) + np.random.normal(0, 0.01, len(in_age)), (in_mass - bf_mass), s=15., edgecolors='none', alpha=0.5, c=in_av - bf_av)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlabel('Input $\log_{10}$(t [yr]) - MATCH $\log_{10}$(t [yr])', fontsize=20)
plt.ylabel('Input $\log_{10}$(M [$M_{\odot}$]) - MATCH $\log_{10}$(M [$M_{\odot}$])', fontsize=20)
plt.xlim([-2.6, 2.6])
plt.ylim([-2.6, 2.6])
cbar = plt.colorbar(pad=0)
cbar.set_label('Input $A_{V}$ [mag] - MATCH $A_{V}$ [mag]', fontsize=18)
#cbar.ax.tick_params(labelsize=15)
for t in cbar.ax.get_yticklabels():
    t.set_fontsize(15)





#age vs av

plt.close('all')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')


plt.axhline(y=0, color='black')
plt.axhline(y=1, color='black', ls='dashed')
plt.axhline(y=-1, color='black', ls='dashed')
plt.axvline(x=0, color='black')
plt.axvline(x=1, color='black', ls='dashed')
plt.axvline(x=-1, color='black', ls='dashed')
plt.scatter((in_age - bf_age) + np.random.normal(0, 0.01, len(in_age)), (in_av - bf_av), s=15., edgecolors='none', alpha=0.5, c=in_mass - bf_mass)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlabel('Input $\log_{10}$(t [yr]) - MATCH $\log_{10}$(t [yr])', fontsize=20)
plt.ylabel('Input $A_{V}$ [mag] - MATCH $A_{V}$ [mag]', fontsize=20)
plt.xlim([-2.6, 2.6])
plt.ylim([-2.6, 2.6])
cbar = plt.colorbar(pad=0)
cbar.set_label('Input $\log_{10}$(M [$M_{\odot}$]) - MATCH $\log_{10}$(M [$M_{\odot}$])', fontsize=18)
#cbar.ax.tick_params(labelsize=15)
for t in cbar.ax.get_yticklabels():
    t.set_fontsize(15)

plt.show()














