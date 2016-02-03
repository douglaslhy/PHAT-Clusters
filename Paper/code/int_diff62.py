import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

'''Difference (optimal - integrated) in the number between integrated and 'best' results

Reads:
'../clust_results709.txt'
Output:
int_diff62.png'''

def get_bins(start, stop, step):
    bins = np.arange(start, stop, step)
    return bins

def bin_width(bins):
    width = np.zeros(len(bins)-1)
    for i in range(len(bins)-1):
        width[i] = bins[i+1] - bins[i]
    return width

#calc bins
age_bins = get_bins(6.2,11,0.2)
age_start_bins = get_bins(6.2,10.7,0.2)
age_mid_bins = get_bins(6.3, 10.8,0.2)

#calc bins
mass_bins = get_bins(2, 5.8, 0.2)
mass_start_bins = get_bins(2, 5.6, 0.2)
mass_mid_bins = get_bins(2.1, 5.6, 0.2)

#read in clust results
results_table = '../clust_results709.txt'

t = np.genfromtxt(results_table)
ap_id = t[:,0].astype(int)
R_ap = t[:,1]
n_stars = t[:,2]
n_bg = t[:,3]
cat = np.genfromtxt(results_table, usecols=(4), dtype=None)
best = np.genfromtxt(results_table, usecols=(5), dtype=None)
int_age = t[:,6]
int_age_P16 = t[:,7]
int_age_P84 = t[:,8]
int_mass = t[:,9]
int_mass_P16 = t[:,10]
int_mass_P84 = t[:,11]
int_av = t[:,12]
int_av_P16 = t[:,13]
int_av_P84 = t[:,14]
cmd_age = t[:,15]
cmd_age_P16 = t[:,16]
cmd_age_P84 = t[:,17]
cmd_mass = t[:,18]
cmd_mass_P16 = t[:,19]
cmd_mass_P84 = t[:,20]
cmd_av = t[:,21]
cmd_av_P16 = t[:,22]
cmd_av_P84 = t[:,23]

clust_age=np.zeros(len(cat))
clust_mass=np.zeros(len(cat))
for i in range(len(cat)):
    if best[i] == 'M':
        clust_age[i] = cmd_age[i]
        clust_mass[i] = cmd_mass[i]
    elif best[i] == 'I':
        clust_age[i] = int_age[i]
        clust_mass[i] = int_mass[i]



plt.close('all')

h1,x1,y1=np.histogram2d(clust_mass, clust_age, bins=([mass_bins, age_bins]))
h2,x2,y2=np.histogram2d(int_mass, int_age, bins=([mass_bins, age_bins]))
diff=h1-h2
# scale by number of clust in each bin
scaled_diff=diff/np.sqrt((h1+h2)/2.)

plt.imshow(scaled_diff,interpolation='none',origin='lower', cmap=plt.cm.seismic,vmin=-7.112)
plt.xticks([0.5, 3.5, 8.5, 13.5, 18.5], ['6.4', '7.0', '8.0', '9.0', '10.0'], fontsize=20)
plt.yticks([0.5, 4.5, 9.5, 14.5], ['2.2', '3.0', '4.0', '5.0'], fontsize=20)
plt.xlabel('log($t$ [yr])', fontsize=25)
plt.ylabel(r'log($M$ [$M_{\odot}])', fontsize=25)
cbar = plt.colorbar(pad=0)
#cbar.set_label('Scaled Number Density', fontsize=18)
cbar.set_label('Fiducial - Integrated Scaled Number', fontsize=18)
cbar.ax.tick_params(labelsize=18)

plt.savefig('int_diff62.png')
