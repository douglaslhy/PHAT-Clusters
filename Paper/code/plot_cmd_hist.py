import numpy as np
import pylab as plt
import os, sys, shutil
import pyfits
import matplotlib.patches as mpatches
from matplotlib import cm
from matplotlib import rc
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
from matplotlib.ticker import NullFormatter, MaxNLocator


'''plot cmd of all clusters, color-coded by which are bad fits, with histograms
Reads:  
by_eye/apdata_cluster_v2.fits
clust_results709.txt

Output:
plot_cmd_hist1.png
'''


def get_bins(start, stop, step):
    bins = np.arange(start, stop, step)
    return bins


def bin_width(bins):
    width = np.zeros(len(bins)-1)
    for i in range(len(bins)-1):
        width[i] = bins[i+1] - bins[i]
    return width



#read in clust results
results_table = 'clust_results709.txt'

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
C11_age_gyr = t[:,24]
C11_mass = t[:,25]
C11_EBV = t[:,26]
C11_Z = t[:,27]
C11_SigmaZ = t[:,28]

C11_age = np.log10(C11_age_gyr * 1e9)
C11_av = C11_EBV * 3.1


int_age_width = int_age_P84 - int_age_P16
int_mass_width = int_mass_P84 - int_mass_P16
int_av_width = int_av_P84 - int_av_P16
cmd_age_width = cmd_age_P84 - cmd_age_P16
cmd_mass_width = cmd_mass_P84 - cmd_mass_P16
cmd_av_width = cmd_av_P84 - cmd_av_P16

M = np.where(best == 'M')
I = np.where(best == 'I')
O = np.where(best == 'O')

# read in cluster photometry
phot_file = 'by_eye/apdata_cluster_v2.fits'

phot_file_open = pyfits.open(phot_file)
phot_data = phot_file_open[1].data
phot_id = phot_data.field('id')
mag_336 = phot_data.field('mag336')
mag_475 = phot_data.field('mag475')
mag_814 = phot_data.field('mag814')
mag_275 = phot_data.field('mag275')
mag_110 = phot_data.field('mag110')
mag_160 = phot_data.field('mag160')
cfrac = phot_data.field('newclstfrac')
 

#bad is for either NS or N classification
bad = np.where((cat == 'NS') | (cat == 'N'))
N = np.where(cat == 'N')
NS = np.where(cat == 'NS')
N_list = ap_id[N]
NS_list = ap_id[NS]



def plot_scatter_cmd_hist(ap_id, N_list, NS_list, xdata, ydata, xbin_min, xbin_max, ybin_min, ybin_max, xbin_width, ybin_width, xlabel, ylabel):

    '''create CMD with histograms along x and y axes'''

    plt.close('all')

    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left+width+0.02

    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]

    # start with a rectangular Figure
    plt.figure(1, figsize=(10,10))

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    axScatter = plt.axes(rect_scatter)
    axHistx = plt.axes(rect_histx)
    axHisty = plt.axes(rect_histy)

    # remove inner axes numbers of histograms
    nullfmt = NullFormatter()
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)

    axScatter.scatter(xdata+np.random.normal(0,0.01,len(xdata)), ydata, s=20, edgecolors='none',c='black',alpha=0.7)


    N_place=[]
    NS_place=[]
    for i in range(len(ap_id)):
        if ap_id[i] in N_list:
            axScatter.scatter(xdata[i]+np.random.normal(0,0.01), ydata[i], s=40, edgecolors='none', c='r')    # with noise
            N_place.append(i)
        elif ap_id[i] in NS_list:
            axScatter.scatter(xdata[i]+np.random.normal(0,0.01), ydata[i], s=40, edgecolors='none', c='green')    # with noise
            NS_place.append(i)


    axScatter.set_xlabel(xlabel, fontsize=20)
    axScatter.set_ylabel(ylabel, fontsize=20)

    axScatter.set_xlim([xbin_min, xbin_max])
    axScatter.set_ylim([ybin_max, ybin_min])  

    ticklabels = axScatter.get_xticklabels()
    for label in ticklabels:
        label.set_fontsize(18)

    ticklabels = axScatter.get_yticklabels()
    for label in ticklabels:
        label.set_fontsize(18)

    #histograms
    x_bins=get_bins(xbin_min, xbin_max, xbin_width)
    y_bins = get_bins(ybin_min, ybin_max, ybin_width)
    axHistx.hist(xdata, bins=x_bins, color='black')
    axHistx.hist(xdata[NS_place], bins=x_bins, color='green')
    axHistx.hist(xdata[N_place], bins=x_bins, color='red')

    axHisty.hist(ydata, bins=y_bins, orientation='horizontal', color='black')
    axHisty.hist(ydata[NS_place], bins=y_bins, orientation='horizontal', color='green')
    axHisty.hist(ydata[N_place], bins=y_bins, orientation='horizontal', color='red')


    axHistx.set_xlim([xbin_min, xbin_max])
    axHisty.set_ylim([ybin_max, ybin_min])

    ticklabels = axHistx.get_yticklabels()
    for label in ticklabels:
        label.set_fontsize(12)

    ticklabels = axHisty.get_xticklabels()
    for label in ticklabels:
        label.set_fontsize(12)

    axHisty.xaxis.set_major_locator(MaxNLocator(4))
    axHistx.yaxis.set_major_locator(MaxNLocator(4))


    plt.savefig('plot_cmd_hist1.png')
    return NS_place



xdata = mag_475 - mag_814
ydata = mag_814
xbin_min = -2
xbin_max = 4
ybin_min = 12
ybin_max = 24
xbin_width = 0.5
ybin_width = 0.5
xlabel = 'F475W - F814W'
ylabel = 'F814W'


NS_place = plot_scatter_cmd_hist(ap_id, N_list, NS_list, xdata, ydata, xbin_min, xbin_max, ybin_min, ybin_max, xbin_width, ybin_width, xlabel, ylabel)
