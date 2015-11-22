import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import rc
import pyfits
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
from matplotlib.ticker import NullFormatter, MaxNLocator

'''plot scatter plots, with histograms along each axis, color-code flagged clusters

Reads in:
../clust_results709.txt
unres_results512.txt
../apdata_cluster_v2.fits

Outputs:
hist_flag611.png
plot_flag_ndet.png
'''

def get_bins(start, stop, step):
    bins = np.arange(start, stop, step)
    return bins


def bin_width(bins):
    width = np.zeros(len(bins)-1)
    for i in range(len(bins)-1):
        width[i] = bins[i+1] - bins[i]
    return width


def random_sample(N, p2, p16, p50, p84, p98):
    '''return median for N random samples of a pdf'''

    a = np.array([p2, p16, p50, p84, p98])
    prob = np.array([.02, .14, .68, .14, .02])
    rand_picks = np.random.choice(a, size=N, replace=True, p=prob) 
    return np.mean(rand_picks), np.median(rand_picks)



#read in clust results
t = np.genfromtxt('../clust_results709.txt')
clust_id = t[:,0].astype(int)
R_ap = t[:,1]
n_stars = t[:,2]
n_bg = t[:,3]
cat = np.genfromtxt('../clust_results709.txt', usecols=(4), dtype=None)
best = np.genfromtxt('../clust_results709.txt', usecols=(5), dtype=None)
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

flag = np.where(cat == 'N')
flag_list = clust_id[flag]
bad_match = np.where((cmd_age <= 8.5) & (best == 'I'))
bad_match_list = clust_id[bad_match]
bad_int = np.where((cmd_age > 8.5) & (best == 'M'))
bad_int_list = clust_id[bad_int]


int_age_width = int_age_P84 - int_age_P16
int_mass_width = int_mass_P84 - int_mass_P16
int_av_width = int_av_P84 - int_av_P16
cmd_age_width = cmd_age_P84 - cmd_age_P16
cmd_mass_width = cmd_mass_P84 - cmd_mass_P16
cmd_av_width = cmd_av_P84 - cmd_av_P16

M = np.where(best == 'M')
I = np.where(best == 'I')
O = np.where(best == 'O')


clust_age=np.zeros(len(cat))
clust_mass=np.zeros(len(cat))
clust_av=np.zeros(len(cat))
for i in range(len(cat)):
    if best[i] == 'M':
        clust_age[i] = cmd_age[i]
        clust_mass[i] = cmd_mass[i]
        clust_av[i] = cmd_av[i]
    elif best[i] == 'I':
        clust_age[i] = int_age[i]
        clust_mass[i] = int_mass[i]
        clust_av[i] = int_av[i]
    elif best[i] == 'O':
        clust_age[i] = C11_age[i]
        clust_mass[i] = C11_mass[i]
        clust_av[i] = C11_av[i]


#get unres results
unres=np.loadtxt('unres_results512.txt')
unres_id = unres[:,0]
int_mag_475 = unres[:,3]
int_mag_814 = unres[:,4]
int_color_48 = unres[:,9]
unres_mag_475 = unres[:,12]
unres_mag_814 = unres[:,13]
unres_color_48 = unres[:,18]
mag_diff_475_in = unres_mag_475 - int_mag_475
mag_diff_814_in = unres_mag_814 - int_mag_814
color_diff_48_in = np.abs(unres_color_48 - int_color_48)

#get apdata results
apdata_file=pyfits.open('../apdata_cluster_v2.fits')
apdata=apdata_file[1].data
ap_id  = apdata.field('ID')
ap_mag275 = apdata.field('MAG275')
ap_det_275 = apdata.field('DETSIG275')
ap_mag336 = apdata.field('MAG336')
ap_det_336 = apdata.field('DETSIG336')
ap_mag475 = apdata.field('MAG475')
ap_det_475 = apdata.field('DETSIG475')
ap_mag814 = apdata.field('MAG814')
ap_det_814 = apdata.field('DETSIG814')
ap_mag110 = apdata.field('MAG110')
ap_det_110 = apdata.field('DETSIG110')
ap_mag160 = apdata.field('MAG160')
ap_det_160 = apdata.field('DETSIG160')


#find diff for clust
mag_diff_475 = np.zeros(len(clust_id))
mag_diff_814 = np.zeros(len(clust_id))
color_diff_48 = np.zeros(len(clust_id))
for i in range(len(clust_id)):
    w=np.where(unres_id == clust_id[i])
    if len(w[0]) > 0:
        mag_diff_475[i] = mag_diff_475_in[w]
        mag_diff_814[i] = mag_diff_814_in[w]
        color_diff_48[i] = color_diff_48_in[w]




#find photometry for clust
det_275 = np.zeros(len(clust_id))
det_336 = np.zeros(len(clust_id))
det_475 = np.zeros(len(clust_id))
det_814 = np.zeros(len(clust_id))
det_110 = np.zeros(len(clust_id))
det_160 = np.zeros(len(clust_id))
n_det = np.zeros(len(clust_id))
for i in range(len(clust_id)):
    w2=np.where(ap_id == clust_id[i])
    if len(w2[0]) > 0:
        det_275[i] = ap_det_275[w2]
        if det_275[i] >= 3: 
            n_det[i] = n_det[i] + 1
        det_336[i] = ap_det_336[w2]
        if det_336[i] >= 3: 
            n_det[i] = n_det[i] + 1
        det_475[i] = ap_det_475[w2]
        if det_475[i] >= 3: 
            n_det[i] = n_det[i] + 1
        det_814[i] = ap_det_814[w2]
        if det_814[i] >= 3: 
            n_det[i] = n_det[i] + 1
        det_110[i] = ap_det_110[w2]
        if det_110[i] >= 3: 
            n_det[i] = n_det[i] + 1
        det_160[i] = ap_det_160[w2]
        if det_160[i] >= 3: 
            n_det[i] = n_det[i] + 1


def plot_scatter_hist(xdata, ydata, xbin_min, xbin_max, ybin_min, ybin_max, xbin_width, ybin_width, xlabel, ylabel, outfile, noise=0):

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

    #remove inner axes numbers of histograms
    nullfmt = NullFormatter()
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)

    if noise == 0:
        axScatter.scatter(xdata, ydata, s=20, edgecolors='none',c='black',alpha=0.7)

        flag_place=[]
        for i in range(len(clust_id)):
            if clust_id[i] in flag_list:
                axScatter.scatter(xdata[i], ydata[i], s=40, edgecolors='none', c='r', alpha=0.5)
                flag_place.append(i)

    else:
        #axScatter.scatter(xdata+np.random.normal(0,0.02,len(xdata)), ydata, s=20, edgecolors='none',c='black',alpha=0.7)
        axScatter.scatter(xdata+np.random.normal(0,0.02,len(xdata)), ydata+np.random.normal(0,0.05,len(ydata)), s=20, edgecolors='none',c='black',alpha=0.7)


        flag_place=[]
        bad_match_place=[]
        bad_int_place=[]
        wrong_place=[]
        for i in range(len(clust_id)):
            if clust_id[i] in flag_list:
                #axScatter.scatter(xdata[i]+np.random.normal(0,0.02), ydata[i], s=40, edgecolors='none', c='r')
                axScatter.scatter(xdata[i]+np.random.normal(0,0.02), ydata[i]+np.random.normal(0,0.05), s=40, edgecolors='none', c='r')
                flag_place.append(i)
            elif clust_id[i] in bad_match_list:
                axScatter.scatter(xdata[i], ydata[i]+np.random.normal(0,0.05), s=40, edgecolors='none', c='green')
                bad_match_place.append(i)
                wrong_place.append(i)
            elif clust_id[i] in bad_int_list:
                axScatter.scatter(xdata[i], ydata[i]+np.random.normal(0,0.05), s=40, edgecolors='none', c='green')
                bad_int_place.append(i)
                wrong_place.append(i)



    axScatter.set_xlabel(xlabel, fontsize=20)
    axScatter.set_ylabel(ylabel, fontsize=20)

    axScatter.set_xlim([xbin_min, xbin_max])
    axScatter.set_ylim([ybin_min, ybin_max])

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
    axHistx.hist(xdata[wrong_place], bins=x_bins, color='green')
    axHistx.hist(xdata[flag_place], bins=x_bins, color='red')
    #axHistx.hist(xdata[bad_match_place], bins=x_bins, color='green')
    #axHistx.hist(xdata[bad_int_place], bins=x_bins, color='green')

    axHisty.hist(ydata, bins=y_bins, orientation='horizontal', color='black')
    axHisty.hist(ydata[wrong_place], bins=y_bins, orientation='horizontal', color='green')
    axHisty.hist(ydata[flag_place], bins=y_bins, orientation='horizontal', color='red')
    #axHisty.hist(ydata[bad_match_place], bins=y_bins, orientation='horizontal', color='green')
    #axHisty.hist(ydata[bad_int_place], bins=y_bins, orientation='horizontal', color='green')

    axHistx.set_xlim([xbin_min, xbin_max])
    axHisty.set_ylim([ybin_min, ybin_max])

    ticklabels = axHistx.get_yticklabels()
    for label in ticklabels:
        label.set_fontsize(12)

    ticklabels = axHisty.get_xticklabels()
    for label in ticklabels:
        label.set_fontsize(12)

    axHisty.xaxis.set_major_locator(MaxNLocator(4))
    axHistx.yaxis.set_major_locator(MaxNLocator(4))


    #plt.show()
    plt.savefig(outfile)
    return flag_place


good = np.where(int_mass > 0)

# xdata = bf_age
# ydata = int_age

# xbin_min = 6.2
# xbin_max = 10.4
# ybin_min = xbin_min
# ybin_max = xbin_max
# xbin_width = 0.2
# ybin_width = 0.2

# xlabel = 'MATCH log(age [yr])'
# ylabel = 'Integrated log(age [yr])'


# xdata = bf_mass
# ydata = int_mass

# xbin_min = 1.8
# xbin_max = 6.0
# ybin_min = xbin_min
# ybin_max = xbin_max
# xbin_width = 0.2
# ybin_width = 0.2

# xlabel = 'MATCH log(mass [$M_{\odot}$])'
# ylabel = 'Integrated log(mass [$M_{\odot}$])'

# xdata = bf_age
# ydata = cfrac

# xbin_min = 6.2
# xbin_max = 10.4
# ybin_min = 0.62
# ybin_max = 1.01
# xbin_width = 0.2
# ybin_width = 0.02

# xlabel = 'MATCH log(age [yr])'
# ylabel = 'Cluster Frac'


# xdata = int_mass
# ydata = n_det

# xbin_min = 1.6
# xbin_max = 6.0
# ybin_min = -0.5
# ybin_max = 7
# xbin_width = 0.2
# ybin_width = 1.

#xlabel = 'Integrated log(M [$M_{\odot}$])'
#ylabel = 'Number of Filter Detections'


#cmd age vs int age, color-coded by unacceptable fits in both methods in red, and clusters where the non default fit is better are green

xdata = cmd_age
ydata = int_age

xbin_min = 6.2
xbin_max = 10.4
ybin_min = 6.2
ybin_max = 10.4
xbin_width = 0.2
ybin_width = 0.2

xlabel = 'MATCH log(t [yr])'
ylabel = 'Integrated log(t [yr])'

flag_place = plot_scatter_hist(xdata, ydata, xbin_min, xbin_max, ybin_min, ybin_max, xbin_width, ybin_width, xlabel, ylabel, outfile='hist_flag_new.png', noise=1)



#int mass vs number of filter detections, color-coded by unacceptable fits in both methods in red, and clusters where the non default fit is better are green


xdata = int_mass
ydata = n_det

xbin_min = 1.6
xbin_max = 6.0
ybin_min = -0.5
ybin_max = 7
xbin_width = 0.2
ybin_width = 1.

xlabel = 'Integrated log(M [$M_{\odot}$])'
ylabel = 'Number of Filter Detections'

flag_place = plot_scatter_hist(xdata, ydata, xbin_min, xbin_max, ybin_min, ybin_max, xbin_width, ybin_width, xlabel, ylabel, outfile='plot_flag_ndet_new.png', noise=1)



