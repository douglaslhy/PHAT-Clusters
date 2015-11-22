#import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)


'''plots match vs int age results w/ uncertainties as a 2d histogram

reads results from ../data/clust_results709.txt

Output:  plot_int_match_age.png
'''


def get_bins(start, stop, step):
    bins = np.arange(start, stop, step)
    return bins


def bin_width(bins):
    width = np.zeros(len(bins)-1)
    for i in range(len(bins)-1):
        width[i] = bins[i+1] - bins[i]
    return width


######

#read in clust results
t = np.genfromtxt('../data/clust_results709.txt')
ap_id = t[:,0].astype(int)
int_age = t[:,6]
int_age_low = t[:,7]
int_age_high = t[:,8]
int_mass = t[:,9]
int_mass_low = t[:,10]
int_mass_high = t[:,11]
int_av = t[:,12]
int_av_low = t[:,13]
int_av_high = t[:,14]
age = t[:,15]
age_low = t[:,16]
age_high = t[:,17]
mass = t[:,18]
mass_low = t[:,19]
mass_high = t[:,20]
av = t[:,21]
av_low = t[:,22]
av_high = t[:,23]




#alt ; is comment section

# n_max = np.zeros(len(h))
# n_min = np.zeros(len(h))

# for i in range(len(h)):
#     n_max[i] = np.max([h[i], h_low[i], h_high[i]])
#     n_min[i] = np.min([h[i], h_low[i], h_high[i]])

# width = bin_width(bins)

# n_mass = h / width
# n_mass_low = h_low / width
# n_mass_high = h_high / width

# n_max2 = n_max / width
# n_min2 = n_min / width

mass_bins=get_bins(2,5.6,0.1)
age_bins=get_bins(6,10.3,0.1)
start_bins_age = get_bins(6,10.2,0.1)


h,x_edges,y_edges=np.histogram2d(int_age,age,age_bins)
h_low,x_edges_low,y_edges_low=np.histogram2d(int_age_low,age_low,age_bins)
h_high,x_edges_high,y_edges_high=np.histogram2d(int_age_high,age_high,age_bins)

new_h = (h+h_low+h_high)/3.
#extent = [x_edges[0], x_edges[-1], y_edges[0], yedges[-1] ]


#create arrays for bin midpoints
x2=np.zeros(len(x_edges)-1)
for i in range(len(x_edges)-1):
    x2[i] = (x_edges[i] + x_edges[i+1])/2.
y2=x2


plt.close('all')
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')


plt.imshow(new_h,origin='lower',interpolation='none', cmap=plt.cm.Blues,vmin=-2)
plt.colorbar()

sdev=np.std(new_h)
#levels=np.arange(sdev,16,sdev)
levels=np.arange(1,16,3)
plt.contour(range(43), range(43), new_h, levels=levels, colors='black')

xx=np.arange(0,43,1)
yy=np.arange(0,43,1)
plt.plot(xx,yy,c='black')

ax.set_xticks([0, 10, 20, 30, 40])
ax.set_xticklabels(['6.0', '7.0', '8.0', '9.0', '10.0'])
ax.set_yticks([0, 10, 20, 30, 40])
ax.set_yticklabels(['6.0', '7.0', '8.0', '9.0', '10.0'])

#plt.xlim([6.5,10.0])
#plt.ylim([-9.5,-4.0])
plt.xlabel('MATCH log(age [yr])', fontsize=20)
plt.ylabel('Integrated log(age [yr])', fontsize=20)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)



#plt.show()
plt.savefig('plot_int_match_age.png')
