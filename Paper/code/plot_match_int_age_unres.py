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

'''plot int age vs match age, with points sized by diff in unresolved flux

Reads in:
../data/unres_results512.txt 
../data/clust_results709.txt

Output:  int_match_age_unres.png
'''


#get unres results
unres=np.loadtxt('../data/unres_results512.txt')
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


#read in clust results
t = np.genfromtxt('../data/clust_results709.txt')
clust_id = t[:,0].astype(int)
int_age = t[:,6]
int_age_low = t[:,7]
int_age_high = t[:,8]
int_mass = t[:,9]
int_mass_low = t[:,10]
int_mass_high = t[:,11]
int_av = t[:,12]
int_av_low = t[:,13]
int_av_high = t[:,14]
bf_age = t[:,15]
age_low = t[:,16]
age_high = t[:,17]
mass = t[:,18]
mass_low = t[:,19]
mass_high = t[:,20]
av = t[:,21]
av_low = t[:,22]
av_high = t[:,23]



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


good = np.where(int_mass > 0)


xdata = bf_age
ydata = color_diff_48

xbin_min = 6.2
xbin_max = 10.4
ybin_min = -0.1
ybin_max = 5.0
xbin_width = 0.3
ybin_width = 0.1

xlabel = 'MATCH log(age [yr])'
ylabel = 'F475W - F814W Color Difference'


plt.scatter(bf_age+np.random.normal(0,0.02,len(bf_age)), int_age, s=mag_diff_475*20, edgecolors='none', alpha=0.7)
plt.ylim([6.2, 10.4])
plt.xlim([6.2, 10.4])
plt.xlabel('MATCH log(age[yr])', fontsize=25)
plt.ylabel('Integrated log(age[yr])', fontsize=25)


#plt.show()
plt.savefig('int_match_age_unres.png')





