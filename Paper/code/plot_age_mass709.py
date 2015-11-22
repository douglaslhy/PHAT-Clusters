import numpy as np
import pylab as plt
import os, sys, shutil
import pyfits
from matplotlib import cm
from matplotlib import rc
from scipy import interpolate
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)


'''
plot Age vs mass for final cluster results, color-coded by Av, excludes clusters labelled 'not sure'
Reads:
clust_results709.txt

Output:
age_mass_comp1018.png
'''


#read in clust results
t = np.genfromtxt('clust_results709.txt')
ap_id = t[:,0].astype(int)
R_ap = t[:,1]
n_stars = t[:,2]
n_bg = t[:,3]
cat = np.genfromtxt('clust_results709.txt', usecols=(4), dtype=None)
best = np.genfromtxt('clust_results709.txt', usecols=(5), dtype=None)
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



#age-mass plot, color is which method used

plt.close('all')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')


plt.scatter(cmd_age[M] + np.random.normal(0, 0.02, len(cmd_age[M])), cmd_mass[M], color='blue', edgecolors='none', alpha=0.5)
plt.scatter(int_age[I], int_mass[I], color='red', edgecolors='none', alpha=0.5)
plt.xlim([6.2, 10.4])
plt.ylim([1.6, 5.8])
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlabel('$\log_{10}$($t$ [yr])', fontsize=20)
plt.ylabel('$\log_{10}$($M$ [$M_{\odot}$])', fontsize=20)






#age-mass plot, color is log(sigma,match / sigma,int)

sigma_age_diff = np.log10(cmd_age_width / int_age_width)
sigma_mass_diff = np.log10(cmd_mass_width / int_mass_width)
sigma_av_diff = np.log10(cmd_av_width / int_av_width)

plt.close('all')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')


plt.scatter(cmd_age[M] + np.random.normal(0, 0.02, len(cmd_age[M])), cmd_mass[M], c=sigma_av_diff[M], edgecolors='none', alpha=0.5, cmap=cm.seismic, vmin=-1.6, vmax=1.6)
plt.scatter(int_age[I], int_mass[I], c=sigma_av_diff[I], edgecolors='none', alpha=0.5, cmap=cm.seismic, vmin=-1.6, vmax=1.6)
plt.xlim([6.2, 10.4])
plt.ylim([1.8, 5.8])
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlabel('$\log_{10}$($t$ [yr])', fontsize=20)
plt.ylabel('$\log_{10}$($M$ [$M_{\odot}$])', fontsize=20)
cbar = plt.colorbar(pad=0)
cbar.set_label('$\log_{10}$ ($\sigma_{A_{V},MATCH}$ / $\sigma_{A_{V},Integrated}$)', fontsize=18)


#age-mass plot, point sizes are fraction of bg stars
bg_frac = n_bg / n_stars

plt.close('all')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')


plt.scatter(cmd_age[M] + np.random.normal(0, 0.02, len(cmd_age[M])), cmd_mass[M], s=bg_frac[M]*20, edgecolors='none', alpha=0.5)
plt.scatter(int_age[I], int_mass[I], s=bg_frac[I]*20, edgecolors='none', alpha=0.5)
plt.xlim([6.2, 10.4])
plt.ylim([1.8, 5.8])
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlabel('$\log_{10}$($t$ [yr])', fontsize=20)
plt.ylabel('$\log_{10}$($M$ [$M_{\odot}$])', fontsize=20)
#cbar = plt.colorbar(pad=0)
#cbar.set_label('$\log_{10}$ ($\sigma_{A_{V},MATCH}$ / $\sigma_{A_{V},Integrated}$)', fontsize=18)



#match-int age plot, point sizes are fraction of bg stars
bg_frac = n_bg / n_stars

plt.close('all')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

line1 = [5,6,7,8,9,10,11]
line2 = [5.5,6.5,7.5,8.5,9.5,10.5,11.5]

plt.scatter(cmd_age + np.random.normal(0, 0.02, len(cmd_age)), int_age, s=bg_frac*40, edgecolors='none', alpha=0.5)
plt.plot(line1, line1, color='black')
plt.plot(line1, line2, color='black', ls='dashed')
plt.plot(line2, line1, color='black', ls='dashed')
plt.xlim([6.2, 10.4])
plt.ylim([6.2, 10.4])
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlabel('MATCH $\log_{10}$($t$ [yr])', fontsize=20)
plt.ylabel('Integrated $\log_{10}$($t$ [yr])', fontsize=20)
#cbar = plt.colorbar(pad=0)
#cbar.set_label('$\log_{10}$ ($\sigma_{A_{V},MATCH}$ / $\sigma_{A_{V},Integrated}$)', fontsize=18)


#plt.show()



#match-int age, cc mass

plt.close('all')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

line1 = [5,6,7,8,9,10,11]
line2 = [5.5,6.5,7.5,8.5,9.5,10.5,11.5]

plt.scatter(cmd_age + np.random.normal(0, 0.02, len(cmd_age)), int_age, c=cmd_mass, edgecolors='none', alpha=0.5, vmin=2.2)
plt.plot(line1, line1, color='black')
plt.plot(line1, line2, color='black', ls='dashed')
plt.plot(line2, line1, color='black', ls='dashed')
plt.xlim([6.2, 10.4])
plt.ylim([6.2, 10.4])
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlabel('MATCH $\log_{10}$($t$ [yr])', fontsize=20)
plt.ylabel('Integrated $\log_{10}$($t$ [yr])', fontsize=20)
cbar = plt.colorbar(pad=0)
#cbar.set_label('$\log_{10}$ ($\sigma_{A_{V},MATCH}$ / $\sigma_{A_{V},Integrated}$)', fontsize=18)


#plt.show()




#age-mass plot, color is Av
plt.close('all')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')



plt.scatter(int_age[I], int_mass[I], c=int_av[I], edgecolors='none', cmap=cm.rainbow, alpha=0.5)
plt.scatter(C11_age[O], C11_mass[O], c=C11_av[O], edgecolors='none', cmap=cm.rainbow, alpha=0.5)
plt.scatter(cmd_age[M] + np.random.normal(0, 0.02, len(cmd_age[M])), cmd_mass[M], c=cmd_av[M], edgecolors='none', cmap=cm.rainbow, alpha=0.5)

plt.xlim([6.2, 10.4])
plt.ylim([1.8, 7.0])
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlabel('$\log_{10}$($t$ [yr])', fontsize=20)
plt.ylabel('$\log_{10}$($M$ [$M_{\odot}$])', fontsize=20)
cbar = plt.colorbar(pad=0)
cbar.set_label('$A_{V}$ [mag]', fontsize=18)

#plt.show()


#don't plot NS clusters

M_sure = np.where((cat != 'NS') & (cat != 'N') & (best == 'M'))
I_sure = np.where((cat != 'NS') & (cat != 'N') & (best == 'I'))


plt.close('all')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')



plt.scatter(int_age[I_sure], int_mass[I_sure], c=int_av[I_sure], edgecolors='none', cmap=cm.rainbow, alpha=0.5)
plt.scatter(C11_age[O], C11_mass[O], c=C11_av[O], edgecolors='none', cmap=cm.rainbow, alpha=0.5)
plt.scatter(cmd_age[M_sure] + np.random.normal(0, 0.02, len(cmd_age[M_sure])), cmd_mass[M_sure], c=cmd_av[M_sure], edgecolors='none', cmap=cm.rainbow, alpha=0.5)

#completeness lines
x1=[6.5, 7.0]
y1=[2.65, 2.65]
x2=[7.0, 7.5]
y2=[2.7, 2.7]
x3=[7.5, 8.0]
y3=[2.77, 2.77]
x4=[8.0, 8.5]
y4=[2.87, 2.87]
x5=[8.5, 9.0]
y5=[3.15, 3.15]
x6=[9.0, 9.5]
y6=[3.7, 3.7]
x7=[9.5, 10.0]
y7=[4.2, 4.2]
#plt.plot(x1, y1, ls='dashed', color='black', lw=2)
#plt.plot(x2, y2, ls='dashed', color='black', lw=2)
#plt.plot(x3, y3, ls='dashed', color='black', lw=2)
#plt.plot(x4, y4, ls='dashed', color='black', lw=2)
#plt.plot(x5, y5, ls='dashed', color='black', lw=2)
#plt.plot(x6, y6, ls='dashed', color='black', lw=2)
#plt.plot(x7, y7, ls='dashed', color='black', lw=2)

#interpolate
med_x = np.array([6.75, 7.25, 7.75, 8.25, 8.75, 9.25, 9.75])
med_y = np.array([2.65, 2.7, 2.77, 2.87, 3.15, 3.7, 4.2])
f = interpolate.interp1d(med_x, med_y)
xnew = np.arange(6.75, 10.25, 0.5)
ynew = f(xnew)

plt.plot(xnew, ynew, lw=2, color='black')

plt.xlim([6.2, 10.4])
plt.ylim([1.8, 7.0])
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlabel('$\log_{10}$($t$ [yr])', fontsize=20)
plt.ylabel('$\log_{10}$($M$ [$M_{\odot}$])', fontsize=20)
cbar = plt.colorbar(pad=0)
cbar.set_label('$A_{V}$ [mag]', fontsize=18)

#plt.show()
plt.savefig('age_mass_comp1018.png')
