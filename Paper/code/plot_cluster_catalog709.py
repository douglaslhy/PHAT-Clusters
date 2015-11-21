import numpy as np
import pyfits
import pyregion
import matplotlib.pyplot as plt
from astLib.astWCS import *
from matplotlib import rc
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
plt.ion()

'''plots NUV and 24um images of M31 along with cluster sample, 
color-coded by age, mass, and Av

Reads data from:
apdata_cluster_v2.fits    (in ../data directory)
clust_results709.txt      (in ../data directory)
m31_mips24.fits           (NOT on github)
m31_nuv.fits              (NOT on github)

Outputs:
plot_cluster_catalog2.png, cluster_catalog_age.png, cluster_catalog_mass.png, cluster_catalog_av.png
'''


def read_img(filename):
    '''read in image data from fits file'''
    f = pyfits.open(filename)
    image = f[0].data
    header = f[0].header
    wcs = WCS(header, mode='pyfits')
    return image, wcs
    

def read_table(filename):
    '''read in table data from fits file'''
    f = pyfits.open(filename)
    table = f[1].data
    return table


def deg_to_pix(coord):
    '''converts coordinate from degrees to pixels'''
    arc = coord * 3600.		#convert from degrees to arcsec
    pix = arc / 0.05		#convert from arcsec to pixels
    return pix


def read_reg(region):
    '''read in ra, dec, and radius data from region file'''
    
    reg=pyregion.open(region)
    n=len(reg)
    racoord=np.empty(n)
    deccoord=np.empty(n)
    r1=np.empty(n)
    r2=np.empty(n)
    angle=np.empty(n)
    for i in range(n):
        racoord[i]=(reg[i].coord_list[0])
        deccoord[i]=(reg[i].coord_list[1])
        r1[i]=(reg[i].coord_list[2])
	r2[i]=(reg[i].coord_list[3])
	angle[i]=(reg[i].coord_list[4])
        
    return racoord, deccoord, r1, r2, angle


#Read in images
image, wcs = read_img('m31_mips24.fits')
img,wcs1=read_img('m31_nuv.fits')


#Set data path
data_path = '../data/'


#clust data w/ ra, dec
clust_file = pyfits.open(data_path + 'apdata_cluster_v2.fits')
clust = clust_file[1].data
clust_ra = clust.field('RA')
clust_dec = clust.field('DEC')

#read in clust results - these are in the same order as from the apdata file
t = np.genfromtxt('clust_results709.txt')
ap_id = t[:,0].astype(int)
R_ap = t[:,1]
n_stars = t[:,2]
n_bg = t[:,3]
cat = np.genfromtxt(data_path + 'clust_results709.txt', usecols=(4), dtype=None)
best = np.genfromtxt(data_path + 'clust_results709.txt', usecols=(5), dtype=None)
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

#masks for best method
M = np.where(best == 'M')
I = np.where(best == 'I')
O = np.where(best == 'O')


#cluster locations in pix
clust_x = np.zeros(len(ap_id))
clust_y = np.zeros(len(ap_id))
for i in range(len(ap_id)):
    coord = WCS.wcs2pix(wcs1, clust_ra[i], clust_dec[i])
    clust_x[i] = coord[0]
    clust_y[i] = coord[1]

rot = WCS.getRotationDeg(wcs)


min_max = WCS.getImageMinMaxWCSCoords(wcs1)
coord_min = WCS.wcs2pix(wcs1, 11.75, 41.25)
coord_med = WCS.wcs2pix(wcs1, 11.0, 41.75)
coord_max = WCS.wcs2pix(wcs1, 10.25, 42.25)
coord2 = WCS.wcs2pix(wcs1, 9.5, 40.75)
coord3 = WCS.wcs2pix(wcs1, 8.75, 40.25)


#phat footprint
poly=(11.950405,42.194671,11.9056,42.200102,11.861626,42.205351,11.836935,42.208257,11.844819,42.243896,11.852662,42.280627,11.861011,42.317612,11.816533,42.323005,11.772255,42.328223,11.727688,42.333478,11.683322,42.338506,11.638951,42.343777,11.594082,42.349192,11.585808,42.311927,11.578047,42.275489,11.570231,42.239467,11.553165,42.241316,11.508989,42.246432,11.464946,42.251582,11.420197,42.257136,11.412109,42.220167,11.404102,42.183353,11.396181,42.14711,11.353673,42.152416,11.345636,42.115388,11.337822,42.078774,11.329952,42.04268,11.31188,42.044603,11.26803,42.049937,11.223261,42.055311,11.21517,42.01822,11.207371,41.981492,11.199445,41.945387,11.156283,41.950624,11.111362,41.955955,11.103418,41.918907,11.095623,41.88196,11.087774,41.846159,11.052493,41.850496,11.044456,41.813442,11.036638,41.776698,11.028904,41.740735,11.023533,41.741428,11.016237,41.707129,10.986565,41.717002,10.966182,41.682993,10.94596,41.649033,10.92576,41.615328,10.892984,41.62633,10.872748,41.592147,10.852463,41.558066,10.83248,41.524656,10.800777,41.535318,10.780408,41.500963,10.76034,41.467229,10.740293,41.433543,10.708758,41.444023,10.688466,41.409807,10.668335,41.375848,10.648432,41.342355,10.636775,41.346465,10.616344,41.31208,10.596349,41.278161,10.576099,41.243925,10.616911,41.230206,10.657015,41.216836,10.697586,41.203296,10.73814,41.189899,10.778194,41.176435,10.818991,41.162904,10.819407,41.163425,10.857908,41.150508,10.897981,41.137003,10.938523,41.123692,10.97877,41.11021,11.018794,41.096768,11.059768,41.0831,11.079806,41.117303,11.100071,41.15119,11.120008,41.184762,11.132299,41.180704,11.152751,41.214894,11.172736,41.248717,11.192951,41.282535,11.224673,41.271766,11.245271,41.305886,11.265333,41.339849,11.285482,41.373339,11.317095,41.362597,11.337701,41.396856,11.357839,41.430802,11.378064,41.464379,11.411506,41.45308,11.432052,41.487427,11.452195,41.521147,11.472426,41.554603,11.481575,41.551519,11.502046,41.585593,11.522471,41.619768,11.543125,41.653832,11.502245,41.667619,11.461762,41.681023,11.428491,41.692284,11.466793,41.687656,11.510731,41.682505,11.554869,41.67718,11.563112,41.714302,11.570797,41.750699,11.578839,41.787092,11.614883,41.782587,11.623102,41.8196,11.631046,41.856249,11.639203,41.892532,11.682915,41.887317,11.727454,41.881819,11.73561,41.918823,11.743493,41.955671,11.751582,41.991684,11.770534,41.989632,11.814509,41.984418,11.858892,41.97887,11.867344,42.015755,11.875173,42.052594,11.883361,42.089376,11.925919,42.08403,11.934072,42.120911,11.942022,42.157742)

poly_arr=np.array(poly)
poly_coords=poly_arr.reshape((127,2))

phat_coord=np.zeros([127,2])
for i in range(127):
    phat_coord[i][0], phat_coord[i][1] = WCS.wcs2pix(wcs1, poly_coords[i][0], poly_coords[i][1])

plt.plot(phat_coord[:,0],phat_coord[:,1])



####plot cluster positions w/ NUV image as background

plt.close('all')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.imshow(img,origin='lower',cmap=plt.get_cmap('Greys'),vmin=0.001,vmax=0.02)

plt.scatter(clust_x, clust_y,  s=20, alpha=0.8, edgecolors='none', c='blue') 
plt.plot(phat_coord[:,0],phat_coord[:,1], color='black', lw=2)

plt.xticks((coord_min[0], coord_med[0], coord_max[0], coord2[0], coord3[0]), (11.75, 11.0, 10.25, 9.5, 8.75))
plt.yticks((coord_min[1], coord_med[1], coord_max[1], coord2[1], coord3[1]), (41.25, 41.75, 42.25, 40.75, 40.25))
plt.xlim([100, 5000])
plt.ylim([600, 6200])
plt.xlabel('RA (degrees)', fontsize=30)
plt.ylabel('Dec (degrees)', fontsize=30)
plt.yticks(fontsize=22)
plt.xticks(fontsize=22)


plt.savefig('plot_cluster_catalog2.png')



####plot cluster positions color-coded by age w/ NUV image as background

plt.close('all')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.imshow(img,origin='lower',cmap=plt.get_cmap('Greys'),vmin=0.001,vmax=0.02)

plt.scatter(clust_x[O], clust_y[O], s=20, edgecolors='none', alpha=0.8, c=C11_age[O], cmap=plt.get_cmap('seismic'), vmin=6.0, vmax=10.0)
plt.scatter(clust_x[I], clust_y[I], s=20, edgecolors='none',  alpha=0.8, c=int_age[I], cmap=plt.get_cmap('seismic'), vmin=6.0, vmax=10.0)
plt.scatter(clust_x[M], clust_y[M], s=20, edgecolors='none',  alpha=0.8, c=cmd_age[M], cmap=plt.get_cmap('seismic'), vmin=6.0, vmax=10.0)

cbar = plt.colorbar(pad=0)
cbar.set_label('$\log_{10}$($t$ [yr])', fontsize=25)
cbar.ax.tick_params(labelsize=20)
plt.xticks((coord_min[0], coord_med[0], coord_max[0]), (11.75, 11.0, 10.25))
plt.yticks((coord_min[1], coord_med[1], coord_max[1]), (41.25, 41.75, 42.25))
plt.xlim([300, 3600])
plt.ylim([2550, 6000])
plt.xlabel('RA (degrees)', fontsize=30)
plt.ylabel('Dec (degrees)', fontsize=30)
plt.yticks(fontsize=22)
plt.xticks(fontsize=22)

plt.savefig('cluster_catalog_age.png')


####plot cluster positions color-coded by mass w/ NUV image as background

plt.close('all')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.imshow(img,origin='lower',cmap=plt.get_cmap('Greys'),vmin=0.001,vmax=0.02)


plt.scatter(clust_x[O], clust_y[O], s=20, alpha=0.8, edgecolors='none', c=C11_mass[O], cmap=plt.get_cmap('jet'), vmin=2.5, vmax=4.6)
plt.scatter(clust_x[I], clust_y[I], s=20, alpha=0.8, edgecolors='none', c=int_mass[I], cmap=plt.get_cmap('jet'), vmin=2.5, vmax=4.6)
plt.scatter(clust_x[M], clust_y[M], s=20, alpha=0.8, edgecolors='none', c=cmd_mass[M], cmap=plt.get_cmap('jet'), vmin=2.5, vmax=4.6)

cbar = plt.colorbar(pad=0)
cbar.set_label('$\log_{10}$($M$ [$M_{\odot}$])', fontsize=25)
cbar.ax.tick_params(labelsize=20)
plt.xticks((coord_min[0], coord_med[0], coord_max[0]), (11.75, 11.0, 10.25))
plt.yticks((coord_min[1], coord_med[1], coord_max[1]), (41.25, 41.75, 42.25))
plt.xlim([300, 3600])
plt.ylim([2550, 6000])
plt.xlabel('RA (degrees)', fontsize=30)
plt.ylabel('Dec (degrees)', fontsize=30)
plt.yticks(fontsize=22)
plt.xticks(fontsize=22)

plt.savefig('cluster_catalog_mass.png')


####plot cluster positions color-coded by mass w/ 24um image as background

plt.close('all')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.imshow(image,origin='lower',cmap=plt.get_cmap('Greys'),vmin=0.001,vmax=0.02)


plt.scatter(clust_x[O], clust_y[O], s=20, alpha=0.8, edgecolors='none', c=C11_av[O], cmap=plt.get_cmap('Reds'), vmin=2.0, vmax=3.0)
plt.scatter(clust_x[I], clust_y[I], s=20, alpha=0.8, edgecolors='none', c=int_av[I], cmap=plt.get_cmap('Reds'), vmin=2.0, vmax=3.0)
plt.scatter(clust_x[M], clust_y[M], s=20, alpha=0.8, edgecolors='none', c=cmd_av[M], cmap=plt.get_cmap('Reds'), vmin=2.0, vmax=3.0)

cbar = plt.colorbar(pad=0)
cbar.set_label('$A_{V}$ [mag]', fontsize=20)
cbar.ax.tick_params(labelsize=18)
plt.xticks((coord_min[0], coord_med[0], coord_max[0]), (11.75, 11.0, 10.25))
plt.yticks((coord_min[1], coord_med[1], coord_max[1]), (41.25, 41.75, 42.25))
plt.xlim([300, 3600])
plt.ylim([2550, 6000])
plt.xlabel('RA (degrees)', fontsize=30)
plt.ylabel('Dec (degrees)', fontsize=30)
plt.yticks(fontsize=22)
plt.xticks(fontsize=22)

plt.savefig('cluster_catalog_av.png')




