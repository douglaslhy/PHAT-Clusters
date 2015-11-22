import numpy as np
import sys
import matplotlib.pyplot as plt
import pyfits

'''creates table of all clusters' ap id, position, match age, mass, Av, and int age, mass, and Av'''


def get_results(dir_name, clust_id):

    '''finds fit results and returns them
    input:  dir_name to find results, eg '.', 'apold_profiles', 'SB_cuts'; clust_id
    returns:  bf_age, bf_av, bf_mass, P1_age, P2_age, P3_age, P1_mass, P2_mass, P3_mass, P1_av, P2_av, P3_av'''

    bf = np.genfromtxt(dir_name +'/ap' +str(clust_id)+ '/ap' +str(clust_id)+ '_bestfit.txt')
    bf_age = bf[0]
    bf_av = bf[1]
    bf_mass = bf[2]

    P_age = np.genfromtxt(dir_name +'/ap' +str(clust_id)+ '/ap' +str(clust_id)+ '_Log_Ageptiles.txt')
    P1_age = P_age[0]
    P2_age = P_age[1]
    P3_age = P_age[2]

    P_mass = np.genfromtxt(dir_name +'/ap' +str(clust_id)+ '/ap' +str(clust_id)+ '_Log_Mptiles.txt')
    P1_mass = P_mass[0]
    P2_mass = P_mass[1]
    P3_mass = P_mass[2]

    P_av = np.genfromtxt(dir_name +'/ap' +str(clust_id)+ '/ap' +str(clust_id)+ '_Avptiles.txt')
    P1_av = P_av[0]
    P2_av = P_av[1]
    P3_av = P_av[2]

    return bf_age, bf_av, bf_mass, P1_age, P2_age, P3_age, P1_mass, P2_mass, P3_mass, P1_av, P2_av, P3_av


def get_nstars(dir_name, clust_id):
    '''calculates number of stars in photometry file
    input: dir_name to find results, eg '.', 'apold_profiles', 'SB_cuts'; clust_id
    returns:  n_stars '''

    phot = np.genfromtxt(dir_name +'/ap' + str(clust_id) + '/ap' + str(clust_id) + '_phot.dat')
    nstars = len(phot)
    return nstars


done_clust = np.genfromtxt('done_list_415.txt', dtype=None)
sb = np.genfromtxt('SBcut_list.txt',dtype=None)
old = np.genfromtxt('apold_list.txt',dtype=None)

#integrated results from Morgan
int_data = pyfits.open('mf2015_obs_stats.fits')
tbdata = int_data[1].data
int_id_full = tbdata.field('ID_1')
int_age_full = tbdata.field('logage_E')
int_age_P1_full = tbdata.field('logage_p_16')
int_age_P3_full = tbdata.field('logage_p_84')
int_mass_full = tbdata.field('logmass_E')
int_mass_P1_full = tbdata.field('logmass_p_16')
int_mass_P3_full = tbdata.field('logmass_p_84')
int_av_full = tbdata.field('av_E')
int_av_P1_full = tbdata.field('av_p_16')
int_av_P3_full = tbdata.field('av_p_84')

#initialize arrays for match results
int_age=np.zeros(len(done_clust))
int_age_P1=np.zeros(len(done_clust))
int_age_P3=np.zeros(len(done_clust))
int_mass=np.zeros(len(done_clust))
int_mass_P1=np.zeros(len(done_clust))
int_mass_P3=np.zeros(len(done_clust))
int_av=np.zeros(len(done_clust))
int_av_P1=np.zeros(len(done_clust))
int_av_P3=np.zeros(len(done_clust))
bf_age=np.zeros(len(done_clust))
bf_mass=np.zeros(len(done_clust))
bf_av=np.zeros(len(done_clust))
P1_age=np.zeros(len(done_clust))
P1_mass=np.zeros(len(done_clust))
P1_av=np.zeros(len(done_clust))
P2_age=np.zeros(len(done_clust))
P2_mass=np.zeros(len(done_clust))
P2_av=np.zeros(len(done_clust))
P3_age=np.zeros(len(done_clust))
P3_mass=np.zeros(len(done_clust))
P3_av=np.zeros(len(done_clust))
ap_id=np.zeros(len(done_clust))
nstars=np.zeros(len(done_clust))

clust_name=[]

#find results by id and store in arrays
for i in range(len(done_clust)):
#for i in range(3):

    clust_id = int(done_clust[i].split('ap')[1])
    ap_id[i] = clust_id
    print ap_id[i]
    
    #find matching cluster by id in clust data table
    x = np.where(int_id_full == clust_id)

    if len(x[0]) > 0:
        int_age[i] = int_age_full[x]
        int_age_P1[i] = int_age_P1_full[x]
        int_age_P3[i] = int_age_P3_full[x]
        int_mass[i] = int_mass_full[x]
        int_mass_P1[i] = int_mass_P1_full[x]
        int_mass_P3[i] = int_mass_P3_full[x]
        int_av[i] = int_av_full[x]
        int_av_P1[i] = int_av_P1_full[x]
        int_av_P3[i] = int_av_P3_full[x]



    #read in match results
    try:

	#read in apold results if apid in apold list -- 'ap1' in sb

	if done_clust[i] in old:
	    bf_age[i], bf_av[i], bf_mass[i], P1_age[i], P2_age[i], P3_age[i], P1_mass[i], P2_mass[i], P3_mass[i], P1_av[i], P2_av[i], P3_av[i] = get_results('apold_profiles', clust_id)
            nstars[i] = get_nstars('apold_profiles', clust_id)

	elif done_clust[i] in sb:
	    bf_age[i], bf_av[i], bf_mass[i], P1_age[i], P2_age[i], P3_age[i], P1_mass[i], P2_mass[i], P3_mass[i], P1_av[i], P2_av[i], P3_av[i] = get_results('.', clust_id)
            nstars[i] = get_nstars('.', clust_id)

	else:
	    bf_age[i], bf_av[i], bf_mass[i], P1_age[i], P2_age[i], P3_age[i], P1_mass[i], P2_mass[i], P3_mass[i], P1_av[i], P2_av[i], P3_av[i] = get_results('.', clust_id)
            nstars[i] = get_nstars('.', clust_id)

    except IOError:
        print str(int_id)+ '_bestfit.txt does not exist'


#save results to file
dt=np.dtype([('clust_id', 'i'), ('int_age', 'd'), ('int_age_P1', 'd'), ('int_age_P3', 'd'), ('int_mass', 'd'), ('int_mass_P1', 'd'), ('int_mass_P3', 'd'), ('int_av', 'd'), ('int_av_P1', 'd'), ('int_av_P3', 'd'), ('bf_age', 'd'), ('P1_age', 'd'), ('P2_age', 'd'), ('P3_age', 'd'), ('bf_mass', 'd'), ('P1_mass', 'd'), ('P2_mass', 'd'), ('P3_mass', 'd'), ('bf_av', 'd'), ('P1_av', 'd'), ('P2_av', 'd'), ('P3_av', 'd'), ('nstars', 'd')])

a=np.zeros(len(done_clust), dt)

a['clust_id'] = ap_id
a['int_age'] = int_age
a['int_age_P1'] = int_age_P1
a['int_age_P3'] = int_age_P3
a['int_mass'] = int_mass
a['int_mass_P1'] = int_mass_P1
a['int_mass_P3'] = int_mass_P3
a['int_av'] = int_av
a['int_av_P1'] = int_av_P1
a['int_av_P3'] = int_av_P3
a['bf_age'] = bf_age
a['bf_mass'] = bf_mass
a['bf_av'] = bf_av
a['P1_age'] = P1_age
a['P1_mass'] = P1_mass
a['P1_av'] = P1_av
a['P2_age'] = P2_age
a['P2_mass'] = P2_mass
a['P2_av'] = P2_av
a['P3_age'] = P3_age
a['P3_mass'] = P3_mass
a['P3_av'] = P3_av
a['nstars'] = nstars

np.savetxt('int_match_results_noSB.txt', a, '%15s')



