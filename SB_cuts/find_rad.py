import numpy as np

'''Reads SB profiles from apclust_profiles/apXXX_profile.dat
 and reports max radius less than given SB

Inputs:  
done_list_1215.txt:  list of names of already run clusters
apclst_profiles/apXXX_profile.dat:  ap cluster profiles where 
    columns are radius (pix), sb475 (mag/arcsec^2), sb814  (mag/arcsec^2)

Output:  results were recorded manually in SBpixrad.txt

To run in Python:
%run find_rad.py
'''


#read in list of done clusters
clust = np.genfromtxt('done_list_1215.txt', dtype=None)

sb475_lim = 20.
sb814_lim = 18.

#loop over all profiles
rad475 = np.zeros(len(clust))
rad814 = np.zeros(len(clust))
clust_id = np.zeros(len(clust)).astype(int)

for i in range(len(clust)):


    try:

        #read in cluster profile
        prof = np.genfromtxt('apclst_profiles/' +clust[i]+ '_profile.dat')
        rad = prof[:,0]
        sb475 = prof[:,1]
        sb814 = prof[:,2]
        lim475 = np.where(sb475 < sb475_lim)
        lim814 = np.where(sb814 < sb814_lim)

        if len(lim475[0] > 0):
            rad475[i] = np.max(rad[lim475])
            print rad475[i]
        if len(lim814[0] > 0):
            rad814[i] = np.max(rad[lim814])
	    clust_id[i] = int(clust[i].split('ap')[1])
            print rad814[i]


    except IOError:
        print clust[i]+ '_profile.dat does not exist'

