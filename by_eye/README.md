#Bye Eye Check of Best Fit

The cluster results table, cluster_results709.txt contains the results of both the by-eye check and for 100+ older 
clusters in Nelson Caldwell's 2011 table.

The N_stars column is just the number of stars in the photometry file, and for the N_bg column I divided the number 
of stars in the background file by 10 to scale it to the aperture size, giving a predicted value for number of background stars.

The Best_method column represents which method I think should be used based on my by-eye check (by_eye_results526.txt). 

M:  Match fit should be used 

I:  Integrated fit should be used

O:  the values in the C11 columns (C11_Age_Gyr, C11_logMass, C11_E(B-V), C11_Z, and C11_SigmaZ) should be used.  
These are from Nelson's 2011 paper.

N:  neither Match nor integrated fits looked acceptable (only 14 clusters)

X:  no Match or integrated results available (only 21 clusters, most of these were the ones Cliff manually added)
