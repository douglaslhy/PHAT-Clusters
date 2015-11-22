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

------------------------------
The Flag column is the classification I gave each cluster in the by-eye check.  These values mean:

B = Both integrated and Match fits look acceptable

MAY = Match fit looks acceptable, integrated fit not acceptable, cluster is young (log(Match Age) <= 8.5)

MAO = Match fit looks acceptable, integrated fit not acceptable, cluster is older (log(Match Age) > 8.5)

IAY = Integrated fit looks acceptable, Match fit not acceptable, cluster is young

IAO = Integrated fit looks acceptable, Match fit not acceptable, cluster is older

NS = I'm not sure what is acceptable, either because the CMD does not have many stars or the background is very crowded, or something else weird is going on

NSMB = I'm not sure if either fit is doing well, but the Match fit looks like it's doing better than the integrated fit, a lot of these are the 1 Gyr Match fits where the integrated fits don't look good

NSIB = I'm not sure, but the integrated fit looks better than the Match fit

N = Neither fit looks acceptable, not many of these

MAMO = Only Match results available, the Match fit looks acceptable

MUMO = Only Match results available, the Match fit does not look acceptable

XXX = No Match results available, most of these are the clusters at the very end of the catalog, I still need to check if any of these have integrated results

From these values, I chose the best method and put that in the best method column.  When Flag = B or NS, I used the 'default' method, which is Match for ages <= 8.5 and int for older.
