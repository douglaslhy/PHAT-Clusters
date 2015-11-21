These files are originally from the chex directory /media/Raid5/lori/Match/Dec2014_runs/SB_cuts

First, find_rad.py was run to find the appropriate pix radius corresponding to a SB of 18 mag/arcsec^2 in F814W.  This value was chosen from 
looking at Nelson's choices for cutoff radius for his old clusters.
SBpixrad.txt was manually created from the output.

This pixel radius was used in cut_center_match.py to create new phot and fake files with only stars outside of this radius to run through MATCH.
