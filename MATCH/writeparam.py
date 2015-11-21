#based off code from Dan Weisz
#run_match.py calls this

import numpy as np
import sys

def main(clustname):

    #i = np.int(sys.argv[1])

    #generate match parameter files

    # line 1
    imf = -1.
    #imf = np.arange(0.01, 5.02, 0.1)
    dmod0 = 24.47
    dmod1 = 24.47
    dmod_delta = 0.1
    av0 = 0.00
    av = 3.00
    av_delta = 0.05
    line1 = ('%4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f') % (imf, dmod0, dmod1, dmod_delta, av0, av, av_delta)


    #line 2
    zmin = -0.2
    zmax = 0.1
    z_delta = 0.1
    line2 = ('%4.2f %4.2f %4.2f') % (zmin, zmax, z_delta)


    # line 3

    bf = 0.35
    baddetect = 0.000001
    line3 = ('%4.2f %4.2f %4.2f') % (bf, baddetect, baddetect)

    # line 4
    ncmd = 1
    line4 = ('%1i') % (ncmd)


    # line 5

    mbin = 0.1
    cbin = 0.1
    smooth = 3
    cmin = -0.5
    cmax = 3.5
    filter1 = 'WFC475W'
    filter2 = 'WFC814W'

    line5 = ('%4.2f %4.2f %1i %4.2f %4.2f') % (mbin, cbin, smooth, cmin, cmax)
    line5 = line5 + ' ' + filter1 + ',' + filter2

    # line 6 & 7

    f1_brightlim = 15.
    f1_faintlim = 27.
    f2_brightlim = 14.5
    f2_faintlim = 26.5
    line6 = ('%4.2f %4.2f') % (f1_brightlim, f1_faintlim)
    line6 = line6 + ' ' + filter1
    line7 = ('%4.2f %4.2f') % (f2_brightlim, f2_faintlim)
    line7 = line7 + ' ' + filter2

    # line 8

    include = 0
    exclude = 0

    line8 = ('%1i %1i') % (include, exclude)

    # time bins

    timebins0 = np.arange(6.60, 10.14, 0.1)
    timebins1 = np.append(np.arange(6.70, 10.14, 0.1), 10.15)

    # line 9 

    tbins = len(timebins0)
    line9 = ('%2i') % (tbins)

    print line1
    print line2
    print line3
    print line4
    print line5
    print line6
    print line7
    print line8
    print line9


    for j in range(len(timebins0)):
        print ('%4.2f %4.2f') % (timebins0[j], timebins1[j])



    # background

    background1 = '-1 3 -1.{0}/bg.dat'.format(clustname)
    background2 = '-1 1 -1'

    print background1
    print background2

if __name__=='__main__':
    sys.exit(main(sys.argv[1], sys.argv[2]))
