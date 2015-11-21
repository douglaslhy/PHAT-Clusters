import numpy as np
import pyfits
import  matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter


'''Create new photometry files for stars and ASTs that are outside the SB cut radius
Files needed:  
SBpixrad.txt: contains clust_id and radius in pixels at which to cut
apdata-cluster_6phot_v4.fits: contains ra, dec, ap_rad for each cluster
reads phot file from ../apXXX/apXXX_phot.fits
reads fake file from ../apXXX/apXXX.dst.fake.fits

Outputs:
apXXXSBcut.png:  plot showing which stars were removed
apXXXSBcut_fake.png:  plot showing which ASTs were removed
apXXX_phot.dat:  new photometry file containing only stars outside the SB cut radius
apXXX_fake.dat:  new AST file containing only stars outside the SB cut radius
SBcut_list.txt:  list of clusters w/ SB cut applied

To run in Python:
%run cut_center_match.py
'''


def get_rad(SBradfile):
    '''read in radius at which to cut - units are pix'''
    SBrad = np.genfromtxt(SBradfile)
    clust_id = SBrad[:,0].astype(int)
    rad = SBrad[:,1]

    return clust_id, rad


def get_phot(phot_file):
    '''read in photometry file for cluster - ra units are degrees'''
    phot = pyfits.open(phot_file)
    phot_data = phot[1].data
    ra = phot_data.field('RA')
    dec = phot_data.field('DEC')
    phot_f475 = phot_data.field('F475W_VEGA')
    phot_f814 = phot_data.field('F814W_VEGA')
    return ra, dec, phot_f475, phot_f814

def get_fake(fake_file):
    '''read in fake file for cluster - ra units are degrees'''
    fake = pyfits.open(fake_file)
    fake_data = fake[1].data
    fake_ra = fake_data.field('RA_J2000')
    fake_dec = fake_data.field('DEC_J2000')
    phot_f475_in = fake_data.field('F475W_IN')
    phot_f814_in = fake_data.field('F814W_IN')
    phot_f475_out = fake_data.field('F475W_VEGA')
    phot_f814_out = fake_data.field('F814W_VEGA')
    return fake_ra, fake_dec, phot_f475_in, phot_f814_in, phot_f475_out, phot_f814_out

def find_cen(clust_file):
    '''find cluster center'''
    clust_data = pyfits.open(clust_file)
    clust_table = clust_data[1].data
    ap_id = clust_table.field('ID')
    ap_ra = clust_table.field('RA')
    ap_dec = clust_table.field('DEC')
    ap_rad = clust_table.field('APRAD')		#in arcsec?
    return ap_id, ap_ra, ap_dec, ap_rad

def dist_from_cen(ra, dec, ra_cen, dec_cen):
    '''find each star's dist from center - units are degrees'''
    radiff = np.abs(ra - ra_cen)
    decdiff = np.abs(dec - dec_cen)
    dist = np.sqrt((radiff * np.cos((dec) / 180.*np.pi))**2. + decdiff**2.)		#separation in degrees
    return dist

def pix_to_deg(pix):
    '''convert pixels to degrees'''
    deg = pix * 0.05 / 3600.
    return deg

def select_stars(rad, dist):
    '''select stars within two annuli - returns place array'''
    in_stars = np.where(dist >= rad)
    return in_stars

def print_stars(phot_f475, phot_f814, out_file):
    '''print only stars within annuli to new phot file'''
    phot = np.zeros([len(phot_f475), 2])
    phot[:,0] = phot_f475
    phot[:,1] = phot_f814
    np.savetxt(out_file, phot, fmt='%f')

def rec_diff(out_mag, in_mag):
    '''tests for recovery for asts
    outputs out-in, if recovered, 99.999 if not recovered'''

    out_in = out_mag - in_mag
    rec = np.abs(out_in) <= 0.75		#test for recovery
    diff = np.zeros(len(out_mag))
    diff[rec] = out_in[rec]			#recovered stars = out - in
    diff[~rec] = 99.999			#unrecovered stars = 99.999
    return diff

def print_fake_stars(phot_f475_in, phot_f814_in, phot_f475_out, phot_f814_out, out_file):
    '''print only stars within annuli to new fake file'''
    diff_f475 = rec_diff(phot_f475_out, phot_f475_in)
    diff_f814 = rec_diff(phot_f814_out, phot_f814_in)
    fake = np.zeros([len(phot_f475_out), 4])
    fake[:,0] = phot_f475_out
    fake[:,1] = phot_f814_out
    fake[:,2] = diff_f475			#out - in
    fake[:,3] = diff_f814			#out - in
    np.savetxt(out_file, fake, fmt='%f')

def plot_rad(ap_in, ra_cen, dec_cen, ra, dec, in_stars, figname):
    '''plot cluster stars along with two annuli -- degrees'''

    plt.close('all')
    plt.axes()
    circle1=plt.Circle((ra_cen, dec_cen), radius=ap_in, fill=False)
    plt.gca().add_patch(circle1)
    plt.axis('scaled')
    plt.scatter(ra, dec)
    plt.scatter(ra[in_stars], dec[in_stars], color='g')
    plt.xlim([np.min(ra),np.max(ra)])
    plt.ylim([np.min(dec), np.max(dec)])
    plt.gca().xaxis.set_major_formatter(FormatStrFormatter('%1.4f'))
    plt.gca().yaxis.set_major_formatter(FormatStrFormatter('%1.4f'))
    plt.savefig(figname)


#read in radius to cut
clust_id, rad_pix = get_rad('SBpixrad.txt')

#convert to degrees
rad = pix_to_deg(rad_pix)

#read in cluster centers from ap data table
ap_id, ap_ra, ap_dec, ap_rad = find_cen('apdata-cluster_6phot_v4.fits')

#loop over all clusters in SBpixrad file
no_ap_data=[]
for i in range(0, len(clust_id)):

    try:
        #read in cluster phot and fake stars
        path = '../ap' + str(clust_id[i])+ '/ap' + str(clust_id[i])
        ra, dec, phot_f475, phot_f814 = get_phot(path + '_phot.fits')
        fake_ra, fake_dec, phot_f475_in, phot_f814_in, phot_f475_out, phot_f814_out = get_fake(path + '.dst.fake.fits')

        #find table entry that matches clust_name
        place = np.where(ap_id == clust_id[i])	

        #find center of this cluster
        ra_cen = ap_ra[place]
        dec_cen = ap_dec[place]	


        if len(ra_cen > 0):
            #find distance from center for real and fake stars (in deg)
            dist = dist_from_cen(ra, dec, ra_cen, dec_cen)
            fake_dist = dist_from_cen(fake_ra, fake_dec, ra_cen, dec_cen)

            #keep only stars outside of cut radius
            in_stars = select_stars(rad[i], dist)
            fake_in_stars = select_stars(rad[i], fake_dist)

            #print number of total and kept stars
            print 'Total stars in cluster: ' +str(len(ra))
            print 'Stars kept: ' +str(len(ra[in_stars]))
            print 'Total fake stars: ' +str(len(fake_ra))
            print 'Fake stars kept: ' +str(len(fake_ra[fake_in_stars]))

            #plot all stars and stars kept
            plot_rad(rad[i], ra_cen, dec_cen, ra, dec, in_stars, 'ap' +str(clust_id[i])+ 'SBcut.png')
            plot_rad(rad[i], ra_cen, dec_cen, fake_ra, fake_dec, fake_in_stars, 'ap' +str(clust_id[i])+ 'SBcut_fake.png')

            #save kept stars to file
            print_stars(phot_f475[in_stars], phot_f814[in_stars], 'ap' +str(clust_id[i])+ '_phot.dat')
            print_fake_stars(phot_f475_in[fake_in_stars], phot_f814_in[fake_in_stars], phot_f475_out[fake_in_stars], phot_f814_out[fake_in_stars], 'ap' +str(clust_id[i])+ '_fake.dat')

            #append to file
            with open("SBcut_list.txt", "a") as myfile:
                myfile.write('ap' +str(clust_id[i])+'\n')

        else:
            no_ap_data.append(clust_id[i])

    except IOError:
        print path + '_phot.fits does not exist'










