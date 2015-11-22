import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

'''plot jpg cluster image, with ap radius and bg radii

Reads in image apXXX_color.jpg  (not on github)

To run in python:
%run plot_jpg_rad.py ap81 1.91
1st argument:  clust name
2nd argument:  ap radius in arcseconds'''


def arcsec_to_pix(arcsec):
    '''convert arcseconds to pixels'''
    pix = arcsec / 0.05
    return pix


clust = sys.argv[1]
Rap_arcsec = float(sys.argv[2])

path = '../Paper3/'

img = mpimg.imread(path +clust + '_color.jpg')

#calc ap radius
Rap_pix = arcsec_to_pix(Rap_arcsec)

#calc bg radii
bg1 = arcsec_to_pix(1.2*Rap_arcsec)
bg2 = arcsec_to_pix(3.4*Rap_arcsec)

#define center of image
cen = 301/2.

#make array for circle
arr = np.linspace(0, 2*np.pi, 100)

#make plot
plt.imshow(img)
plt.plot(Rap_pix*np.cos(arr)+cen, Rap_pix*np.sin(arr)+cen, c='g', lw=2)
plt.plot(bg1*np.cos(arr)+cen, bg1*np.sin(arr)+cen, c='r', lw=2)
plt.plot(bg2*np.cos(arr)+cen, bg2*np.sin(arr)+cen, c='r', lw=2)
#plt.axhline(y=25/301., xmin=25/301., xmax=85/301., lw=2, color='white')
plt.axhline(y=15, xmin=0.1875, xmax=0.3375, lw=2, color='white')
plt.text(50, 35, '3"', size=15, color='white')


plt.show()





