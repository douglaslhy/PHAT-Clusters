#! /usr/bin/env python

#code from Alexia Lewis, based on IDL code from Cliff Johnson
#run_match.py calls this

import sys
import math
import numpy as np
import pylab as plt

from matplotlib import colors, cm
from sys import argv
from linecache import getline
from decimal import Decimal
from scipy.stats import mode
from pdb import set_trace

def main(name):

	'''plots residuals for pcXXX_out.cmd
	outputs pcXXX_out_WFC.png
	to call in ipython:  %run pg.py pc1
	or from cmd line:  python pg.py pc1'''

	def mycolormap(colorfile):
	    z = np.arange(0,256)
	    n = len(z)
	    z1 = min(z)
	    zn = max(z)
	    x0 = (z-z1)/float(zn-z1)

	    R,G,B = np.loadtxt(colorfile, unpack=True)
	    R = R[::-1]
	    G = G[::-1]
	    B = B[::-1]

	    cmap_dict = {}
	    cmap_dict['red'] = [(x0[i],R[i],R[i]) for i in range(len(R))]
	    cmap_dict['green']=[(x0[i],G[i],G[i]) for i in range(len(G))]
	    cmap_dict['blue'] =[(x0[i],B[i],B[i]) for i in range(len(B))]
	    mymap = colors.LinearSegmentedColormap('mymap',cmap_dict)
	    return mymap


	def plot_exclude_box(n, i, ax=None):
	    f = open(args.exclude, 'r')
	    lines = f.readlines()
	    f.close()

	    if n==1: ll = 7
	    if n==2: ll = 10
	    if n==3: ll = 13

	    exline = lines[ll+i].split()

	    if len(exline) > 2:
		x = [exline[1], exline[3], exline[5], exline[7], exline[1]]
		y = [exline[2], exline[4], exline[6], exline[8], exline[2]]
	    ##x = [1.25,5,5,1.25,1.25]
	    ##y = [ylow, ylow, yhigh, yhigh, ylow]
		ax.plot(x,y,'k--',lw=2)

	def save_plot(file, camera):
	    plotname = file+'_'+camera+'.png'
	    plt.savefig(plotname, dpi=100)
	    
	##############################################################################
	infile = name+'/'+name+'_out.cmd'
	ncmd = 1

	f = open(infile, 'r')
	lines = f.readlines()
	f.close()

	nums = lines[0].split()
	a1,a2,a3,a4,a5 = nums[0:5]

	line_n = 0
	for i_cmd in range(0,ncmd):
	    dims = lines[line_n+1].split()
	    xdim,ydim = int(dims[0]), int(dims[1])
	    nrows = xdim*ydim
	    col = lines[line_n+2].rstrip('\n')
	    mag = lines[line_n+3].rstrip('\n')

	    linelen = len(lines[line_n+4].split())
	    cmd1 = lines[(line_n+4):(line_n+4)+nrows]
	    cmd = np.zeros([nrows,linelen], dtype=float)

	    for i in range(0,len(cmd1)):
		newline = (cmd1[i].split())
		for j in range(0,len(newline)):
		    cmd[i][j] = float(newline[j])

	    #c1=mag, c2=color, c3=number observed, c4=number model, c5=difference, c6=significance
	    c1,c2,c3,c4,c5,c6,c7 = \
		       cmd[:,0],cmd[:,1],cmd[:,2],cmd[:,3],cmd[:,4],cmd[:,5],cmd[:,6]

	    line_n = (line_n+4)+nrows-1

	## ind = np.where(c7 < 0)
	## if np.count_zeros(c7 < 0) > 0:
	##    c3[ind] = 0.0
	##    c4[ind] = 0.0
	##    c5[ind] = 0.0
	##    c6[ind] = 0.0

	    xlab = col
	    ylab = mag

	    xmin = min(c2)
	    xmax = max(c2)
	    #ymin = min(c1)
	    ymin = 20.05
	    ymax = max(c1)
	    zmax = (max(c3))**0.4

	    img  = np.zeros((xdim,ydim),dtype=float)
	    img2 = np.zeros((xdim,ydim),dtype=float)
	    img3 = np.zeros((xdim,ydim),dtype=float)
	    img4 = np.zeros((xdim,ydim),dtype=float)

	    tot = len(c1)-1
	    
	    c5norm=max(abs(c5))-mode(c5)[0]
	    c6norm=max(abs(c6))-mode(c6)[0]

	    worstc5 = np.where(abs(c5) == max(abs(c5)))
	    worstc5pos = np.where(c5 == max(c5))
	    worstc5neg = np.where(c5 == min(c5))
	    worstc6 = np.where(abs(c6) == max(abs(c6)))
	    minc6 = np.where(c6 == min(c6))
	    maxc6 = np.where(c6 == max(c6))
	    worstdiff = max(abs(c5[worstc5]/c3[worstc5]))
	    worstdiffpos = max(c5[worstc5pos]/c3[worstc5pos])
	    worstdiffneg = max(c5[worstc5neg]/c3[worstc5neg])
	    worstfrac = max(c5[worstc6]/c3[worstc6])
	    minfrac   = max(c5[minc6]/c3[minc6])
	    maxfrac   = max(c5[maxc6]/c3[maxc6])
	    worstsig  = max(c6[worstc6])
	    worstsig  = max(c6[worstc6])

	    for j in range(0,ydim):
		for i in range(0,xdim):
		    element = ((j+1)+(i*ydim))-1
		    img[i,j] = (c3[element])**0.5
		    img2[i,j] = (c4[element])**0.5
		    img3[i,j] = c5[element] * (127/c5norm) + 127
		    img4[i,j] = c6[element] * (127/c6norm) + 127


	    fig = plt.figure()
	    plt.rc("axes", linewidth=1.5)

	    binary = cm.binary

	    cmap1 = binary
	    cmap2 = binary




	    ax1 = fig.add_subplot(221)
	    cnorm = colors.Normalize(vmin=(img.min())**2, vmax=(img.max())**2)
	    im = ax1.imshow(img[50:134,:], cmap=cmap1,aspect='auto',extent=(xmin,xmax,ymax,ymin),
		            interpolation=None, norm=cnorm)
	    ax1.imshow(img[50:134,:], cmap=cmap1, aspect='auto', extent=(xmin,xmax,ymax,ymin),
		       interpolation='nearest')
	    cbar = fig.colorbar(im, ax=ax1)
	    for t in cbar.ax.get_yticklabels():
		t.set_fontsize(12)
	    ax1.set_xlim(xmin,xmax)
	    ax1.set_ylim(ymax,ymin)
	    ax1.set_ylabel(ylab)
	    ax1.set_title('Observed', size=12)

	    
	    ax2 = fig.add_subplot(222, sharex=ax1, sharey=ax1)
	    cnorm  = colors.Normalize(vmin=(img2.min())**2, vmax=(img2.max())**2)
	    im = ax2.imshow(img2[50:134,:], cmap=cmap1, aspect='auto',
		            extent=(xmin,xmax,ymax,ymin), interpolation='nearest',
		            norm=cnorm)
	    ax2.imshow(img2[50:134,:], cmap=cmap1, aspect='auto',
		       extent=(xmin,xmax,ymax,ymin), interpolation='nearest')
	    cbar = fig.colorbar(im)
	    for t in cbar.ax.get_yticklabels():
		t.set_fontsize(12)
	    ax2.set_xlim(xmin,xmax)
	    ax2.set_ylim(ymax,ymin)
	    ax2.set_title('Modeled', size=12)

	    ax3 = fig.add_subplot(223, sharex=ax1, sharey=ax1)
	    cnorm = colors.Normalize(vmin=min(c5), vmax=max(c5))
	    im = ax3.imshow(img3[50:134,:], cmap=cmap2, aspect='auto',
		            extent=(xmin,xmax,ymax,ymin), interpolation=None,
		            norm=cnorm)
	    ax3.imshow(img3[50:134,:], cmap=cmap2, aspect='auto',
		       extent=(xmin,xmax,ymax,ymin), interpolation='nearest')
	    cbar = fig.colorbar(im)
	    for t in cbar.ax.get_yticklabels():
		t.set_fontsize(12)
	    ax3.set_xlim(xmin,xmax)
	    ax3.set_ylim(ymax,ymin)
	    ax3.set_xlabel(xlab)
	    ax3.set_ylabel(ylab)
	    ax3.set_title('Residuals', size=12)
	    
	    ax4 = fig.add_subplot(224, sharex=ax1, sharey=ax1)
	    cnorm = colors.Normalize(vmin=min(c6), vmax=max(c6))
	    im = ax4.imshow(img4[50:134,:], cmap=cmap2, aspect='auto',
		            extent=(xmin,xmax,ymax,ymin), interpolation=None,
		            norm=cnorm)
	    ax4.imshow(img4[50:134,:], cmap=cmap2, aspect='auto',
		       extent=(xmin,xmax,ymax,ymin), interpolation='nearest')
	    cbar = fig.colorbar(im)
	    for t in cbar.ax.get_yticklabels():
		t.set_fontsize(12)
	    ax4.set_xlim(xmin,xmax)
	    ax4.set_ylim(ymax,ymin)
	    ax4.set_xlabel(xlab)
	    ax4.set_title('Significance of Residuals', size=12)

	    
	    if col == 'WFC475W-WFC814W': cam='WFC'
	    if col == 'IR110W-IR160W': cam='IR'
	    if col == 'UVIS275W-UVIS336W': cam='UVIS'

	    save_plot(infile, cam)

if __name__=='__main__':
        main(sys.argv[1])
