import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
import scipy.optimize as so
import matplotlib.cm as cm
import pyfits
from matplotlib import rc
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

'''comparison plots for syn clusters using synindex_v3.fits and synindex_6band.dat'''

#read in input data
input = pyfits.open('synindex_v3.fits')
input_data = input[1].data
synid1 = input_data.field('synid')
batch1 = input_data.field('syniter')
name1 = input_data.field('name')
in_age1 = input_data.field('logage')
in_mass1 = input_data.field('logmass')
in_av1 = input_data.field('av')

synid2 = np.loadtxt('synindex_6band.dat', usecols=[0]).astype(int)
batch2 = np.loadtxt('synindex_6band.dat', usecols=[1]).astype(int)
name2 = np.genfromtxt('synindex_6band.dat', usecols=(2), dtype=None)
in_age2 = np.loadtxt('synindex_6band.dat', usecols=[7]).astype(float)
in_mass2 = np.loadtxt('synindex_6band.dat', usecols=[8]).astype(float)
in_av2 = np.loadtxt('synindex_6band.dat', usecols=[9]).astype(float)

#concatenate arrays
synid = np.concatenate([synid1, synid2])
batch = np.concatenate([batch1, batch2])
name = np.concatenate([name1, name2])
in_age = np.concatenate([in_age1, in_age2])
in_mass = np.concatenate([in_mass1, in_mass2])
in_av = np.concatenate([in_av1, in_av2])

#create arrays for match results
bf_age = np.zeros(len(name))
bf_mass = np.zeros(len(name))
bf_av = np.zeros(len(name))




old_list=[]

#loop over input data to find match results
field=np.zeros(len(name)).astype('str')

for i in range(len(name)):
    field[i] = name[i].split('-')[1]

    try:
        if field[i] == 'F01':
            clust = name[i].split('_')[0]
            best_fit = np.loadtxt('B17-F01-' +str(batch[i])+'/' +clust+ '/' +clust+ '_bestfit.txt')

            bf_age[i] = best_fit[0]
            bf_av[i] = best_fit[1]
            bf_mass[i] = best_fit[2]

        elif field[i] == 'F11':

            clust = name[i].split('_')[0]
            best_fit = np.loadtxt('B06-F11-' +str(batch[i])+'/' +clust+ '/' +clust+ '_bestfit.txt')

            bf_age[i] = best_fit[0]
            bf_av[i] = best_fit[1]
            bf_mass[i] = best_fit[2]

    except IOError:
        print str(name[i])+ '_bestfit.txt does not exist'



def plot_ratio(x, num, denom, xlim, ylim, xlabel, ylabel, title):
	'''plot x vs num/denom'''

	plt.scatter(x, num/denom)
	plt.axhline(y=1, color='black', linestyle='--', lw=3)
	plt.xlim(xlim)
	plt.ylim(ylim)
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	plt.title(title)

def plot_diff(x, y1, y2, xlim, ylim, xlabel, ylabel, title):
	'''plot x vs y1-y2'''

	plt.scatter(x, y1-y2)
	plt.axhline(y=0, color='black', linestyle='--', lw=3)
	plt.xlim(xlim)
	plt.ylim(ylim)
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	plt.title(title)


def scatter_plot(subplot, xlim1, xlim2, X, Y, color, xlabel, ylabel, title, label, loc):
    '''plot scatter plot w/ density contours
    inputs:  
	subplot:  number of subplot, eg 321
	xlim1, xlim2:  x limits for plotting
	X, Y:  X and Y data to be plotted
	color:  color array for points
	xlable, ylabel, title:  string labels
    '''

    ax=fig.add_subplot(subplot, aspect='equal')


    eq1 = np.linspace(xlim1, xlim2, 10)
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    #plt.scatter(X, Y, c=color, s=15., edgecolors='none', cmap=cm.cool, alpha=0.7)
    plt.scatter(X, Y, c=color, s=15., edgecolors='none', cmap=cm.cool, alpha=0.5, label=label)
    plt.plot(eq1, eq1, linestyle='--', color='black')
    plt.ylim(xlim1, xlim2)
    plt.xlim(xlim1, xlim2)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.xlabel(xlabel, fontsize=20)
    plt.ylabel(ylabel, fontsize=20)
    #plt.title(title)
    #plt.colorbar(pad=0)
    plt.legend(loc=loc, scatterpoints=1, prop={'size':11})




def make_hist(subplot, bins, D1, D2, label1, label2, xlim, ylim, xlabel, ylabel):
    '''make comparison histograms for 2 distributions
    inputs:
        subplot:  number of subplot, eg 323
        bins:  numpy array of bin edges, eg, np.arange(0, 10.5, 0.33)
	D1, D2:  distributions to plot
	label1, label2:  string of labels for legend
	xlim, ylim:  list of limits
	xlabel, ylabel:  string for label
    '''
    ax=fig.add_subplot(subplot)
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.hist(D1, bins=bins, histtype='stepfilled',label=label1)
    plt.hist(D2, bins=bins, histtype='stepfilled', alpha=0.5,color='r',label=label2)
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.xlabel(xlabel, fontsize=20)
    plt.ylabel(ylabel, fontsize=20)
    plt.legend()



def grid(x, y, z, resX=100, resY=100):
    xi=np.linspace(min(x), max(x), resX)
    yi = np.linspace(min(y), max(y), resY)
    Z=griddata(x,y,z,xi,yi)
    X,Y= np.meshgrid(xi,yi)
    return X,Y,Z


def find_confidence_interval(x, pdf, confidence_level):
    return pdf[pdf > x].sum() - confidence_level
 
def density_contour(subplot, xlim1, xlim2, xdata, ydata, nbins_x, nbins_y, xlabel, ylabel, title, **contour_kwargs):
    """ Create a density contour plot.

    Parameters
    ----------
    subplot : 3 digit number
    xlim1, xlim2 : xlimits
    xdata : numpy.ndarray
    ydata : numpy.ndarray
    nbins_x : int
    Number of bins along x dimension
    nbins_y : int
    Number of bins along y dimension
    xlabel, ylabel, title : strings
    ax : matplotlib.Axes (optional)
    If supplied, plot the contour to this axis. Otherwise, open a new figure
    contour_kwargs : dict
    kwargs to be passed to pyplot.contour()
    """

    ax=fig.add_subplot(subplot)
    eq1 = np.linspace(xlim1, xlim2, 10)

    H, xedges, yedges = np.histogram2d(xdata, ydata, bins=(nbins_x,nbins_y), normed=True)
    x_bin_sizes = (xedges[1:] - xedges[:-1]).reshape((1,nbins_x))
    y_bin_sizes = (yedges[1:] - yedges[:-1]).reshape((nbins_y,1))

    pdf = (H*(x_bin_sizes*y_bin_sizes))

    one_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.68))
    two_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.95))
    three_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.99))
    levels = [one_sigma, two_sigma]

    X, Y = 0.5*(xedges[1:]+xedges[:-1]), 0.5*(yedges[1:]+yedges[:-1])
    Z = pdf.T

    plt.scatter(xdata, ydata, c='b',s=2.)
    plt.plot(eq1, eq1, linestyle='--', color='black')

    plt.contourf(X, Y, Z, levels=levels, origin="lower", **contour_kwargs)

    plt.ylim(xlim1, xlim2)
    plt.xlim(xlim1, xlim2)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    #plt.colorbar(pad=0)


good = np.where(bf_mass > 0)
name_good = name[good]
in_age=in_age[good]
bf_age=bf_age[good]
in_mass=in_mass[good]
bf_mass=bf_mass[good]
in_av=in_av[good]
bf_av=bf_av[good]

plt.close('all')
fig=plt.figure()
#scatter_contour(231, 6.5, 9, 10, in_age[good] + np.random.normal(0, 0.025, len(in_age[good])), bf_age[good], [10,20,30], '$\log_{10}$(input age)', '$\log_{10}$(match age)', 'age comparison')
#scatter_contour(232, 2.5, 4, 10, in_mass[good] + np.random.normal(0, 0.025, len(in_mass[good])), bf_mass[good], [10,20,30], '$\log_{10}$(input mass)', '$\log_{10}$(match mass)', 'mass comparison')
#scatter_contour(233, -0.1, 1.0, 10, in_av[good] + np.random.normal(0, 0.025, len(in_av[good])), bf_av[good], [80,100,120,140,160], 'input av', 'match av', 'av comparison')

#density_contour(231, 6.5, 9, in_age[good], bf_age[good], 10, 10, '$\log_{10}$(input age)', '$\log_{10}$(match age)', 'age comparison')
#density_contour(232, 2.5, 4, in_mass[good], bf_mass[good], 10, 10, '$\log_{10}$(input mass)', '$\log_{10}$(match mass)', 'mass comparison')
#density_contour(233, -0.1, 1.0, in_av[good], bf_av[good], 10, 10, 'input av', 'match av', 'av comparison')


avdiff=in_av-bf_av
#color=avdiff[good]
color=in_mass

#make four mass bins for color-code
bin1 = np.where((in_mass > 0) & (in_mass < 3))	 	#purple
bin2 = np.where((in_mass >= 3) & (in_mass < 3.5)) 	#blue
bin3 = np.where((in_mass >= 3.5) & (in_mass < 4)) 	#green
bin4 = np.where((in_mass >= 4) & (in_mass < 5)) 	#red
bin5 = np.where(in_mass >= 5)				#orange

label1 = r'$\log_{10}$(M/$M_{\odot}$) \textless\ 3'
label2 = r'$\log_{10}$(M/$M_{\odot}$) \textless\ 3.5'
label3 = r'$\log_{10}$(M/$M_{\odot}$) \textless\ 4'
label4 = r'$\log_{10}$(M/$M_{\odot}$) \textless\ 4.5'

scatter_plot(231, 6.0, 10.0, in_age[bin1] + np.random.normal(0, 0.03, len(in_age[bin1])), bf_age[bin1] + np.random.normal(0, 0.03, len(in_age[bin1])), 'purple', 'Input $\log_{10}$(t [yr])', 'MATCH $\log_{10}$(t [yr])', ' ', label1, 4)
scatter_plot(232, 2.0, 4.6, in_mass[bin1] + np.random.normal(0, 0.03, len(in_age[bin1])), bf_mass[bin1] + np.random.normal(0, 0.03, len(in_mass[bin1])), 'purple', 'Input $\log_{10}$(M [$M_{\odot}$])', 'MATCH $\log_{10}$([$M_{\odot}$])', ' ', label1, 4)
scatter_plot(233, -0.1, 2.2, in_av[bin1] + np.random.normal(0, 0.03, len(in_age[bin1])), bf_av[bin1] + np.random.normal(0, 0.03, len(in_av[bin1])), 'purple', 'Input $A_{V}$ [mag]', 'MATCH $A_{V}$ [mag]', ' ', label1, 2)

scatter_plot(231, 6.0, 10.0, in_age[bin2] + np.random.normal(0, 0.03, len(in_age[bin2])), bf_age[bin2] + np.random.normal(0, 0.03, len(in_age[bin2])), 'blue', 'Input $\log_{10}$(t [yr])', 'MATCH $\log_{10}$(t [yr])', ' ', label2, 4)
scatter_plot(232, 2.0, 4.6, in_mass[bin2] + np.random.normal(0, 0.03, len(in_age[bin2])), bf_mass[bin2] + np.random.normal(0, 0.03, len(in_mass[bin2])), 'blue', 'Input $\log_{10}$(M [$M_{\odot}$])', 'MATCH $\log_{10}$([$M_{\odot}$])', ' ', label2, 4)
scatter_plot(233, -0.1, 2.2, in_av[bin2] + np.random.normal(0, 0.03, len(in_age[bin2])), bf_av[bin2] + np.random.normal(0, 0.03, len(in_av[bin2])), 'blue', 'Input $A_{V}$ [mag]', 'MATCH $A_{V}$ [mag]', ' ', label2, 2)

scatter_plot(231, 6.0, 10.0, in_age[bin3] + np.random.normal(0, 0.03, len(in_age[bin3])), bf_age[bin3] + np.random.normal(0, 0.03, len(in_age[bin3])), 'green', 'Input $\log_{10}$(t [yr])', 'MATCH $\log_{10}$(t [yr])', ' ', label3, 4)
scatter_plot(232, 2.0, 4.6, in_mass[bin3] + np.random.normal(0, 0.03, len(in_age[bin3])), bf_mass[bin3] + np.random.normal(0, 0.03, len(in_mass[bin3])), 'green', 'Input $\log_{10}$(M [$M_{\odot}$])', 'MATCH $\log_{10}$([$M_{\odot}$])', ' ', label3, 4)
scatter_plot(233, -0.1, 2.2, in_av[bin3] + np.random.normal(0, 0.03, len(in_age[bin3])), bf_av[bin3] + np.random.normal(0, 0.03, len(in_av[bin3])), 'green', 'Input $A_{V}$ [mag]', 'MATCH $A_{V}$ [mag]', ' ', label3, 2)

scatter_plot(231, 6.0, 10.0, in_age[bin4] + np.random.normal(0, 0.03, len(in_age[bin4])), bf_age[bin4] + np.random.normal(0, 0.03, len(in_age[bin4])), 'red', 'Input $\$\log_{10}$_{10}$(t [yr])', 'MATCH $\log_{10}$(t [yr])', ' ', label4, 4)
scatter_plot(232, 2.0, 4.6, in_mass[bin4] + np.random.normal(0, 0.03, len(in_age[bin4])), bf_mass[bin4] + np.random.normal(0, 0.03, len(in_mass[bin4])), 'red', 'Input $\log_{10}$(M [$M_{\odot}$])', 'MATCH $\log_{10}$([$M_{\odot}$])', ' ', label4, 4)
scatter_plot(233, -0.1, 2.2, in_av[bin4] + np.random.normal(0, 0.03, len(in_age[bin4])), bf_av[bin4] + np.random.normal(0, 0.03, len(in_av[bin4])), 'red', 'Input $A_{V}$ [mag]', 'MATCH $A_{V}$ [mag]', ' ', label4, 2)



make_hist(234, np.arange(0,10.0,0.33), in_age, bf_age, 'Input $\log_{10}$(t [yr])', 'MATCH $\log_{10}$(t [yr])', [6, 10], [0, 300], '$\log_{10}$(t [yr])', 'Number')
make_hist(235, np.arange(2.2,4.6,0.2), in_mass, bf_mass, 'Input $\log_{10}$(M [$M_{\odot}$])', 'MATCH $\log_{10}$([$M_{\odot}$])', [2.2, 4.6], [0, 520], '$\log_{10}$([$M_{\odot}$])', 'Number')
make_hist(236, np.arange(0,2,0.2), in_av, bf_av, 'Input $A_{V}$ [mag]', 'MATCH $A_{V}$ [mag]', [0, 2.2], [0, 600], '$A_{V}$ [mag]', 'Number')


plt.show()

sqdiff=(in_age-bf_age)**2.
sum_sqdiff=np.sum(sqdiff)
rms=(sum_sqdiff/len(in_age))**0.5

######################################################################################


#plt.close('all')
#fig=plt.figure()

#ax=fig.add_subplot(131)
#plot_ratio(bf_age[good], in_age[good], bf_age[good], [6.5,9.0], [0.7, 1.3], 'match age', 'input age / match age', 'Age Comparison')
#plot_ratio(bf_age[good], in_mass[good], bf_mass[good], [6.5,9.0], [0.7, 1.4], 'match age', 'input mass / match mass', 'Age and Mass Comparison')
#plot_diff(bf_age[good], in_av[good], bf_av[good], [6.5,9.0], [-1, 1], 'match age', 'input av - match av', 'Age and Av Comparison')

#ax=fig.add_subplot(132)
#plot_ratio(bf_mass[good], in_age[good], bf_age[good], [1.5,5.5], [0.6, 1.4], 'match mass', 'input age / match age', 'Age and Mass Comparison')
#plot_ratio(bf_mass[good], in_mass[good], bf_mass[good], [2.0,5.0], [0.7, 1.4], 'match mass', 'input mass / match mass', 'Mass Comparison')
#plot_diff(bf_mass[good],  in_av[good], bf_av[good], [2.0,5.0], [-1, 1], 'match mass', 'input av - match av', 'Mass and Av Comparison')

#ax=fig.add_subplot(133)
#plot_ratio(bf_av[good], in_age[good], bf_age[good], [0.,1.0], [0.6, 1.4], 'match av', 'input age / match age', 'Age and Av Comparison')
#plot_ratio(bf_av[good], in_mass[good], bf_mass[good], [0.,1.0], [0.7, 1.4],'match av', 'input mass / match mass',  'Av and Mass Comparison')
#plot_diff(bf_av[good],  in_av[good], bf_av[good], [0.,1.0], [-1, 1], 'match av', 'input av - match av', 'Av Comparison')

