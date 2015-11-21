import sys
import numpy as np
import pyfits

def main(name):
	"""writes console output to fits table
	inputs:  cluster name (string)
        output:  apXXX/apXXX_table.fits
	example call:  >>> import writefits
		       >>> writefits.main('ap1') """

        #read in output data from console file -- this file has 10 lines of header and 2 lines of footer
	matchoutput = np.genfromtxt(name+'/console_'+name+'.txt', skip_header=10, skip_footer=2)

        #calculate mass from the sfr and age columns
	mass = np.log10((matchoutput[:,6]*(10**(matchoutput[:,3]+0.1) - 10**matchoutput[:,3]))+1)
	n=len(mass)

        #create pyfits columns from matchoutput columns
	col1=pyfits.Column(name='Av', format='E', array=matchoutput[:,0])
	col2=pyfits.Column(name='IMF', format='E', array=matchoutput[:,1])	
	col3=pyfits.Column(name='dmod', format='E', array=matchoutput[:,2])
	col4=pyfits.Column(name='age', format='E', array=matchoutput[:,3])
	col5=pyfits.Column(name='Z', format='E', array=matchoutput[:,4])
	col6=pyfits.Column(name='fit', format='E', array=matchoutput[:,5])
	col7=pyfits.Column(name='SFR', format='E', array=matchoutput[:,6])
	col8=pyfits.Column(name='mass', format='E', array=mass)
	cols=pyfits.ColDefs([col1, col2, col3, col4, col5, col6, col7, col8])
	
        #create new table
	tbhdu=pyfits.new_table(cols)
	m=np.arange(n,8)
	hdu=pyfits.PrimaryHDU(m)
	thdulist=pyfits.HDUList([hdu,tbhdu])

	#make header with comment
 	with open(name+'/console_'+name+'.txt','r') as file:
  		data=file.readlines()
		headerdata=data[0:10]
	h1=headerdata[1].split(" ")
	list=[]
	for i in range(10):
		d=str(data[i])
		list.append(d)
	hstring= " ".join(list)
	hstring=hstring.replace("\n", "")
	prihdr=thdulist[0].header
	prihdr.add_comment(hstring)

        #write table to file
	thdulist.writeto(name+'/'+name+'_table.fits')

if __name__=='__main__':
        main(sys.argv[1])

