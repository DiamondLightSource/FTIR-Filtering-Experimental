'''
Created on 9 Jul 2015

@author: flb41892
'''

from myClasses import fft 
from myClasses import fftfilter
import h5py
import pylab
import matplotlib.pyplot as plt
import pandas as pd
#############################################
#USER CODE BELOW

#VARIABLES

xl = pd.ExcelFile('/home/flb41892/data/ExcelFringes/' + '1.xlsx')

data = xl.parse('Sheet1', index_col = None, na_values = {'NA'})
data1 = xl.parse('a0', index_col = None, na_values = {'NA'})

interi = data.as_matrix()
interf = data1.as_matrix()

ind = interi[:,0]
initial = interi[:,1]
final =interf[:,1]


ft = fft()
ft2 = fftfilter()
f = h5py.File("/home/flb41892/data/fringe/Spectrum3.nxs","r") #file to analyse

#s = f["entry1/instrument/interferometer/sample_interferogram_changed_scan"][...] #signal on which to perform FFT and find the single channel spectrum
s = initial
res = f['/entry1/instrument/interferometer/opus_parameters/acquisition_changed/resolution'][...] #This extracts the resolution of your scan. 

highfold = f['/entry1/instrument/interferometer/opus_parameters/instrument_changed/high_folding_limit'][...]
# com = f["entry1/instrument/interferometer/sample_changed_scan"][...]# this is the FT of the same file, as performed by Opus commercial software

refer = f['/entry1/instrument/interferometer/reference_changed_scan'][...] #reference scan
renergy = f['/entry1/instrument/interferometer/reference_changed_energy'][...] # energy axis of reference scan 

fw = 1100 #defines the filter width. The filter is composed of a pair of symmetric, inverse Blackman Harris 3 (-67bD side lobe) filters, used to eliminate secondary fringes
#The filters are zero filled in the middle to sweep out a variable region of the spectrum. Select number of points you want to eliminate with the 'fw' parameter
fmin = 0.0 #the value of this parameter determines the minimum value of the filter, eg if fmin =0.0, some of the points will be completely eliminated. if fmin>0, points will be dampened only.
dv = 27.5 # dv = half of period of oscillatory fringes in absorbance spectrum/intensity spectrum, along the energy axis. 
#            NB.  Needs to be in units of cm-1!
#When you input dv, the program will use that information to position the inverse Blackman Harris filter to eliminate
#the oscillations in your spectrum. 


#######################################
#COMMANDS

#Non oscillatory data

#command that outputs 2 arrays, 1st array is the wavenumber axis [cm-1] and the second array is the single channel spectrum
schannel = ft.singleChannel(s,highfold) #use this function if you have a single, raw interferogram

schannel2 = ft.singleChannel2(s,highfold) #use this function if you have a double sided interferogram 
#NB. the high folding limit must be in cm-1

#absorb = ft.absorbance(schannel2, refer, renergy)
#########################################
#Oscillatory data

schanneloscil = ft2.singleChannel(s, fw, fmin,highfold,dv) #use this function if you have a single, oscillatory, raw interferogram
#         

schannel2oscil = ft2.singleChannel2(s, fw, fmin,highfold,dv) #use this function if you have a double sided, oscillatory, interferogram 
#NB. the high folding limit must be in cm-1
absorboscil = ft2.absorbance(schanneloscil, refer, renergy)


# example plotting tool below
pylab.plot(schannel2oscil[1],schannel2oscil[0],'g')
pylab.plot(schannel2[1],schannel2[0],'r')
pylab.xlabel("wavenumber [cm-1] ", fontsize = 17 )
pylab.ylabel("Intensity Spectrum [AU]", fontsize = 17)
plt.show()