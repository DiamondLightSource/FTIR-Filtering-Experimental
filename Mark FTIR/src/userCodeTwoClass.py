
from myClasses1 import fft 
from myClasses1 import fftfilter
import h5py
import pylab
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams.update( {'font.size':22})
#############################################
#USER CODE BELOW

#VARIABLES
ft = fft()
ft2 = fftfilter()
f = h5py.File("/home/flb41892/data/markmapping /linear result v2.hdf5","r") # load you nexus file here
f1 = h5py.File("/home/flb41892/data/markmapping /4 scan 3 00um.0.hdf5","r")
s = f["/SampleInterferogram"][...] #signal on which to perform FFT
#com = f["/Data/Sample/yfolded"][...]# this is the FT of the same file, as performed by opus
refer = f['/ReferenceInterferogram'][...] #reference scan

dim = np.shape(s)
spectra = np.zeros((dim[0],dim[1],389))
absorbance = np.zeros((dim[0],dim[1],389))
#reference = np.zeros((dim[0],dim[1],512))


for i in range(0,dim[0]):
    for j in range(0,dim[1]):
        s1 = s[i,j,:]
        #com1 = com[i,j,:]
        highfold = f1['/Parameters/Instrument/High Folding Limit'][...]
        
        zerofill = f1['/Parameters/FT/Zero Filling Factor'][...]
        zerofill =np.asarray(zerofill, int)
        
        refer1 = refer[i,j,:]
        #renergy =  f["entry1/instrument/interferometer/reference_energy"][...] # energy axis of reference scan 
        #absenergy = f["entry1/instrument/interferometer/ratio_absorbance_energy"][...] # energy axis of reference scan 
        #ymax = f['/Parameters/SampleDataInterferogram/Y - Maximum'][...] #max amplitude of interferogram processed by Opus
        ymax = s1.max()
        #yscaling = f['/Parameters/SampleDataInterferogram/Y - Scaling Factor'][...] #scaling factor that Opus applies to each intererigram before processing it.
        yscaling = 1.0
        ymaxspect = f1['/Parameters/SampleData/Y - Maximum'][...]#scaling factor that Opus applies to the final spectrum before plotting it.
        axis  = f1["/Data/Sample/x"][...] #energy axis from Opus
        
        fw = 40 #defines the filter width. The filter is composed of a pair of symmetric, inverse Blackman Harris 3 (-67bD side lobe) filters, used to eliminate secondary fringes
        #The filters are zero filled in the middle to sweep out a variable region of the spectrum. Select number of points you want to eliminate with the 'fw' parameter
        fmin = 0.0 #the value of this parameter determines the minimum value of the filter, eg if fmin =0.0, some of the points will be completely eliminated. if fmin>0, points will be dampened only.
        dv = 5 # dv = half of period of oscillatory fringes in absorbance spectrum/intensity spectrum, along the energy axis. 
        #            NB.  Needs to be in units of cm-1!
        #When you input dv, the program will use that information to position the inverse Blackman Harris filter to eliminate
        #the oscillations in your spectrum. 
        ymaxinterf = s.max() #compute the max height of interferogram (we ll need to pass it as a parameter)
        
        #######################################
        #COMMANDS
        
        #Non oscillatory, single sided data
        
        # single channel function that outputs 2 arrays, 1st array is the  single channel spectrum and the second array is the wavenumber axis [cm-1]
        
        #schannel = ft.singleChannel(s1,highfold,zerofill,ymax,ymaxinterf,yscaling,ymaxspect) #use this function if you have a single, raw interferogram
        #absorb = ft.absorbance(schannel, refer1, highfold,zerofill,ymax,ymaxinterf,yscaling,ymaxspect)
        
        #Non oscillatory, DOUBLE sided data
        
        schannel2 = ft.singleChannel2(s1,highfold,zerofill,ymax,ymaxinterf,yscaling,ymaxspect) #use this function if you have a double sided interferogram 
        #NB. the high folding limit must be in cm-1
        
        absorb2 = ft2.absorbance2(schannel2, refer1, highfold,zerofill,ymax,ymaxinterf,yscaling,ymaxspect) # absorbance function for a double sided sample and reference interferogram
        #########################################
        #Oscillatory data
        
        #Oscillatory, single sided data.
        
        
        #schanneloscil = ft2.singleChannel(s1, fw, fmin,highfold,dv,zerofill,ymax,ymaxinterf,yscaling,ymaxspect) #use this function if you have a single-sided, oscillatory, raw interferogram
        #absorboscil = ft.absorbance(schanneloscil, refer1, highfold,zerofill,ymax,ymaxinterf,yscaling,ymaxspect)
        
        #Oscillatory, double sided data
        
        
        schannel2oscil = ft2.singleChannel2(s1, fw, fmin,highfold,dv,zerofill,ymax,ymaxinterf,yscaling,ymaxspect) #use this function if you have a double sided, oscillatory, interferogram 
        #NB. the high folding limit must be in cm-1
        absorboscil = ft2.absorbance2(schannel2oscil, refer1, highfold,zerofill,ymax,ymaxinterf,yscaling,ymaxspect)
        
        
        
        # example plotting tool below
        x = schannel2 # select which spectrum to plot
        
        #Tool for selecting which region of spectrum to display (default is 500 to 4000cm-1)
        a = np.add(x[1],-axis[0])
        a = np.abs(a)
        a = a.argmin()
        b = np.add(x[1],-axis[axis.size-1])
        b = np.abs(b)
        b = b.argmin()
        final = x[0][a-1:b+1]
        energy = x[1][a-1:b+1]
        absorb2 = absorb2[a-1:b+1]
        #referen = absorb2[1][a-1:b+1] 
        #final = final * com1.max()/final.max()
        spectra[i][j] = final
        absorbance[i][j] = absorb2
        #reference[i][j] = referen


#pylab.plot(schannel2oscil[1],schannel2oscil[0])
pylab.plot(energy,final, label = 'our spectrum')
pylab.plot(energy,com, label = 'Opus spectrum')
pylab.xlabel("wavenumber [cm-1] ", fontsize = 21 )
pylab.ylabel("Intensity Spectrum [AU]", fontsize = 21)
pylab.legend(loc = 'upper right')

plt.show()