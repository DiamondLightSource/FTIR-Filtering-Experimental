'''
Created on 7 Jul 2015

@author: flb41892
'''
import numpy as np
import cmath as m

class fft():
    def __init__(self):
        pass
    def _zerofill(self,single,highfold,zerofill,ymax,ymaxinterf,yscaling):
        single = single - np.mean(single) # eliminate offset of interferogram
        single = single*yscaling*10.0/(ymaxinterf/ymax)
        if highfold <7900.0:
            if 16384<single.size < 32768:
            
                    if zerofill < 4:
                        single = np.concatenate((single,np.zeros(32768-single.size)))
                    if zerofill == 4:
                        single = np.concatenate((single,np.zeros(65536-single.size)))
            
            if 8192<single.size < 16384:
            
                    if zerofill < 4:
                        single = np.concatenate((single,np.zeros(16384-single.size)))
                    if zerofill == 4:
                        single = np.concatenate((single,np.zeros(32768-single.size)))
            
            if 4096<single.size < 8192:
                    if zerofill < 4:
                        single = np.concatenate((single,np.zeros(8192-single.size)))
                    if zerofill == 4:
                        single = np.concatenate((single,np.zeros(16384-single.size)))
            if single.size < 4096:
                    if zerofill < 4:
                        single = np.concatenate((single,np.zeros(4096-single.size)))
                    if zerofill == 4:
                        single = np.concatenate((single,np.zeros(8192-single.size)))
            single = single*2
        if 7900.0<highfold <15800.0:
            if 16384<single.size < 32768:
            
                    if zerofill < 4:
                        single = np.concatenate((single,np.zeros(32768-single.size)))
                    if zerofill == 4:
                        single = np.concatenate((single,np.zeros(65536-single.size)))
            
            if 8192<single.size < 16384:
            
                    if zerofill < 4:
                        single = np.concatenate((single,np.zeros(16384-single.size)))
                    if zerofill == 4:
                        single = np.concatenate((single,np.zeros(32768-single.size)))
            
            if 4096<single.size < 8192:
                    if zerofill < 4:
                        single = np.concatenate((single,np.zeros(8192-single.size)))
                    if zerofill == 4:
                        single = np.concatenate((single,np.zeros(16384-single.size)))
            
            if single.size < 4096:
                    if zerofill < 4:
                        single = np.concatenate((single,np.zeros(4096-single.size)))
                    if zerofill == 4:
                        single = np.concatenate((single,np.zeros(8192-single.size)))
            

        return single
    def _mertz(self,single):
        n = 256 # number of points to select for phase correction about ZPD point
        zeros = single[single.argmax()-n:single.argmax()+n]
        #ramp function (ramp is old triangular fcn, better to use the Black Harris 3 step fct w[t] )
        ramp = np.zeros(2*n)
        ramp[0:n] = np.linspace(0,1,n,endpoint=False)
        ramp[n:] = np.linspace(1,0,n)
        zeros = zeros*ramp #multiply zeros array by ramp fcn to prepare array for phase correction
        #rearrange data, so that right side of data(including peak) is moved to front of array and left hand side
        #is moved to the back
        
        #rotate the 512 long array
        interf = []
        interf[0:n] = zeros[np.argmax(zeros):zeros.size]
        
        interf[n:]=zeros[0:np.argmax(zeros)]
        ones = np.ones(np.size(interf))
        
        ones[25:60] = np.linspace(0.5,0,35, endpoint = False)
        ones[460:500] = np.linspace(0,0.5, 40)
        interf1 = interf * ones
        #zerofill central region to eliminate secondary peak
#         z = 50
#         interf1[np.size(interf1)/2-z:np.size(interf1)/2+z] = np.zeros(2*z)
        #fourier transform 
        output_axis= np.size(interf)
        trans= np.fft.fft(interf, output_axis)
        #reff= np.fft.rfft(ref,output_axis)
        #decompose into real and imaginary parts of fourier spectrum
        freal= np.real(trans)
        fim= np.imag(trans)
        
        #reffr = np.abs(reff)#do same with reference set
        #calculate phase angle
        phi = np.arctan(np.divide(fim,freal))
        cphi = np.cos(phi)
        sphi = np.sin(phi)
        #extend phase arrays to match interferogram arrays(interpolation)
        xp = np.arange(0,2*n)
        x = np.arange(0,2*n,512./single.size)
        cphi2 = np.interp(x,xp,cphi)
        sphi2 = np.interp(x,xp,sphi)
        self.zeros = zeros
        self.interf = interf
        self.interf1 = interf1
        self.ones = ones
        self.freal = freal
        self.angle = cphi2

        return cphi2,sphi2
    def _apod(self,single0): #Black Harris 3 term apodization fcn.
        apodf = np.zeros(single0.size) #61dB
        apodf2 = np.zeros(single0.size) #67 dB 
        for j in range(0,single0.size):
            apodf[j] = 0.44959-.49364*np.cos(2*m.pi*j/single0.size)+.05677*np.cos(4*m.pi*j/single0.size)
            apodf2[j] = 0.42323-.49755*np.cos(2*m.pi*j/single0.size)+.07922*np.cos(4*m.pi*j/single0.size)
        ins = ((np.size(single0)-np.argmax(single0)) - (np.argmax(single0)))/2
        single0 = np.insert(single0,0,np.zeros(ins))
        single0 = single0[:np.size(single0)-ins]
        single0 = single0 *apodf2
        apod_singler = np.zeros(np.size(single0))
        apod_singler[0:single0.size-np.argmax(single0)] = single0[np.argmax(single0):single0.size]
        apod_singler[single0.size-np.argmax(single0):]=single0[0:np.argmax(single0)]
        self.apodf2 = apodf2
        self.single0 = single0
        return apod_singler
    def _fft(self,filtered,angles):
        output_axis1 = np.size(filtered)
        apodfi = np.fft.fft(filtered, output_axis1)
        apodr = np.real(apodfi)
        apodi = np.imag(apodfi)
        finalr = np.multiply(apodr,angles[0]) 
        finali = np.multiply(apodi,angles[1])
        final = np.add(finalr,finali)
        self.apodr = apodr
        return final
    def singleChannel(self,s,highf,zerofill,ymax,ymaxinterf,yscaling):
        ft =fft()
        """
        ###############################
        Converts a SINGLE SIDED INTERFEROGRAM to a single channel spectrum.
        
        It processes a single interferogram.
        
        Takes in 2 arguments :
        
        s = single sided sample interferogram (with or without secondary fringes)
        
        highf = high frequency folding limit (parameter of the scan), can usually be found at 
        '/entry1/instrument/interferometer/opus_parameters/instrument_changed/high_folding_limit'
        
        Outputs:
        2D array
        
        0th array = intensity spectrum:
        Single channel spectrum is computed from the original interferogram, which was zerofilled to next factor of 2,
        apodized using a 3 term Blackman Harris function and phase corrected following the Mertz method.
        
        1st array= associated wavenumber axis.
        
        
        """
        
        
        single = s[:] #in case of bifringent interferogram, take only one peak to analyse (avoids sinusoidal modulations)
        #zero filling(pad until 16,384 if array is below this number and up to 65536 points if array is larger)
        
        single0 = ft._zerofill(single,highf,zerofill,ymax,ymaxinterf,yscaling)
        self.sing = single0
        angles = ft._mertz(single0)
        
        single0 = ft._apod(single0)
        self.sing1 = single0
        
        schannel = ft._fft(single0,angles)
        
        #calculate the axis in frequency space
        #frequency axis
        lmda = 1/highf#cm (this is the laser waveleght you re using)
        k = np.arange(np.size(single0))
        v = np.divide(2*k,lmda*np.size(single0)) # f = k/(N*lambda) where k is range of values from zero to array size,
        self.dv = v[100]-v[99]
        kbig = np.arange(np.size(single0))
        vbig = np.divide(kbig,lmda*np.size(single0))
        
        return schannel,v
    def singleChannel2(self,s,highf,zerofill,ymax,ymaxinterf,yscaling):
        ft =fft()
        """
        ###############################
        Converts a DOUBLE SIDED INTERFEROGRAM to a single channel spectrum.
        
        It processes the 2 single interferograms separately and averages them at the end.
        
        Takes in 2 arguments :
        
        s = double sided sample interferogram (with or without secondary fringes)
        
        highf = high frequency folding limit (parameter of the scan), can usually be found at 
        '/entry1/instrument/interferometer/opus_parameters/instrument_changed/high_folding_limit'
        
        Outputs:
        2D array
        
        0th array = intensity spectrum:
        Single channel spectrum is computed from the original interferogram, which was zerofilled to next factor of 2,
        apodized using a 3 term Blackman Harris function and phase corrected following the Mertz method.
        
        1st array= associated wavenumber axis.
        
        
        """
        
        
        single = s[0.5*s.size:] #in case of bifringent interferogram, take only one peak to analyse (avoids sinusoidal modulations)
        #zero filling(pad until 16,384 if array is below this number and up to 65536 points if array is larger)
        
        single0 = ft._zerofill(single,highf,zerofill,ymax,ymaxinterf,yscaling)
        
        angles = ft._mertz(single0)
        
        single0 = ft._apod(single0)
        
        schannel = ft._fft(single0,angles)
        #calculate the axis in frequency space
        #frequency axis
        lmda = 1/highf#cm (this is the laser waveleght you re using)
        k = np.arange(np.size(single0))
        v = np.divide(2*k,lmda*np.size(single0)) # f = k/(N*lambda) where k is range of values from zero to array size,
        self.dv = v[100]-v[99]
        kbig = np.arange(np.size(single0))
        vbig = np.divide(kbig,lmda*np.size(single0))
        
        schannel2 = ft.singleChannel(s[0.5*s.size:], highf,zerofill,ymax,ymaxinterf,yscaling)
        final = np.add(schannel,schannel2[0])
        final = np.true_divide(final,2)
        self.single2 = single0
                

        return final,v
    def absorbance(self, schannel, refer,highfold,zerofill,ymax,ymaxinterf,yscaling):
        """
        #############################
        User Function absorbance
        NB. WORKS ONLY IF A REFERENCE INTERFEROGRAM IS PROVIDED.
        
        inputs 4 arguments:
        
        schannel = single channel spectrum (no axis necessary, just the data array)
        refer = the reference, single sided interferogram (cannot calculate absorbance if reference interferogram has not been taken)
        highfold = high folding limit 
        zerofill = zerofill factor (eg. 1,2 or 4)
        
        returns:
        
        1D array containing the absorbance spectrum. 
        The energy axis is the same as the single channel axis. (in k numbers)
        
        """
        ft =fft()
        refer = ft.singleChannel(refer,highfold,zerofill,ymax,ymaxinterf,yscaling)
        absorbance = -np.log10(schannel[0]/refer[0])
        return absorbance[0:absorbance.size/2]
class fftfilter(fft):
    def __init__(self):
        pass
    def _filter(self,apod_singler,c,m,dv,highf):

        apod_singler2 = apod_singler
        #find position of the secondary fringe
        dz = 1./dv
        dn = np.round(dz*highf)
        #implement black harris inverse filter 
        a = dn - c/2
        b = dn + c/2
        blh = np.ones(np.size(apod_singler)-c) #remove 2c values from array, because we will add 2c zeros later

        blh[a:b] = np.add(np.ones(b-a),-np.blackman(b-a))
        
        blh = np.insert(blh,np.argmin(blh),np.multiply(m,np.ones(c)))
        c = blh[:np.size(blh)/2]
        c = c[::-1]
        blh[np.size(blh)/2:] = c
        apod_singler2 = np.multiply(apod_singler2,blh)
        self.blh = blh
        self.ap = apod_singler2
        return apod_singler2
    def singleChannel(self,s,fw,fmin,highf,dv,zerofill,ymax,ymaxinterf,yscaling):
        ft =fft()
        ft2=fftfilter()
        """
        ###############################
        Converts a SINGLE INTERFEROGRAM to a single channel spectrum.
        
        Takes in 5 arguments :
        
        s = single sample interferogram (with or without secondary fringes)
        
        fw = filter width, should be zero if interferogram is good quality.
        
        fmin = minimum filter height  (0<fmin<1)
        
        fmin = 0, eliminate all points within filter range.
        fmin = 1, full pass filter.
        
        highf = high frequency folding limit (parameter of the scan), can usually be found at 
        '/entry1/instrument/interferometer/opus_parameters/instrument_changed/high_folding_limit'
        
        dv = half of period of oscillatory fringes in absorbance spectrum/intensity spectrum, along the energy axis. 
           NB.  Needs to be in units of cm-1!
        
        When you input dv, the program will use that to position an inverse Blackman Harris filter to eliminate
        the oscillations in your spectrum. 
        
        Outputs:
        2D array
        
        0th array = single channel spectrum:
        Single channel spectrum is computed from the original interferogram, which was zerofilled to next factor of 2,
        apodized using a 3 term Blackman Harris function and phase corrected following the Mertz method.
        
        1st array= associated wavenumber axis.
        
        
        """
        
        
        single = s[:] #in case of bifringent interferogram, take only one peak to analyse (avoids sinusoidal modulations)
        #zero filling(pad until 16,384 if array is below this number and up to 65536 points if array is larger)
        
        single0 = ft._zerofill(single,highf,zerofill,ymax,ymaxinterf,yscaling)
        
        angles = ft._mertz(single0)
        
        single0 = ft._apod(single0)
        
        filtered= ft2._filter(single0, fw, fmin,dv,highf)
        
        schannel = ft._fft(filtered,angles)
        
        #calculate the axis in frequency space
        #frequency axis
        lmda = 1/highf#cm (this is the laser waveleght you re using)
        k = np.arange(np.size(single0))
        v = np.divide(2*k,lmda*np.size(single0)) # f = k/(N*lambda) where k is range of values from zero to array size,
        self.dv = v[100]-v[99]
        self.single0 = single0
        kbig = np.arange(np.size(single0))
        vbig = np.divide(kbig,lmda*np.size(single0))
#         schannel = schannel/schannel.max()
        return schannel,v
    def singleChannel2(self,s,fw,fmin,highf,dv,zerofill,ymax,ymaxinterf,yscaling):
        ft =fft()
        ft2 = fftfilter()
        """
        ###############################
        Converts a DOUBLE SIDED INTERFEROGRAM to a single channel spectrum.
        
        It processes the 2 single interferograms separately and averages them at the end.
        
        Takes in 4 arguments :
        
        s = double sided sample interferogram (with or without secondary fringes)
        
        fw = filter width, should be zero if interferogram is good quality.
        
        fmin = minimum filter height  (0<fmin<1)
        
        fmin = 0, eliminate all points within filter range.
        fmin = 1, full pass filter.
        
        highf = high frequency folding limit (parameter of the scan), can usually be found at 
        '/entry1/instrument/interferometer/opus_parameters/instrument_changed/high_folding_limit'
        
        dv = half of period of oscillatory fringes in absorbance spectrum/intensity spectrum, along the energy axis. 
           NB.  Needs to be in units of cm-1!
        
        When you input dv, the program will use that to position an inverse Blackman Harris filter to eliminate
        the oscillations in your spectrum. 
        
        Outputs:
        2D array
        
        0th array = intensity spectrum:
        Single channel spectrum is computed from the original interferogram, which was zerofilled to next factor of 2,
        apodized using a 3 term Blackman Harris function and phase corrected following the Mertz method.
        
        1st array= associated wavenumber axis.
        
        
        """
        
        
        single = s[:0.5*s.size] #in case of bifringent interferogram, take only one peak to analyse (avoids sinusoidal modulations)
        #zero filling(pad until 16,384 if array is below this number and up to 65536 points if array is larger)
        
        single0 = ft._zerofill(single,highf,zerofill,ymax,ymaxinterf,yscaling)
        
        angles = ft._mertz(single0)
        
        single0 = ft._apod(single0)
        
        filtered= ft2._filter(single0, fw, fmin,dv,highf)
        
        schannel = ft._fft(filtered,angles)
        
        #calculate the axis in frequency space
        #frequency axis
        lmda = 1/highf#cm (this is the laser waveleght you re using)
        k = np.arange(np.size(single0))
        v = np.divide(2*k,lmda*np.size(single0)) # f = k/(N*lambda) where k is range of values from zero to array size,
        self.dv = v[100]-v[99]
        kbig = np.arange(np.size(single0))
        vbig = np.divide(kbig,lmda*np.size(single0))
        
        schannel2 = ft2.singleChannel(s[0.5*s.size:], fw, fmin, highf, dv,zerofill,ymax,ymaxinterf,yscaling)
        final = np.add(schannel,schannel2[0])
        final = np.true_divide(final,2)
        self.single2 = single0

#        schannel = final/final.max() # normalise wrt highest peak

        return final,v
    def absorbance2(self, schannel, refer,highfold,zerofill,ymax,ymaxinterf,yscaling):
        """
        #############################
        User Function absorbance
        NB. WORKS ONLY IF A REFERENCE INTERFEROGRAM IS PROVIDED.
        
        inputs 4 arguments:
        
        schannel = single channel spectrum (no axis necessary, just the data array)
        refer = the reference, double sided interferogram (cannot calculate absorbance if reference interferogram has not been taken)
        highfold = high folding limit 
        zerofill = zerofill factor (eg. 1,2 or 4)
        
        returns:
        
        1D array containing the absorbance spectrum. 
        The energy axis is the same as the single channel axis. (in k numbers)
        
        """
        ft =fft()
        refer = ft.singleChannel2(refer,highfold,zerofill,ymax,ymaxinterf,yscaling)
        absorbance = -np.log10(schannel[0]/refer[0])
        return absorbance[0:absorbance.size/2]

