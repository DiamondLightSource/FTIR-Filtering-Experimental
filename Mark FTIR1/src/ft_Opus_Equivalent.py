import h5py
import numpy as np
import matplotlib.pyplot as plt
import cmath as m
from scipy import signal
import pylab
from myClasses import fft

# ft = fft()
# 
# 
# f = h5py.File("/home/flb41892/data/Nexus different parameters/NLC on Res2 ZF4 HFL7899 PR32.0.nxs","r")
# s = f["entry1/instrument/interferometer/sample_interferogram_scan"][...] #signal on which to perform FFT
# #ref = f["entry1/instrument/interferometer/reference_scan"][...] #noise signal
# highfold = f['/entry1/instrument/interferometer/opus_parameters/instrument/high_folding_limit'][...]
# 
# zerofill = f['/entry1/instrument/interferometer/opus_parameters/ft/zero_filling_factor'][...]
# zerofill =np.asarray(zerofill, int)
# refer = f['/entry1/instrument/interferometer/reference_scan'][...] #reference scan
# renergy = f['/entry1/instrument/interferometer/reference_energy'][...] #reference energy
# rint = f['/entry1/instrument/interferometer/reference_interferogram_scan'][...] #reference scan
# com = f["entry1/instrument/interferometer/sample_scan"][...]# this is the FT of the same file, as performed by opus
# 
# axis  = f["entry1/instrument/interferometer/sample_energy"][...] #signal on which to perform FFT
# 
# ymax = f['/entry1/instrument/interferometer/opus_parameters/sample_data_interferogram/y_maximum'][...]
# yscaling = f['/entry1/instrument/interferometer/opus_parameters/sample_data_interferogram/y_scaling_factor'][...]
# ymaxspect = f['/entry1/instrument/interferometer/opus_parameters/sample_data/y_maximum'][...]
#n = 13 #choose the index of the interferogram you want to analyse

ft = fft()
f = h5py.File("/home/flb41892/data/markmapping /4 scan 3 00um.0.hdf5","r") # load you nexus file here

s = f["Data/SampleInterferogram/yfolded"][...] #signal on which to perform FFT
s = s[20,20,:]
com = f["/Data/Sample/yfolded"][...]# this is the FT of the same file, as performed by opus
com = com[1,2,:]
highfold = f['/Parameters/Instrument/High Folding Limit'][...]

zerofill = f['/Parameters/FT/Zero Filling Factor'][...]
zerofill =np.asarray(zerofill, float)

refer = f['/Data/ReferenceInterferogram/yfolded'][...] #reference scan
refer = refer[40,30,:]

klaser = f['/Parameters/Instrument/Laser Wavenumber'][...]
#renergy =  f["entry1/instrument/interferometer/reference_energy"][...] # energy axis of reference scan 
#absenergy = f["entry1/instrument/interferometer/ratio_absorbance_energy"][...] # energy axis of reference scan 
ymax = f['/Parameters/SampleDataInterferogram/Y - Maximum'][...] #max amplitude of interferogram processed by Opus
yscaling = f['/Parameters/SampleDataInterferogram/Y - Scaling Factor'][...] #scaling factor that Opus applies to each intererigram before processing it.
ymaxspect = f['/Parameters/SampleData/Y - Maximum'][...]#scaling factor that Opus applies to the final spectrum before plotting it.
axis  = f["/Data/Sample/x"][...] #energy axis from Opus
s = refer
single = s[0:0.5*s.size] #in case of bifringent interferogram, take only one peak to analyse (avoids sinusoidal modulations)
#zero filling(pad until 16,384 if array is below this number and up to 65536 points if array is larger)
single = single - np.mean(single) # eliminate offset of interferogram
single = single*yscaling/(s.max()/ymax)
if highfold < 3950.0:
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
    if 2048<single.size < 4096:
            if zerofill < 4:
                single = np.concatenate((single,np.zeros(4096-single.size)))
            if zerofill == 4:
                single = np.concatenate((single,np.zeros(8192-single.size)))
    if 1024<single.size < 2048:
            if zerofill < 4:
                single = np.concatenate((single,np.zeros(2048-single.size)))
            if zerofill == 4:
                single = np.concatenate((single,np.zeros(4096-single.size)))
    if single.size < 1024:
            if zerofill < 4:
                single = np.concatenate((single,np.zeros(1024-single.size)))
            if zerofill == 4:
                single = np.concatenate((single,np.zeros(2048-single.size)))
    single = single*4
if 3950.0<highfold <7900.0:
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
    if 2048<single.size < 4096:
            if zerofill < 4:
                single = np.concatenate((single,np.zeros(4096-single.size)))
            if zerofill == 4:
                single = np.concatenate((single,np.zeros(8192-single.size)))
    if single.size < 2048:
            if zerofill < 4:
                single = np.concatenate((single,np.zeros(2048-single.size)))
            if zerofill == 4:
                single = np.concatenate((single,np.zeros(4096-single.size)))
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
    if 2048<single.size < 4096:
            if zerofill < 4:
                single = np.concatenate((single,np.zeros(4096-single.size)))
            if zerofill == 4:
                single = np.concatenate((single,np.zeros(8192-single.size)))
    if single.size < 2048:
            if zerofill < 4:
                single = np.concatenate((single,np.zeros(2048-single.size)))
            if zerofill == 4:
                single = np.concatenate((single,np.zeros(4096-single.size)))

#phase correction -Mertz method
n = 256 # number of points to select for phase correction about ZPD point
zeros = np.zeros(2*n)#make array of zeros of same length as the signal to be analysed
zeros[:] = single[np.argmax(single)-n:np.argmax(single)+n]
#ramp function (ramp is old triangular fcn, better to use the Black Harris 3 step fct w[t] )
ramp = np.zeros(2*n)
ramp[0:n] = np.linspace(0,1,n,endpoint=False)
ramp[n:] = np.linspace(1,0,n)
N = 2*n
w = np.zeros(N)
w2 = np.zeros(N)

for j in range(0,N):
    w[j] = 0.44959-.49364*np.cos(2*m.pi*j/N)+.05677*np.cos(4*m.pi*j/N)
    w2[j] = 0.42323-.49755*np.cos(2*m.pi*j/N)+.07922*np.cos(4*m.pi*j/N)
zeros = zeros*ramp #mpllultiply zeros array by ramp fcn to prepare array for phase correction
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

#frequency axis
lmda = 1./highfold#cm
k = np.arange(np.size(single))
v = np.divide(2.0*k,lmda*np.size(single)) # f = k/(N*lambda) where k is range of values from zero to array size,
kbig = np.arange(np.size(single))
vbig = np.divide(kbig,lmda*np.size(single))
#N is number of points in interferogram

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


pw = np.sqrt(np.add(np.square(freal) , np.square(fim)))

frealp = freal*pw
fimp = fim*pw

#apodization using a Black Harris 3 term fcn

apodf = np.zeros(single.size) #61dB
apodf2 = np.zeros(single.size) #67 dB
for j in range(0,single.size):
    apodf[j] = 0.44959-.49364*np.cos(2*m.pi*j/single.size)+.05677*np.cos(4*m.pi*j/single.size)
    apodf2[j] = 0.42323-.49755*np.cos(2*m.pi*j/single.size)+.07922*np.cos(4*m.pi*j/single.size)

ins = ((np.size(single)-np.argmax(single)) - (np.argmax(single)))/2
single = np.insert(single,0,np.zeros(ins))
single = single[:np.size(single)-ins]
single = single *apodf2
apod_singler = np.zeros(np.size(single))
apod_singler[0:single.size-np.argmax(single)] = single[np.argmax(single):single.size]
apod_singler[single.size-np.argmax(single):]=single[0:np.argmax(single)]

apod_singler2 = apod_singler
#apod_singler2[0:1500] =np.zeros(1500) 

#implement black harris inverse filter 

blh = np.ones(100)-np.blackman(100)

c = 100
m = 0.0
np.insert(blh,np.argmin(blh),np.multiply(m,np.ones(c)))
#apod_singler2[-1500:] =np.zeros(1500) 
exp = np.linspace(0,np.size(apod_singler2),np.size(apod_singler2))
exp[0:np.size(exp)] = np.exp(-exp/150)
f = exp[::-1]
exp[np.size(exp)/2:] = f[np.size(exp)/2:]
l = 300
exp[np.size(exp)/2-l:np.size(exp)/2+l] = np.ones(2*l)
#smoothen out middle of pass filter using gaussian fcn
d = signal.gaussian(200,33)
exp[np.argmin(exp):np.argmin(exp)+100]= d[0:np.size(d)/2]
exp[np.argmin(exp)-100:np.argmin(exp)] = d[np.size(d)/2:]

apod_singler2 = np.multiply(apod_singler2,exp)


    #can zerofill most on interferogram (useful to determine where the secondary peaks are)
    
    #output_axis1 = single.size
    #c = 300
    #apod_singler2[c:np.size(apod_singler2)-c] = np.zeros(np.size(apod_singler2)-2*c)
    
    #FFT the interferogram which was previously apodized and rotated
    


#extend phase arrays to match interferogram arrays(interpolation)
xp = np.arange(0,2*n)
x = np.arange(0,2*n,512./single.size)
cphi2 = np.interp(x,xp,cphi)
sphi2 = np.interp(x,xp,sphi)
#power spectrum


output_axis1 = np.size(apod_singler)
apodfi = np.fft.fft(apod_singler, output_axis1)
apodr = np.real(apodfi)
apodi = np.imag(apodfi)
#see difference between non eliminated and eliminated secondary fringes
#with fringe elim:
apodfi2 = np.fft.fft(apod_singler2, output_axis1)
apodr2 = np.real(apodfi2)
apodi2 = np.imag(apodfi2)
#multiply the complex fourier transform by the complex phase correction factor
finalr = np.multiply(apodr,cphi2) 
finali = np.multiply(apodi,sphi2)

finalr2 = np.multiply(apodr2,cphi2) 
finali2 = np.multiply(apodi2,sphi2)

final = np.add(finalr,finali)

final2 = np.add(finalr2,finali2)

#average two sided interferogram results
schannel = ft.singleChannel(s[0.5*s.size:], highfold,zerofill,ymax,s.max(),yscaling,ymaxspect)
#final = np.add(final,schannel[0])
#final = np.true_divide(final,2)

a = np.add(v,-axis[0])
a = np.abs(a)
a = a.argmin()
b = np.add(v,-axis[axis.size-1])
b = np.abs(b)
b = b.argmin()
full = final

final = final[a-1:b+1]  
final = np.multiply(final,com.max()/final.max())
refer = ft.singleChannel2(refer,highfold,zerofill,ymax,s.max(),yscaling,ymaxspect)
absorbance = -np.log10(schannel[0]/refer[0])

#normalisation wrt. highest peak


#plt.plot(axis,final*ymaxspect/final.max())
#plt.plot(axis,final)
#plt.plot(axis,com)
plt.plot(sphi,label ='sine of phase angle at position (40,30)')
pylab.legend(loc = 'upper right')

plt.show()

#comparison with Opus FFT

#dif = np.add(final_transf,-com)

#transmission calculation
#t=trans/reffr 
#absorption
#a = -np.log10(t)




#plot graphs 
#f,(ax1,ax2) = plt.subplots(2)


#ax1.plot(single)
#ax1.set_title("Initial Signal Slice")
#ax2.plot(frealp,'g',fimp,'r')

#ax2.set_title("FFT of Initial Signal Slice")

#pylab.plot(apodf2, label = 'Original Oscillation Plagued Spectrum')
# pylab.xlabel("cell index  ", fontsize = 17 )
# pylab.ylabel("BH3 Filter", fontsize = 17)

#pylab.legend(loc = 'upper right')

'''
Created on 17 Jun 2015

@author: flb41892
'''
