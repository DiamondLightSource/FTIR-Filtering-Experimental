'''
Created on 10 Jun 2015

@author: flb41892
'''
import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from PIL import Image
import io


f = h5py.File("/dls/science/groups/das/ExampleData/OpusData/Nexus/x15 Calu T.2.nxs","r")


dataset_energy = f["entry1/data/ratio_absorbance_energy"][...]




ds = f["entry1/data/ratio_absorbance_scan"][...]

#finding the average spectrum of the file, make sure we load all of the file into memory using "..." to optimise process


total = []

total = ds.mean(0).mean(0) #string together 2 operations to find average of first two indices (as fcn of third)

#plt.plot(total)

# total intensity of image
intensity = ds.sum(2)

#filter out 5% of data on each extremity 

vmin = np.amin(intensity)
vmax = np.amax(intensity)

test = np.sort(intensity.flatten())
vmin = test[np.round(test.size*0.01)]
vmax = test[test.size-np.round(test.size*0.01)-1]

region = 0.5*(vmax-vmin) # select 5% of max to select region to cutout
#above is bad way, does not get rid of data

f,(ax1,ax2) = plt.subplots(1,2)

ax1.plot(dataset_energy,total)
ax1.set_xlabel("Wavenumber [1/cm]")
ax1.set_ylabel("Absorption Spectrum")
ax1.set_title("Average Intensity Spectrum")
ax2.imshow(intensity,vmin = vmin,vmax = vmax)
ax2.set_xlabel("Horizontal Cell Index")
ax2.set_ylabel("Vertical Cell Index")
plt.close(f)
#remove additional whitespace, axis and labels
fig = plt.figure(frameon=False)
ax = plt.Axes(fig,[0.,0.,1.,1.])
ax.set_axis_off()
fig.add_axes(ax)
ax.imshow(intensity,vmin = vmin,vmax = vmax)

#set thumbnail size
size = 600,400

plt.show()
#convert what we see on the screen to bytes 
buf = io.BytesIO()
plt.savefig(buf,format ='png')
buf.seek(0)
im = Image.open(buf)
im.thumbnail(size, Image.ANTIALIAS)
im.save('/home/flb41892/data/ispyb/thumb5.png') # save thumbnail
plt.close()
