import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

dir = '/Volumes/TIM/Lab4/oct7/oct7_-'

image = np.zeros((1024,1024,16))
image2 = np.zeros((1024,1024,16))
image3 = np.zeros((1024,1024,16))

for i in range(0,10):
    image[:,:,i] = fits.getdata(dir+'00'+str(i)+'b.fits')
for i in range(10,16):
    image[:,:,i] = fits.getdata(dir+'00'+str(i)+'b.fits')
    
for i in range(0,10):
    image2[:,:,i] = fits.getdata(dir+'00'+str(i)+'r.fits')
for i in range(10,16):
    image2[:,:,i] = fits.getdata(dir+'00'+str(i)+'r.fits')

for i in range(0,10):
    image3[:,:,i] = fits.getdata(dir+'00'+str(i)+'v.fits')
for i in range(10,16):
    image3[:,:,i] = fits.getdata(dir+'00'+str(i)+'v.fits')

plt.figure(1)
plt.imshow(image
plt.figure(2)

plt.figure(3)
