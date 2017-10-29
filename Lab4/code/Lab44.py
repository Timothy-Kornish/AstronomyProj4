import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import math
import montage_wrapper as montage
import pickle
import scipy.optimize as opt
import aplpy
from photutils import daofind
from photutils import aperture_photometry, CircularAperture
from photutils import CircularAnnulus
from astropy.table import hstack


dir1 = '/Volumes/TIM/20140902/darks-'
dir2 = '/Volumes/TIM/20140902/flat-'
dir3 = '/Volumes/TIM/20140902/Cluster1_Sept2-'


#-----------------
darks = np.zeros((1024,1024,10,2))
for i in range(0,10):
    for j in range(0,2):
        darks[:,:,i,j]= fits.getdata(dir1+'0'+str(i+66)+'d'+str((j+1)*5)+'.fit')
darkavgs= np.zeros((1024,1024,2))
for i in range(0,2):
    darkavgs[:,:,i] = np.mean(darks[:,:,:,i])
darkmedian= np.median(darks[:,:,:,0])
#-----------------

flats_b = np.zeros((1024,1024,18))
flats_v = np.zeros((1024,1024,18))
flats_r = np.zeros((1024,1024,18))

for i in range(0,9):
    flats_r[:,:,i] = fits.getdata(dir2+'00' + str(i+1)+'r.fit')
for i in range(9,18):
    flats_r[:,:,i] = fits.getdata(dir2+'0' + str(i+1)+'r.fit')

for i in range(0,9):
    flats_v[:,:,i] = fits.getdata(dir2+'00' + str(i+1)+'v.fit')
for i in range(9,18):
    flats_v[:,:,i] = fits.getdata(dir2+'0' + str(i+1)+'v.fit')

for i in range(0,9):
    flats_b[:,:,i] = fits.getdata(dir2+'00' + str(i+1)+'b.fit')
for i in range(9,18):
    flats_b[:,:,i] = fits.getdata(dir2+'0' + str(i+1)+'b.fit')

flats_b_med = np.zeros((1024,1024,18))
flats_v_med = np.zeros((1024,1024,18))
flats_r_med = np.zeros((1024,1024,18))

for i in range(0,18):
    flats_b_med[:,:,i]=np.median(flats_b[:,:,(1):((1+i))])
    flats_v_med[:,:,i]=np.median(flats_v[:,:,(1):((1+i))])
    flats_r_med[:,:,i]=np.median(flats_r[:,:,(1):((1+i))])

flats= np.zeros((1024,1024,18,3))
flats[:,:,:,0] = flats_r_med
flats[:,:,:,1] = flats_v_med
flats[:,:,:,2] = flats_b_med

mediansignal = np.zeros((18,3))
for i in range(0,18):
    for j in range(0,3):
        mediansignal[i,j] = np.median(flats[:,:,i,j]-darkmedian)

flats3 = np.zeros((1024,1024,18,3))
for j in range(0,18):
    for i in range(0,3):
        flats3[:,:,j,i] = (flats[:,:,j,i]-darkmedian)/(mediansignal[j,i])
normflat = np.median(flats3)


Cluster_r = np.zeros((1024,1024,25))
Cluster_v = np.zeros((1024,1024,25))
Cluster_b = np.zeros((1024,1024,25))

for i in range(0,10):
    Cluster_r[:,:,i], hd1 = (fits.getdata(dir3+ '0'+str(i+16)+'r.fit', header=True)-darkavgs[:,:,0])/normflat
    Cluster_v[:,:,i], hd2 = (fits.getdata(dir3+ '0'+str(i+16)+'v.fit', header=True)-darkavgs[:,:,0])/normflat
    Cluster_b[:,:,i], hd3 = (fits.getdata(dir3+ '0'+str(i+16)+'b.fit', header=True)-darkavgs[:,:,1])/normflat
    
for i in range(10,25):
    Cluster_r[:,:,i], hd1 = (fits.getdata(dir3+ '0'+str(i+16)+'r.fit', header=True)-darkavgs[:,:,0])/normflat
    Cluster_v[:,:,i], hd2 = (fits.getdata(dir3+ '0'+str(i+16)+'v.fit', header=True)-darkavgs[:,:,0])/normflat
    Cluster_b[:,:,i], hd3 = (fits.getdata(dir3+ '0'+str(i+16)+'b.fit', header=True)-darkavgs[:,:,1])/normflat

fits.writeto('Cluster_r',Cluster_r, header=hd1)
fits.writeto('Cluster_v',Cluster_v, header=hd2)
fits.writeto('Cluster_b',Cluster_b, header=hd3)




"""
#input and output directories, B-filter;
input_dir_B = dir + '/clust_B_new'
output_dir_B = dir + '/mosaic_B'
montage.mosaic(input_dir_B, output_dir_B, background_match=True, combine='median')

#R-filter
input_dir_R = dir + '/clust_R_new'
output_dir_R = dir + '/mosaic_R'
montage.mosaic(input_dir_R, output_dir_R, background_match=True, combine='median')

#V-filter
input_dir_V = dir + '/clust_V_new'
output_dir_V = dir + '/mosaic_V'
montage.mosaic(input_dir_V, output_dir_V, background_match=True, combine='median')
"""