import numpy as np
from astropy.io import fits



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

flat = np.zeros((1024,1024,18,3))
for j in range(0,18):
    for i in range(0,3):
        flat[:,:,j,i] = (flats[:,:,j,i]-darkmedian)/(mediansignal[j,i])
normflat = np.median(flat)


Cluster_r = np.zeros((1024,1024,25))
Cluster_v = np.zeros((1024,1024,25))
Cluster_b = np.zeros((1024,1024,25))

for i in range(0,10):
    Cluster_r[:,:,i], hd1 = (fits.getdata(dir3+ '0'+str(i+16)+'r.fit', header=True))
    Cluster_v[:,:,i], hd2 = (fits.getdata(dir3+ '0'+str(i+16)+'v.fit', header=True))
    Cluster_b[:,:,i], hd3 = (fits.getdata(dir3+ '0'+str(i+16)+'b.fit', header=True))
    
    Cluster_r[:,:,i] = (Cluster_r[:,:,i] -darkavgs[:,:,0])/normflat
    Cluster_v[:,:,i] = (Cluster_v[:,:,i] -darkavgs[:,:,0])/normflat
    Cluster_b[:,:,i] = (Cluster_b[:,:,i] -darkavgs[:,:,1])/normflat
    
    fits.writeto('/Volumes/TIM/Lab4/Data_r/Cluster_r_0'+str(i)+'.fits',Cluster_r[:,:,i], header=hd1)
    fits.writeto('/Volumes/TIM/Lab4/Data_v/Cluster_v_0'+str(i)+'.fits',Cluster_v[:,:,i], header=hd2)
    fits.writeto('/Volumes/TIM/Lab4/Data_b/Cluster_b_0'+str(i)+'.fits',Cluster_b[:,:,i], header=hd3)

for i in range(10,25):
    Cluster_r[:,:,i], hd1 = (fits.getdata(dir3+ '0'+str(i+16)+'r.fit', header=True))
    Cluster_v[:,:,i], hd2 = (fits.getdata(dir3+ '0'+str(i+16)+'v.fit', header=True))
    Cluster_b[:,:,i], hd3 = (fits.getdata(dir3+ '0'+str(i+16)+'b.fit', header=True))
    
    Cluster_r[:,:,i] = (Cluster_r[:,:,i] -darkavgs[:,:,0])/normflat
    Cluster_v[:,:,i] = (Cluster_v[:,:,i] -darkavgs[:,:,0])/normflat
    Cluster_b[:,:,i] = (Cluster_b[:,:,i] -darkavgs[:,:,1])/normflat
    
    fits.writeto('/Volumes/TIM/Lab4/Data_r/Cluster_r_0'+str(i)+'.fits',Cluster_r[:,:,i], header=hd1)
    fits.writeto('/Volumes/TIM/Lab4/Data_v/Cluster_v_0'+str(i)+'.fits',Cluster_v[:,:,i], header=hd2)
    fits.writeto('/Volumes/TIM/Lab4/Data_b/Cluster_b_0'+str(i)+'.fits',Cluster_b[:,:,i], header=hd3)



