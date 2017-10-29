import math
import numpy as np
from astropy.io import fits
import pylab as plt
import scipy.optimize as opt

dire = '/Volumes/TIM/Astro362/'

#####

#darks flats and fits

darks = np.zeros((1024,1024,9,4))
for i in range(0,9):
    for j in range(0,4):
        darks[:,:,i,j]= fits.getdata(dire+'Darks-00'+str(i+1)+'d'+str(j+1)+'.fits')
darkavgs = np.zeros((1024,1024,4))
for i in range(0,4):
    darkavgs[:,:,i] = np.mean(darks[:,:,:,i], axis =2)
    
flats = np.zeros((1024,1024,60))
for i in range(1,10):
    flats[:,:,i-1] = fits.getdata(dire+'Flats-00'+str(i)+'.fits')
for i in range(10,60):
    flats[:,:,i-1] = fits.getdata(dire+'Flats-0'+str(i)+'.fits')

flatmedians = np.zeros((1024,1024,3))
for i in range (0,3):
    flatmedians[:,:,i] = np.median(flats[:,:,(20*1):(20*(1+i))], axis =2)   
darkmedian = np.median(darks[:,:,:,0], axis=2)
mediansignals= np.zeros((3))
for i in range(0,3):
    mediansignals[i] = np.median(flatmedians[:,:,i]-darkmedian)
flats3 = np.zeros((1024,1024,3))
for i in range(0,3):
    flats3[:,:,i] = (flatmedians[:,:,i]-darkmedian)/(mediansignals[i])
normflat = np.median(flats3,axis=2)

##########

#landolt stars

#GSC4781575
ld1 = np.zeros((1024,1024,9))
for i in range(0,9):
    ld1[:,:,i] = (fits.getdata(dire+'GSC4781575-00'+str(i+1)+'.fits')-darkavgs[:,:,3])/normflat

#GSC447541
ld2 = np.zeros((1024,1024,9))
for i in range(0,9):
    ld2[:,:,i] = (fits.getdata(dire+'GSC447541-0'+str(i+10)+'.fits')-darkavgs[:,:,3])/normflat

#GSC5112062
ld3 = np.zeros((1024,1024,9))
for i in range(0,4):
    ld2[:,:,i] = (fits.getdata(dire+'GSC5112062-00'+str(i+1)+'.fits')-darkavgs[:,:,3])/normflat
ld2[:,:,4] = (fits.getdata(dire+'GSC5112062-006.fits')-darkavgs[:,:,3])/normflat
for i in range(5,7):
    ld2[:,:,i] = (fits.getdata(dire+'GSC5112062-00'+str(i+3)+'.fits')-darkavgs[:,:,3])/normflat
for i in range(7,9):
    ld2[:,:,i] = (fits.getdata(dire+'GSC5112062-0'+str(i+3)+'.fits')-darkavgs[:,:,3])/normflat


#GSC0047801575
ld4 = np.zeros((1024,1024,9))
for i in range(0,9):
    ld4[:,:,i] = (fits.getdata(dire+'GSC0047801575-00'+str(i+1)+'.fits')-darkavgs[:,:,3])/normflat

############

# Science Stars

#GSC10511778
sc1= np.zeros((1024,1024,9))
for i in range(0,9):
    sc1[:,:,i] = (fits.getdata(dire+'GSC10511778-00'+str(i+1)+'.fits')-darkavgs[:,:,1])/normflat

#GSC1084548
sc2= np.zeros((1024,1024,9))
for i in range(1,9):
    sc2[:,:,i-1] = (fits.getdata(dire+'GSC1084548-00'+str(i+1)+'.fits')-darkavgs[:,:,3])/normflat
    sc2[:,:,8] = (fits.getdata(dire+'GSC1084548-010.fits')-darkavgs[:,:,3])/normflat

#GSC31841316
sc3= np.zeros((1024,1024,9))
for i in range(0,7):
    sc3[:,:,i] = (fits.getdata(dire+'GSC31841316-00'+str(i+1)+'.fits')-darkavgs[:,:,1])/normflat
for i in range(7,9):
    sc3[:,:,i] = (fits.getdata(dire+'GSC31841316-0'+str(i+5)+'.fits')-darkavgs[:,:,1])/normflat

##########

# clicking 

def onclick(event):
    global xc, yc
    xc, yc = event.xdata, event.ydata
    print 'x = %d, y=%d'%(xc, yc)

    global coords
    coords.append((xc, yc))

    if len(coords) == 1:
        fig.canvas.mpl_disconnect(cid)
        plt.close(1)
    return

###########

#science star locations

sciences = np.zeros((1024,1024,9,3))
sciences[:,:,:,0] = sc1
sciences[:,:,:,1] = sc2
sciences[:,:,:,2] = sc3
    
sciences_starloc = np.zeros((9,2,3))
for i in range(0,9):
    for j in range(0,3):
        fig = plt.figure(1)
        plt.gray()
        ax = fig.add_subplot(111)
        ax.imshow(sciences[:,:,i:j],vmin=0.25*np.median(sciences[:,:,i:j]),vmax=4*np.median(sciences[:,:,i:j]),interpolation='nearest',origin = 'lower')

        coords = []

        # this command calls the function "onclick" that we defined earlier
        cid = fig.canvas.mpl_connect('button_press_event', onclick)
        plt.show(1)
        sciences_starloc[i,:,j] =xc,yc
############

# Landolt star locations

landolts =np.zeros((1024,1024,9,4))
landolts[:,:,:,0] = ld1
landolts[:,:,:,1] = ld2
landolts[:,:,:,2] = ld3  
landolts[:,:,:,3] = ld4

landolts_loc = np.zeros((9,2,4))
for i in range(0,9):
    for j in range(0,4):
        fig = plt.figure(1)
        plt.gray()
        ax = fig.add_subplot(111)
        ax.imshow(landolts[:,:,i:j],vmin=0.25*np.median(landolts[:,:,i:j]),vmax=4*np.median(landolts[:,:,i:j]),interpolation='nearest',origin = 'lower')

        coords = []

        # this command calls the function "onclick" that we defined earlier
        cid = fig.canvas.mpl_connect('button_press_event', onclick)
        plt.show(1)
        landolts_loc[i,:,j] = xc,yc
        
########### 
#2D - Gaussian

def twoD_Gaussian((x, y), amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    xo = float(xo)
    yo = float(yo)    
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                            + c*((y-yo)**2)))
    return g.ravel()
    
# parameters
initial_guess = (20000,15,15,5,5,0,30)
xsize = 30
ysize = 30
x = np.arange(xsize)
y = np.arange(ysize)
      
# subimages

science_sub = np.zeros((30,30,9,3))
for i in range(0,9):
    for j in range(0,3):
        science_sub[:,:,i,j] = sciences[sciences_starloc[i,1,j]-15:sciences_starloc[i,1,j]+15,sciences_starloc[i,0,j]-15:sciences_starloc[i,0,j]+15,i,j]    

science_sigma =np.zeros((9))
science_centers = np.zeros((2,9,3))
for i in range(0,9):
    for j in range(0,3):
        params, covar = opt.curve_fit(twoD_Gaussian, (x, y), science_sub[:,:,i:j].reshape(xsize*ysize), p0=initial_guess)
        science_centers[:,i:j]= params[1:3] +sciences_starloc[i,:,j]-15
        science_sigma[i] = (math.fabs(params[3])+math.fabs(params[4]))/2
        data_fitted = (twoD_Gaussian((x, y), *params))
        
        fig, ax = plt.subplots(1, 1)
        ax.hold(True)
        ax.imshow(science_sub[:,:,i,j], cmap=plt.cm.gray, origin='bottom',extent=(x.min(), x.max(), y.min(), y.max()))
        ax.contour(x, y, data_fitted.reshape(30,30), 8, colors='r')
        plt.xlabel('pixels')
        plt.ylabel('pixels')
        plt.show()
        
###########

# center for landolts

landolt_sub = np.zeros((30,30,9,4))
landolt_center = np.zeros((2,9,4))
for i in range (0,9):
    for j in range(0,4):
        landolt_sub[:,:,i,j] = landolts[landolts_loc[i,1,j]-15:landolts_loc[i,1,j]+15,landolts_loc[i,0,j]-15:landolts_loc[i,0,j]+15,i,j]
        
landolt_sigma= np.zeros((9))
landolt_centers= np.zeros((2,9,4))
        
for i in range(0,9):
    for j in range(0,3):
        params, covar = opt.curve_fit(twoD_Gaussian, (x, y), landolt_sub[:,:,i:j].reshape(xsize*ysize), p0=initial_guess)
        landolt_centers[:,i:j]= params[1:3] +landolts_loc[i,:,j]-15
        landolt_sigma[i] = (math.fabs(params[3])+math.fabs(params[4]))/2
        data_fitted = (twoD_Gaussian((x, y), *params))
        
        fig, ax = plt.subplots(1, 1)
        ax.hold(True)
        ax.imshow(landolt_sub[:,:,i,j], cmap=plt.cm.gray, origin='bottom',extent=(x.min(), x.max(), y.min(), y.max()))
        ax.contour(x, y, data_fitted.reshape(30,30), 8, colors='r')
        plt.xlabel('pixels')
        plt.ylabel('pixels')
        plt.show()        
        
        
###########

gain = 1.39
RN = 8.41

xarray = np.zeros((1024,1024))+np.arrange(1024)
yarray = xarray.T

distance = ((xarray-sciences_starloc[0,0,1])**2 +(yarray-sciences_starloc[0,1,1])**2)**0.5
N_star = np.zeros((20))
n_aperature = np.zeros((20))
for i in range(0,20):  
    n_aperature[i]=(distance<[i]).sum()
    N_star[i]=np.sum(sc2[:,:,0][distance<i])

FWHM = 2.3548*science_sigma[0]
inner_r=6*FWHM
outer_r = 8*FWHM
annulus = (distance < outer_r) & (distance > inner_r)
sky_counts= np.median(sc2[:,:,0][annulus])

snr = np.zeros((20))
for i in range(0,20):
    snr[i] = gain*N_star[i]/(((gain*N_star[i])+(n_aperature[i]*gain*sky_counts)+(n_aperature[i]*RN**2))**0.5)      
        
ranges = np.zeros((20))
plt.plot(ranges,snr,'go')
plt.plot(ranges,snr,'r--')
plt.xlabel('apparatus radius')
plt.ylabel('mean counts')
plt.show()
        
        
        
        
        