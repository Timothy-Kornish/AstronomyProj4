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


dir = ('/Volumes/TIM/Lab4')
dir1 = '/Volumes/TIM/Lab4/oct7/oct7_-'
npixelsx = 1024
npixelsy= 1024

cluster_r= np.zeros((1024,1024,25))
for i in range(0,25):
    cluster_r[:,:,i] = fits.getdata(dir+'/Data_r/Cluster_r_0'+str(i)+'.fits')

cluster_v= np.zeros((1024,1024,25))
for i in range(0,25):
    cluster_v[:,:,i] = fits.getdata(dir+'/Data_v/Cluster_v_0'+str(i)+'.fits')

cluster_b= np.zeros((1024,1024,25))
for i in range(0,25):
    cluster_b[:,:,i] = fits.getdata(dir+'/Data_b/Cluster_b_0'+str(i)+'.fits')

"""
#input and output directories, B-filter;
input_dir_B = dir + '/Data_r'
output_dir_B = dir + '/mosaic_B'
montage.mosaic(input_dir_B, output_dir_B, background_match=True, combine='median')

#R-filter
input_dir_R = dir + '/Data_r'
output_dir_R = dir + '/mosaic_R'
montage.mosaic(input_dir_R, output_dir_R, background_match=True, combine='median')

#V-filter
input_dir_V = dir + '/Data_v'
output_dir_V = dir + '/mosaic_V'
montage.mosaic(input_dir_V, output_dir_V, background_match=True, combine='median')
"""

#------------------
#Determine photometric zeropoint in each filter:
# stars are




ld1_b = np.zeros((1024,1024,5))
ld2_b = np.zeros((1024,1024,5))
ld3_b = np.zeros((1024,1024,5))

ld1_r = np.zeros((1024,1024,5))
ld2_r = np.zeros((1024,1024,5))
ld3_r = np.zeros((1024,1024,5))

ld1_v = np.zeros((1024,1024,5))
ld2_v = np.zeros((1024,1024,5))
ld3_v = np.zeros((1024,1024,5))

# Blue filter

for i in range(0,5):
    ld1_b[:,:,i] = fits.getdata(dir1+'00'+str(i)+'b.fits')

for i in range(0,4):
    ld2_b[:,:,i] = fits.getdata(dir1+'00'+str(i+6)+'b.fits')
ld2_b[:,:,5] = fits.getdata(dir1+'010b.fits')

for i in range(0,5):
    ld3_b[:,:,i] = fits.getdata(dir1+'0'+str(i+11)+'b.fits')

# Red filter

for i in range(0,5):
    ld1_r[:,:,i] = fits.getdata(dir1+'00'+str(i)+'r.fits')

for i in range(0,4):
    ld2_r[:,:,i] = fits.getdata(dir1+'00'+str(i+6)+'r.fits')
ld2_r[:,:,5] = fits.getdata(dir1+'010r.fits')

for i in range(0,5):
    ld3_r[:,:,i] = fits.getdata(dir1+'0'+str(i+11)+'r.fits')

# Vis filter

for i in range(0,5):
    ld1_v[:,:,i] = fits.getdata(dir1+'00'+str(i)+'v.fits')

for i in range(0,4):
    ld2_v[:,:,i] = fits.getdata(dir1+'00'+str(i+6)+'v.fits')
ld2_v[:,:,5] = fits.getdata(dir1+'010v.fits')

for i in range(0,5):
    ld3_v[:,:,i] = fits.getdata(dir1+'0'+str(i+11)+'v.fits')

############### mouse clicker ####################


def onclick(event):
    global xc, yc
    xc, yc = event.xdata, event.ydata
    print 'c = %d, %d'%(xc, yc)

    global coords
    coords.append((xc, yc))

    if len(coords) == 1:
        fig.canvas.mpl_disconnect(cid)
        plt.close(i)
    return

#for-loop to use the mouse click on Landolt images:

nlanimages = 5
nlanstars = 3
lanstarlocB = np.zeros((2,nlanimages,nlanstars))

#B-filter
L_image_B = np.zeros((npixelsx,npixelsy,nlanimages,nlanstars))
L_image_B[:,:,:,0] = ld1_b
L_image_B[:,:,:,1] = ld2_b
L_image_B[:,:,:,2] = ld3_b

for i in range(1,6):
    for j in range(1,4):
        fig = plt.figure(i)
        plt.gray()
        ax = fig.add_subplot(111)
        ax.imshow(L_image_B[:,:,i-1,j-1], vmin=0.25*np.median(L_image_B[:,:,i-1,j-1]),vmax=4*np.median(L_image_B[:,:,i-1,j-1]),interpolation='nearest')

        coords = []

        cid = fig.canvas.mpl_connect('button_press_event', onclick)
        plt.show(i)

        lanstarlocB[:,i-1,j-1] = xc, yc 

pickle.dump(lanstarlocB, open(dir + 'lanstarlocB.txt', 'wb')) 
lanstarloc_Btxt = pickle.load(open(dir + 'lanstarlocB.txt', 'rb'))



#R-filter
lanstarlocR = np.zeros((2,nlanimages,nlanstars))


L_image_R = np.zeros((npixelsx,npixelsy,nlanimages,nlanstars))
L_image_R[:,:,:,0] = ld1_r
L_image_R[:,:,:,1] = ld2_r
L_image_R[:,:,:,2] = ld3_r

for i in range(1,6):
    for j in range(1,4):
        fig = plt.figure(i)
        plt.gray()
        ax = fig.add_subplot(111)
        ax.imshow(L_image_R[:,:,i-1,j-1], vmin=0.25*np.median(L_image_R[:,:,i-1,j-1]),vmax=4*np.median(L_image_R[:,:,i-1,j-1]),interpolation='nearest')

        coords = []

        cid = fig.canvas.mpl_connect('button_press_event', onclick)
        plt.show(i)

        lanstarlocR[:,i-1,j-1] = xc, yc #storage array, saved below

pickle.dump(lanstarlocR, open(dir + 'lanstarlocR.txt', 'wb')) #done!
lanstarloc_Rtxt = pickle.load(open(dir + 'lanstarlocR.txt', 'rb'))




#V-filter
lanstarlocV = np.zeros((2,nlanimages,nlanstars))


L_image_V = np.zeros((npixelsx,npixelsy,nlanimages,nlanstars))
L_image_V[:,:,:,0] = ld1_v
L_image_V[:,:,:,1] = ld2_v
L_image_V[:,:,:,2] = ld3_v

for i in range(1,6):
    for j in range(1,4):
        fig = plt.figure(i)
        plt.gray()
        ax = fig.add_subplot(111)
        ax.imshow(L_image_V[:,:,i-1,j-1], vmin=0.25*np.median(L_image_V[:,:,i-1,j-1]),vmax=4*np.median(L_image_V[:,:,i-1,j-1]),interpolation='nearest')

        coords = []

        cid = fig.canvas.mpl_connect('button_press_event', onclick)
        plt.show(i)

        lanstarlocV[:,i-1,j-1] = xc, yc #storage array, saved below

pickle.dump(lanstarlocV, open(dir + 'lanstarlocV.txt', 'wb')) #done
lanstarloc_Vtxt = pickle.load(open(dir + 'lanstarlocV.txt', 'rb'))


########### 2D Gaussian #############

def twoD_Gaussian((x,y), amplitude, x0, y0, sigma_x, sigma_y, theta, offset):
    x0 = float(x0)
    y0 = float(y0)
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-x0)**2) + 2*b*(x-x0)*(y-y0) + c*((y-y0)**2)))
    return g.ravel()

#B filter

xsize = 30
ysize = 30 #for subarray around stars
x = np.arange(xsize)
y = np.arange(ysize)
x, y = np.meshgrid(x,y)

Lparam_arrayB = np.zeros((5,3,2))
Limage_coordsB = lanstarloc_Btxt
Linitial_guessB = (1000,xsize/2,ysize/2,5,5,0,30)
lancentersB = np.zeros((2,5,3))
lansubB = np.zeros((xsize,ysize,5,3))
lan_FWHM_B = np.zeros((5,3))
for i in range(0,5): #5 images
    for j in range(0,3): #3 stars
        lansubB[:,:,i,j] = L_image_B[Limage_coordsB[1,i,j]-15:Limage_coordsB[1,i,j]+15,Limage_coordsB[0,i,j]-15:Limage_coordsB[0,i,j]+15,i,j]
        LparamsB, LcovarB = opt.curve_fit(twoD_Gaussian, (x,y), lansubB[:,:,i,j].reshape(xsize*ysize), p0=Linitial_guessB)
        lan_FWHM_B[i,j] = math.fabs(2.3548*np.mean(LparamsB[3:5]))
        Lparam_arrayB[i,j,:] = LparamsB[1:3]
        lancentersB[:,i,j] = LparamsB[1:3] + lanstarloc_Btxt[:,i,j]-15

Ldata_fittedB = (twoD_Gaussian((x,y), *LparamsB)).reshape(ysize,xsize)


#R filter

Lparam_arrayR = np.zeros((5,3,2))
Limage_coordsR = lanstarloc_Rtxt
Linitial_guessR = (1000,xsize/2,ysize/2,5,5,0,30)
lancentersR = np.zeros((2,5,3))
lansubR = np.zeros((xsize,ysize,5,3))
lan_FWHM_R = np.zeros((5,3))
for i in range(0,5): #5 images
    for j in range(0,3): #3 stars
        lansubR[:,:,i,j] = L_image_R[Limage_coordsR[1,i,j]-15:Limage_coordsR[1,i,j]+15,Limage_coordsR[0,i,j]-15:Limage_coordsR[0,i,j]+15,i,j]
        LparamsR, LcovarR = opt.curve_fit(twoD_Gaussian, (x,y), lansubR[:,:,i,j].reshape(xsize*ysize), p0=Linitial_guessR)
        lan_FWHM_R[i,j] = math.fabs(2.3548*np.mean(LparamsR[3:5]))
        Lparam_arrayR[i,j,:] = LparamsR[1:3]
        lancentersR[:,i,j] = LparamsR[1:3] + lanstarloc_Rtxt[:,i,j]-15

Ldata_fittedR = (twoD_Gaussian((x,y), *LparamsR)).reshape(ysize,xsize)

#V filter

Lparam_arrayV = np.zeros((5,3,2))
Limage_coordsV = lanstarloc_Vtxt
Linitial_guessV = (500,xsize/2,ysize/2,5,5,0,30)
lancentersV = np.zeros((2,5,3))
lansubV = np.zeros((xsize,ysize,5,3))
lan_FWHM_V = np.zeros((5,3))
for i in range(0,5): #5 images
    for j in range(0,3): #3 stars
        lansubV[:,:,i,j] = L_image_V[Limage_coordsV[1,i,j]-15:Limage_coordsV[1,i,j]+15,Limage_coordsV[0,i,j]-15:Limage_coordsV[0,i,j]+15,i,j]
        LparamsV, LcovarV = opt.curve_fit(twoD_Gaussian, (x,y), lansubV[:,:,i,j].reshape(xsize*ysize), p0=Linitial_guessV)
        lan_FWHM_V[i,j] = math.fabs(2.3548*np.mean(LparamsV[3:5]))
        Lparam_arrayV[i,j,:] = LparamsV[1:3]
        lancentersV[:,i,j] = LparamsV[1:3] + lanstarloc_Vtxt[:,i,j]-15

Ldata_fittedV = (twoD_Gaussian((x,y), *LparamsV)).reshape(ysize,xsize)


############## Annulus and Aperature ###############
g= 1.39
RN = 8.4



#For the B-filter first
x_lan = np.zeros((1024,1024)) + np.arange(npixelsx)
y_lan = x_lan.T
nimages = 5
nstars = 3
B_distance = np.zeros((npixelsx,npixelsy,nimages,nstars))
B_annulus = np.zeros((npixelsx,npixelsy,nimages,nstars)).astype(bool)
B_annulus_inner = np.zeros((nimages,nstars))
B_annulus_outer = np.zeros((nimages,nstars))
B_center = np.zeros((2,nimages,nstars))
for i in range(0,5):
    for j in range(0,3):
        B_center = Limage_coordsB[:,i,j]
        B_annulus_inner[i,j] = 6*lan_FWHM_B[i,j]
        B_annulus_outer[i,j] = 8*lan_FWHM_B[i,j]
        B_distance[:,:,i,j] = ((x_lan-lancentersB[0,i,j])**2 + (y_lan-lancentersB[1,i,j])**2)**0.5
        B_annulus[:,:,i,j] = (B_annulus_inner[i,j]<B_distance[:,:,i,j]) & (B_distance[:,:,i,j]<B_annulus_outer[i,j])
        
skyB = np.zeros((nimages,nstars))
for i in range(0,5):
    for j in range(0,3):
        skyB[i,j] = np.median(L_image_B[:,:,i,j][B_annulus[:,:,i,j]])
skyB_med = np.zeros((3))
for i in range(0,3):
    skyB_med[i] = np.median(skyB[:,i])



N_B = np.zeros((18,5,3)) #NET star counts (sky-subtracted)
n_B_aperture = np.zeros((18,5,3)) #aperture pixel sum
for i in range(1,19):
    for j in range(0,5):
        for k in range(0,3):
            n_B_aperture[i-1,j,k] = (B_distance[:,:,j,k]<(i)).sum()
            N_B[i-1,j,k] = np.sum(L_image_B[:,:,j,k][(B_distance[:,:,j,k]<(i))])-(skyB_med[k]*n_B_aperture[i-1,j,k])


SNR_B = np.zeros((18,5,3))
for i in range(0,18): #0-18 pixels
    for j in range(0,5):
        for k in range(0,3):
            SNR_B[i,j,k] = g*N_B[i,j,k]/(((g*N_B[i,j,k])+(n_B_aperture[i,j,k]*g*skyB_med[k])+(n_B_aperture[i,j,k]*RN**2))**0.5)

plt.figure(1) #this isn't really useful for the actual data analysis, just visualization
plt.plot(np.arange(18)+1, SNR_B[:,0,0], 'bo')
plt.xlabel('Aperture Radius (pixels)',fontsize=15)
plt.ylabel('Signal to Noise Ratio',fontsize=15)


aperture_radius_B = np.zeros((5,3))
maxval_B = np.zeros((15))
for i in range(0,18):
    for j in range(0,5):
        for k in range(0,3):
            if SNR_B[i,j,k] == max(SNR_B[:,j,k]):
                aperture_radius_B[j,k] = i

best_aperture_B =  np.median(aperture_radius_B) # 5.0 (for the B filter)


#For the R filter

R_distance = np.zeros((npixelsx,npixelsy,nimages,nstars))
R_annulus = np.zeros((npixelsx,npixelsy,nimages,nstars)).astype(bool)
R_annulus_inner = np.zeros((nimages,nstars))
R_annulus_outer = np.zeros((nimages,nstars))
R_center = np.zeros((2,nimages,nstars))
for i in range(0,5):
    for j in range(0,3):
        R_center = Limage_coordsR[:,i,j]
        R_annulus_inner[i,j] = 6*lan_FWHM_R[i,j]
        R_annulus_outer[i,j] = 8*lan_FWHM_R[i,j]
        R_distance[:,:,i,j] = ((x_lan-lancentersR[0,i,j])**2 + (y_lan-lancentersR[1,i,j])**2)**0.5
        R_annulus[:,:,i,j] = (R_annulus_inner[i,j]<R_distance[:,:,i,j]) & (R_distance[:,:,i,j]<R_annulus_outer[i,j])
        
skyR = np.zeros((nimages,nstars))
for i in range(0,5):
    for j in range(0,3):
        skyR[i,j] = np.median(L_image_R[:,:,i,j][R_annulus[:,:,i,j]])
skyR_med = np.zeros((3))
for i in range(0,3):
    skyR_med[i] = np.median(skyR[:,i])



N_R = np.zeros((18,5,3)) #NET star counts (sky-subtracted)
n_R_aperture = np.zeros((18,5,3)) #aperture pixel sum
for i in range(1,19):
    for j in range(0,5):
        for k in range(0,3):
            n_R_aperture[i-1,j,k] = (R_distance[:,:,j,k]<(i)).sum()
            N_R[i-1,j,k] = np.sum(L_image_R[:,:,j,k][(R_distance[:,:,j,k]<(i))])-(skyR_med[k]*n_R_aperture[i-1,j,k])


SNR_R = np.zeros((18,5,3))
for i in range(0,18): #0-18 pixels
    for j in range(0,5):
        for k in range(0,3):
            SNR_R[i,j,k] = g*N_R[i,j,k]/(((g*N_R[i,j,k])+(n_R_aperture[i,j,k]*g*skyR_med[k])+(n_R_aperture[i,j,k]*RN**2))**0.5)


aperture_radius_R = np.zeros((5,3))
maxval_R = np.zeros((15))
for i in range(0,18):
    for j in range(0,5):
        for k in range(0,3):
            if SNR_R[i,j,k] == max(SNR_R[:,j,k]):
                aperture_radius_R[j,k] = i

best_aperture_R =  np.median(aperture_radius_R) # 7.0 (for the R filter)




#For the V filter

V_distance = np.zeros((npixelsx,npixelsy,nimages,nstars))
V_annulus = np.zeros((npixelsx,npixelsy,nimages,nstars)).astype(bool)
V_annulus_inner = np.zeros((nimages,nstars))
V_annulus_outer = np.zeros((nimages,nstars))
V_center = np.zeros((2,nimages,nstars))
for i in range(0,5):
    for j in range(0,3):
        V_center = Limage_coordsV[:,i,j]
        V_annulus_inner[i,j] = 6*lan_FWHM_V[i,j]
        V_annulus_outer[i,j] = 8*lan_FWHM_V[i,j]
        V_distance[:,:,i,j] = ((x_lan-lancentersV[0,i,j])**2 + (y_lan-lancentersV[1,i,j])**2)**0.5
        V_annulus[:,:,i,j] = (V_annulus_inner[i,j]<V_distance[:,:,i,j]) & (V_distance[:,:,i,j]<V_annulus_outer[i,j])
        
skyV = np.zeros((nimages,nstars))
for i in range(0,5):
    for j in range(0,3):
        skyV[i,j] = np.median(L_image_V[:,:,i,j][V_annulus[:,:,i,j]])
skyV_med = np.zeros((3))
for i in range(0,3):
    skyV_med[i] = np.median(skyV[:,i])



N_V = np.zeros((18,5,3)) #NET star counts (sky-subtracted)
n_V_aperture = np.zeros((18,5,3)) #aperture pixel sum
for i in range(1,19):
    for j in range(0,5):
        for k in range(0,3):
            n_V_aperture[i-1,j,k] = (V_distance[:,:,j,k]<(i)).sum()
            N_V[i-1,j,k] = np.sum(L_image_V[:,:,j,k][(V_distance[:,:,j,k]<(i))])-(skyV_med[k]*n_V_aperture[i-1,j,k])


SNR_V = np.zeros((18,5,3))
for i in range(0,18): #0-18 pixels
    for j in range(0,5):
        for k in range(0,3):
            SNR_V[i,j,k] = g*N_V[i,j,k]/(((g*N_V[i,j,k])+(n_V_aperture[i,j,k]*g*skyV_med[k])+(n_V_aperture[i,j,k]*RN**2))**0.5)


aperture_radius_V = np.zeros((5,3))
maxval_V = np.zeros((15))
for i in range(0,18):
    for j in range(0,5):
        for k in range(0,3):
            if SNR_V[i,j,k] == max(SNR_V[:,j,k]):
                aperture_radius_V[j,k] = i

best_aperture_V =  np.median(aperture_radius_V) # 6.0 (for the V filter)



############ Zeropoints of filter ##############


#B

exp_time_B = np.array([10,10,5]) #Landolts 1 & 2, 3 in the B filter
B_m_inst = np.zeros((5,3))
B_N_div_t = np.zeros((5,3))
for i in range(0,5):
    for j in range(0,3):
        B_N_div_t[i,j] = N_B[5,i,j]/exp_time_B[j] #5 is the best aperture radius
        B_m_inst[i,j] = -2.5*np.log10(B_N_div_t[i,j])

B_sim = np.array([9.174, 12.489, 12.341])
ZP_B = np.zeros((5,3))
for i in range(0,5):
    for j in range(0,3):
        ZP_B[i,j] = B_sim[j] - B_m_inst[i,j]
mean_ZP_B = np.mean(ZP_B, axis=0)
best_ZP_B = np.mean(mean_ZP_B)

SDOM_ZP_B = np.std(ZP_B, axis=0)/np.sqrt(5)
mean_SDOM_ZP_B = np.mean(SDOM_ZP_B)/np.sqrt(3)

#ZP in B filter is 21.3758854 +/- 0.0106

#R

exp_time_R = np.array([3,3,3]) #Landolts 1 & 2, 3 in the R filter
R_m_inst = np.zeros((5,3))
R_N_div_t = np.zeros((5,3))
for i in range(0,5):
    for j in range(0,3):
        R_N_div_t[i,j] = N_R[7,i,j]/exp_time_R[j] #7 is the best aperture radius
        R_m_inst[i,j] = -2.5*np.log10(R_N_div_t[i,j])

R_sim = np.array([8.844, 9.783, 9.205])
ZP_R = np.zeros((5,3))
for i in range(0,5):
    for j in range(0,3):
        ZP_R[i,j] = R_sim[j] - R_m_inst[i,j]
mean_ZP_R = np.mean(ZP_R, axis=0)
best_ZP_R = np.mean(mean_ZP_R)

SDOM_ZP_R = np.std(ZP_R, axis=0)/np.sqrt(5)
mean_SDOM_ZP_R = np.mean(SDOM_ZP_R)/np.sqrt(3)

#ZP in R filter is 21.77578466 +/- 0.0056

#V

exp_time_V = np.array([3,3,10]) #Landolts 1 & 2, 3 in the V filter
V_m_inst = np.zeros((5,3))
V_N_div_t = np.zeros((5,3))
for i in range(0,5):
    for j in range(0,3):
        V_N_div_t[i,j] = N_V[6,i,j]/exp_time_V[j] #6 is the best aperture radius
        V_m_inst[i,j] = -2.5*np.log10(V_N_div_t[i,j])

V_sim = np.array([8.965, 10.748, 10.382])
ZP_V = np.zeros((5,3))
for i in range(0,5):
    for j in range(0,3):
        ZP_V[i,j] = V_sim[j] - V_m_inst[i,j]
mean_ZP_V = np.mean(ZP_V, axis=0)
best_ZP_V = np.mean(mean_ZP_V)

SDOM_ZP_V = np.std(ZP_V, axis=0)/np.sqrt(5)
mean_SDOM_ZP_V = np.mean(SDOM_ZP_V)/np.sqrt(3)

#ZP in V filter is 21.68807 +/- 0.0058


#################

#################

##################  RBG IMAGING ########################

#################


aplpy.make_rgb_cube([dir + '/Data_r/mosaic_r/mosaic.fits', dir +  '/Data_v/mosaic_v/mosaic.fits', dir+  '/Data_b/mosaic_b/mosaic.fits'], dir + '/rgb_cube.fits')

aplpy.make_rgb_image(dir + '/rgb_cube.fits', dir + '/cluster_rgb.png', embed_avm_tags=False)

fig = aplpy.FITSFigure(dir+'/rgb_cube_2d.fits')
fig.show_rgb(dir+'/cluster_rgb.png')

fig.add_grid()
fig.grid.set_linestyle('dotted')




#############  HR Diagram ################


