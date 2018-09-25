import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import astropy.units as u
from astropy import wcs
from astropy.coordinates import SkyCoord, SkyOffsetFrame, FK5

file_in='fits_files/test_file.fits'

create_dummy_file=True
#
# Central coordinate
ra0=15*(3+(30+ 25.6/60.)/60.)
dec0= (32+(10+ 25.6/60.)/60.)
# PA Angle in Degrees
PA_Angle=45.
# Inclination angle in Degrees
inclination=42.
# Distance in pc
distance=109.
# Stellar Mass in Solar Masses
Mstar = 4.2
#epsilon
epsilon=5e-6*u.deg
if create_dummy_file:
    # create dummy FITS file
    dx = 50
    data= np.zeros( (2*dx+1,2*dx+1))
    # hdu container
    hdu = fits.PrimaryHDU(data)
    pixel=(0.01*(u.arcsec).to(u.deg))
    #
    hdu.header['CTYPE1']='RA---TAN'
    hdu.header['CTYPE2']='DEC--TAN'
    hdu.header['CUNIT1']='deg'
    hdu.header['CUNIT2']='deg'
    hdu.header['CDELT1']=-1*pixel
    hdu.header['CDELT2']=pixel
    hdu.header['CRVAL1']=ra0
    hdu.header['CRVAL2']=dec0
    hdu.header['CRPIX1']=dx*1.0
    hdu.header['CRPIX2']=dx*1.0
    hdu.header['RADESYS']='FK5'
    hdu.header['EQUINOX']=2000.
    hdu.header['BUNIT']='km/s'
    hdu.header['BMAJ']=4*pixel
    hdu.header['BMIN']=4*pixel
    hdu.header['BPA']=0.0
    hdu.writeto( file_in, overwrite=True)


center = SkyCoord(ra0*u.deg, dec0*u.deg, frame='fk5')
# Load WCS 
header=fits.getheader(file_in)
w = wcs.WCS(header)
# Create xy array and then coordinates
x=np.arange(header['naxis1'])
y=np.arange(header['naxis2'])
xx, yy = np.meshgrid(x, y)
world = w.wcs_pix2world(xx.flatten(), yy.flatten(), 1)

radec = SkyCoord(world[0]*u.deg, world[1]*u.deg, frame='fk5')
radec_off = radec.transform_to(center.skyoffset_frame())

lon=radec_off[:].lon
lat=radec_off[:].lat
lon.shape=xx.shape
lat.shape=xx.shape

# plt.ion()
# plt.figure(figsize=(12,6))
# plt.subplot(1, 2, 1)
# plt.imshow(lon.value, origin='lowest')
# plt.title('$\Delta x$')
# #
# plt.subplot(1, 2, 2)
# plt.imshow(lat.value, origin='lowest')
# plt.title('$\Delta y$')
# This plots the offset in x and y
# This is all good now
angle=np.deg2rad(PA_Angle)
c, s = np.cos(angle), np.sin(angle)
# Rotate the axes
lon_PA = c*lon - s*lat
lat_PA = s*lon + c*lat
# Deprojection 
lat_PA /= np.cos(np.deg2rad(inclination))

dep_radius=np.clip(np.sqrt( lat_PA**2 + lon_PA**2), epsilon, None)
Kep_velo = 29.78 * np.sqrt( Mstar / (distance *(dep_radius.to(u.arcsec)).value)) * np.sin(np.deg2rad(inclination)) * np.sign(lon_PA)
#
# Plot
# plt.figure(figsize=(12,6))
# plt.subplot(1, 2, 1)
# plt.imshow(lon_PA.value, origin='lowest')
# plt.contour(lon_PA.value, c_levels=np.linspace(np.min(lon_PA.value),np.max(lon_PA.value),num=10), colors='k')
# plt.plot( [50, 50+20*np.cos(angle +np.deg2rad(90.))],[50,50+20*np.sin(angle+np.deg2rad(90.))], color='red')
# plt.title('$\Delta x_{PA}$')
# #
# plt.subplot(1, 2, 2)
# plt.imshow(lat_PA.value, origin='lowest')
# plt.contour(lat_PA.value, c_levels=np.linspace(np.min(lat_PA.value),np.max(lat_PA.value),num=10), colors='k')
# plt.plot( [50, 50+20*np.cos(angle +np.deg2rad(90.))],[50,50+20*np.sin(angle+np.deg2rad(90.))], color='red')
# plt.title('$\Delta y_{PA}$')
# Plot
plt.figure(figsize=(12,6))
plt.subplot(1, 2, 1)
plt.imshow(dep_radius.value, origin='lowest')
plt.contour(dep_radius.value, c_levels=np.linspace(np.min(dep_radius.value),np.max(dep_radius.value),num=10), colors='k')
plt.plot( [50, 50+20*np.cos(angle +np.deg2rad(90.))],[50,50+20*np.sin(angle+np.deg2rad(90.))], color='red')
plt.title('Deprojected radius, $r_{PA,i}$')
#
plt.subplot(1, 2, 2)
plt.imshow(Kep_velo, origin='lowest', cmap='RdYlBu_r')
plt.contour(Kep_velo, c_levels=np.linspace(np.min(Kep_velo),np.max(Kep_velo),num=10), colors='k')
plt.title('Projected $V_{Kep}$')
