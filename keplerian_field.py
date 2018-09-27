import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import astropy.units as u
from astropy import wcs
from astropy.coordinates import SkyCoord, SkyOffsetFrame, FK5

file_in='fits_files/test_file.fits'
file_out='fits_files/test_Vc.fits'
create_dummy_file=True
save_Vfield=True
do_plot=True

#
# Central coordinate
ra0=15*(3+(30+ 25.6/60.)/60.)
dec0= (32+(10+ 25.6/60.)/60.)
# PA Angle in Degrees
PA_Angle=142.*u.deg
# Inclination angle in Degrees
inclination=42.*u.deg
# Distance in pc
distance=110.02*u.pc
# Outer radius in au
R_out=300.*u.au
# Stellar Mass in Solar Masses
Mstar = 2.2*u.Msun
#epsilon in au
epsilon=10*u.au
# V_lsr
Vc= 5.2*u.km/u.s

if create_dummy_file:
    # create dummy FITS file
    dx = 1500
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
# This plots the offset in x and y
# This is all good now
c, s = np.cos(PA_Angle), np.sin(PA_Angle)
# Rotate the axes
# lon_PA = c*lon - s*lat
# lat_PA = s*lon + c*lat
lat_PA = c*lat - s*lon
lon_PA = s*lat + c*lon
# Deprojection 
lat_PA /= np.cos(inclination)
dep_angle=np.sqrt( lat_PA**2 + lon_PA**2)
angle_PA = np.arctan2(lat_PA, lon_PA)
#
dep_radius=np.clip( dep_angle.to('', equivalencies=u.dimensionless_angles())*distance.to(u.au), epsilon.to(u.au), None)
Kep_velo = 29.78 * np.sqrt( Mstar/u.Msun / (dep_radius/u.au)) * np.sin(inclination) * np.cos(angle_PA)*u.km/u.s
Kep_velo += Vc
Kep_velo[dep_radius > R_out] = np.nan
#
if do_plot:
    # Plot
    axes_extent=[ (lon.to(u.arcsec).value).max(), (lon.to(u.arcsec).value).min(),
                  (lat.to(u.arcsec).value).min(), (lat.to(u.arcsec).value).max()]
    plt.ion()
    plt.figure(figsize=(12,6))
    plt.subplot(1, 2, 1)
    plt.imshow(dep_radius.value, origin='lowest', interpolation='none', extent=axes_extent)
    x0, y0= w.wcs_world2pix( center.ra, center.dec, 1)
    dx_line= (lat.to(u.arcsec).value).max() * 0.2
    plt.contour(dep_radius.value, c_levels=np.linspace(np.min(dep_radius.value),np.max(dep_radius.value),num=10), colors='k', extent=axes_extent)
    plt.plot( [0.0, dx_line*np.sin(PA_Angle)],[0.0,dx_line*np.cos(PA_Angle)], color='gray')
    plt.plot( [0.0, dx_line*np.sin(PA_Angle+90.*u.deg)],[0.0,dx_line*np.cos(PA_Angle+90.*u.deg)], color='red')
    plt.title('Deprojected radius, $r_{PA,i}$')
    #
    plt.subplot(1, 2, 2)
    plt.imshow(Kep_velo.value, origin='lowest', cmap='RdYlBu_r', interpolation='none', extent=axes_extent)
    plt.contour(Kep_velo.value, c_levels=np.linspace(np.nanmin(Kep_velo.value),np.nanmax(Kep_velo.value),num=10), colors='k', extent=axes_extent)
    plt.title('Projected $V_{Kep}$')

if save_Vfield:
    header2=header
    header2['BUNIT']='km/s'
    fits.writeto( file_out, Kep_velo.value, header2)