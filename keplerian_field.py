import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import astropy.units as u

# create dummy FITS file
dx = 50
data= np.zeros( (2*dx+1,2*dx+1))
# hdu container
hdu = fits.PrimaryHDU(data)
ra0=15*(3+(30+ 25.6/60.)/60.)
dec0= (32+(10+ 25.6/60.)/60.)
pixel=(1*(u.arcsec).to(u.deg))
#
hdu.header['CTYPE1']='RA---TAN'
hdu.header['CTYPE2']='DEC--TAN'
hdu.header['CUNIT1']='deg'
hdu.header['CUNIT2']='deg'
hdu.header['CDELT1']=-1*pixel
hdu.header['CDELT2']=pixel
hdu.header['CRVAL1']=ra0
hdu.header['CRVAL2']=dec0
hdu.header['CRPIX1']=dx
hdu.header['CRPIX2']=dx
hdu.header['RADESYS']='FK5'
hdu.header['EQUINOX']=2000.
hdu.header['BUNIT']='km/s'
hdu.header['BMAJ']=4*pixel
hdu.header['BMIN']=4*pixel
hdu.header['BPA']=0.0
