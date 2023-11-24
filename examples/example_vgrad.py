import matplotlib.pyplot as plt

from astropy.wcs import WCS
import astropy.units as u
from astropy.io import fits
from astropy.visualization.wcsaxes import add_beam, add_scalebar

plt.ion()

import numpy as np


from velocity_tools import vgrad

plt.rcParams.update({"text.usetex": True,
                     "font.family": "serif",
                     'xtick.direction': 'in',
                     'ytick.direction': 'in',
                     'xtick.major.size': 6,
                     'ytick.major.size': 6,
                     'xtick.minor.size': 3,
                     'ytick.minor.size': 3})

def plot_setup_wcsaxes(ax, header, label_col='black'):
    scalebar_length = 100*u.au
    scalebar_text="100 au"
    tick_minor_x = 4
    tick_minor_y = 6
    # ticks
    ra_ax = ax.coords[0] # lon
    dec_ax = ax.coords[1] # lat
    ra_ax.set_major_formatter('hh:mm:ss')
    dec_ax.set_major_formatter('dd:mm:ss')
    ra_ax.display_minor_ticks(True)
    dec_ax.display_minor_ticks(True)
    ra_ax.set_minor_frequency(tick_minor_x)
    dec_ax.set_minor_frequency(tick_minor_y)
    ax.tick_params(which='major', length=6)
    ax.tick_params(which='minor', length=3)
    ax.autoscale(enable=False)
    # Add beamsize
    add_beam(ax, header=header, frame=False, pad=0.5, 
             color=label_col, corner='top right')
    # # Scalebar
    length = (scalebar_length / my_distance).to(
        u.deg, u.dimensionless_angles()
        )
    add_scalebar(ax, length, label=scalebar_text, color=label_col, 
        corner='bottom right')
    #
    ax.set_xlabel(r'Right Ascension (J2000)')
    ax.set_ylabel(r'Declination (J2000)')
    return #fig_i



my_figsize = (6.5, 5)
my_distance = 10.0*u.pc

# 111 x 111 image
# pixel of 1 arcsec
# delta x is negative, R.A.
n_pix = 111
x_t = -1* np.arange(n_pix) * u.arcsec#.to(u.deg)
y_t = np.arange(n_pix) * (u.arcsec)#.to(u.deg)
xv, yv = np.meshgrid(x_t, y_t, indexing='xy')
gx = -2.0 * u.km / u.s / u.deg
gy = -2.0 * u.km / u.s / u.deg
v0 = 8.5 * u.km / u.s
ev0 = 0.1 * u.km / u.s
v_t = (xv * gx + yv * gy).to(u.km / u.s) + v0 # km/s
ev_t = np.empty(v_t.shape) * u.km / u.s
ev_t[:] = ev0 # km/s

hdu = fits.PrimaryHDU(v_t.value)
#
hdu.header['CDELT1'] = (-1 * u.arcsec.to(u.deg), 'deg')
hdu.header['CDELT2'] = (1 * u.arcsec.to(u.deg), 'deg')
hdu.header['CRPIX1'] = (len(x_t.value) * 0.5, 'Pixel reference')
hdu.header['CRPIX2'] = (len(y_t.value) * 0.5, 'Pixel reference')
hdu.header['CRVAL1'] = (15*3.5, 'deg')
hdu.header['CRVAL2'] = (31 + 1./3., 'deg')
hdu.header['CTYPE1'] = 'RA---TAN'
hdu.header['CTYPE2'] = 'DEC--TAN'
#
hdu.header['BUNIT'] = ('km / s', 'Vlsr')
hdu.header['BMAJ'] = 2 * hdu.header['CDELT2']
hdu.header['BMIN'] = hdu.header['BMAJ']
hdu.header['BPA'] = 0.0

file_vc = 'test_vel.fits'
file_evc = 'test_vel_err.fits'
hdu.writeto(file_vc, overwrite=True)
hdu.data = ev_t.value
hdu.writeto(file_evc, overwrite=True)

result_all = vgrad.vfit(xv, yv, v_t, ev_t, distance=10.0*u.pc)
# 1 parsec = 1 * u.pc.to(u.au) 
# = (1 * u.pc.to(u.au) / 10. * u.arcsec).to(u.deg)
conv_grad = (1 * u.pc.to(u.au) / 10. * u.arcsec).to(u.deg) / u.pc
input_grad = conv_grad * np.sqrt(gx**2 + gy**2)
input_PA = np.rad2deg(np.arctan2(gx.value, gy.value))

print('VGrad input = {0:4.2f}    VGrad output = {1:4.2f}'.format(input_grad, 
                        result_all['grad']))

result = vgrad.vfit_image(file_vc, file_evc, distance=10.0*u.pc, 
                          n_oversample=7, width=3*u.arcsec)
print('VGrad input = {0:4.2f}    VGrad output = {1:4.2f}'.format(input_grad, 
                        result['header']['VG']))

print('Input PA = {0:5.2f}'.format(input_PA))
#
#

# coordinate system for all images
wcs = WCS(hdu.header)

plt.close('all')
#
# Velocity map
#
fig1 = plt.figure(figsize=my_figsize)
ax1 = plt.subplot(projection=wcs)
im1 = ax1.imshow(fits.getdata(file_vc), cmap='RdYlBu_r', interpolation='none')
plot_setup_wcsaxes(ax1, hdu.header, label_col='black')
ax1.text(0.1, 0.9, 'Vlsr', transform=ax1.transAxes)
fig1.colorbar(im1, ax=ax1)


# This requires to fix display for bad/edge data points
cmap = plt.get_cmap(name='inferno').copy()
cmap.set_bad(color='0.85')

#
# Gradient amplitude map
#
fig2 = plt.figure(figsize=my_figsize)
ax2 = plt.subplot(projection=wcs)
im2 = ax2.imshow(result['grad'], cmap=cmap, interpolation='none',
    vmin=0.99 * input_grad.value, vmax=1.01 * input_grad.value)
plot_setup_wcsaxes(ax2, hdu.header, label_col='black')
ax2.text(0.1, 0.9, 'Gradient Vlsr', transform=ax2.transAxes)
fig2.colorbar(im2, ax=ax2)


fig3 = plt.figure(figsize=my_figsize)
ax3 = plt.subplot(projection=wcs)
im3 = ax3.imshow(result['posang'], cmap=cmap, interpolation='none',
    vmin=input_PA - 5, vmax=input_PA + 5)
plot_setup_wcsaxes(ax2, hdu.header, label_col='black')
ax3.text(0.1, 0.9, 'Gradient PA', transform=ax3.transAxes)
fig3.colorbar(im3, ax=ax3)

#
# Add quivers to velocity map, while undersampling the gradient results
ns = 10
x_map = np.arange(0, hdu.header['NAXIS1'])
y_map = np.arange(0, hdu.header['NAXIS2'])
xv_map, yv_map = np.meshgrid(x_map, y_map)
ra_map, dec_map = wcs.wcs_pix2world(xv_map, yv_map, 0)
#
gd_grad = np.isfinite(result['grad']) & (xv_map % ns == 0) & (yv_map % ns == 0)
ra_grad = ra_map[gd_grad]
dec_grad = dec_map[gd_grad]
grad_abs = result['grad'][gd_grad]
grad_pa = result['posang'][gd_grad]
#
grad_x = grad_abs * np.cos(np.deg2rad(grad_pa + 90))
grad_y = grad_abs * np.sin(np.deg2rad(grad_pa + 90))
"""
I don't know how to plot the quivers in sky unit coordinates
ax1.quiver(ra_grad*u.deg, dec_grad*u.deg, grad_x*u.deg, grad_y*u.deg, scale=10, 
    transform=ax1.get_transform('world'))
"""
quiv = ax1.quiver(xv_map[gd_grad], yv_map[gd_grad], grad_x, grad_y, scale=5, 
    pivot='mid', width=0.003, scale_units='xy', angles='xy', 
    headlength=3.5, headaxislength=3.0)
ax1.quiverkey(quiv, X=0.5*hdu.header['NAXIS1'], Y=1.03*hdu.header['NAXIS2'], 
    U=10, label=r'10 km s$^{-1}$ pc$^{-1}$', coordinates='data')
