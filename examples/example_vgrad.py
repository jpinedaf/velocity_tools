import matplotlib.pyplot as plt
plt.ion()

import numpy as np

import astropy.units as u
from astropy.io import fits

from velocity_tools import vgrad

# 11 x 11 image
# pixel of 1 arcsec
# delta x is negative, R.A.
x_t = -1* np.arange(111) * u.arcsec#.to(u.deg)
y_t = np.arange(111) * (u.arcsec)#.to(u.deg)
xv, yv = np.meshgrid(x_t, y_t, indexing='xy')
gx = -2.0 * u.km / u.s / u.deg
gy = 2.0 * u.km / u.s / u.deg
v0 = 8.5 * u.km / u.s
ev0 = 0.1 * u.km / u.s
v_t = (xv * gx + yv * gy).to(u.km / u.s) + v0 # km/s
ev_t = np.empty(v_t.shape) * u.km / u.s
ev_t[:] = ev0 # km/s

hdu = fits.PrimaryHDU(v_t.value)
hdu.header['CDELT1'] = (-1 * u.arcsec.to(u.deg), 'deg')
hdu.header['CDELT2'] = (1 * u.arcsec.to(u.deg), 'deg')
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
plt.close('all')
fig1, ax1 = plt.subplots(figsize=(8, 8), ncols=1)
pos1 = ax1.imshow(fits.getdata(file_vc), origin='lower', cmap='RdYlBu_r',
                interpolation='none')
ax1.text(0.1, 0.9, 'Vlsr', transform=ax1.transAxes)
fig1.colorbar(pos1, ax=ax1)

fig2, ax2 = plt.subplots(figsize=(8, 8), ncols=1)
pos2 = ax2.imshow(result['grad'], origin='lower', cmap='inferno', interpolation='none',
                    vmin=0.99 * input_grad.value, vmax=1.01 * input_grad.value)
ax2.text(0.1, 0.9, 'Gradient Vlsr', transform=ax2.transAxes)
fig2.colorbar(pos2, ax=ax2)

fig3, ax3 = plt.subplots(figsize=(8, 8), ncols=1)
pos3 = ax3.imshow(result['posang'], origin='lower', cmap='inferno', 
                interpolation='none', vmin=input_PA - 5, vmax=input_PA + 5)
ax3.text(0.1, 0.9, 'Gradient PA', transform=ax3.transAxes)
fig3.colorbar(pos3, ax=ax3)
