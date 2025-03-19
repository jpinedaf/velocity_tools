import numpy as np
import sys
import os

import astropy.units as u
from astropy.io import fits


def vfit(x, y, v, ev, distance=300.0*u.pc):
    """
    .. py:function:: vfit(x, y, v, ev, [distance=300.0*u.pc])

    It calculates the velocity gradient to a group of velocity 
    measuments. 

    :param float x: The array with the offsets in x-axis (usually delta RA) 
    in angle units.
    :param float y: The array with the offsets in y-axis (usually delta Dec) 
    in angle units.
    :param float v: The array with the velocity measurements (usually Vlsr) 
    in velocity units.
    :param float ev: The array with the uncertainty in the velocity 
    measurements (usually Vlsr) in velocity units.
    :float distance: The distance to the region studied with distance 
    units. Default value is 300 pc.
    :return: structure with results from the fit. The gradient is in units 
    of km/s/pc, while the position angle is in degrees measured 
    East-from-North. Velocities are reported in km/s.
    """
    # TODO:
    # -Check input arrays have same length
    #
    if (np.count_nonzero(np.isnan(x))):  sys.exit('There are NaNs in x-axis')
    if (np.count_nonzero(np.isnan(y))):  sys.exit('There are NaNs in y-axis')
    if (np.count_nonzero(np.isnan(v))):  sys.exit('There are NaNs in v-axis')
    if (np.count_nonzero(np.isnan(ev))):  sys.exit('There are NaNs in velocity error')
    if (np.count_nonzero(ev == 0.0)):  sys.exit('There are 0.0 in velocity error')
    #
    if (distance <= 0):  sys.exit('distance cannot be negative')

    x_d = (x.to('', equivalencies=u.dimensionless_angles()) * distance).to(u.pc).value
    y_d = (y.to('', equivalencies=u.dimensionless_angles()) * distance).to(u.pc).value

    # Obtain total weight, and average (x,y,v) to create new variables (dx,dy,dv)
    # which provide a lower uncertainty in the fit.
    #
    
    npts = x.shape
    
    v_d = v.to(u.km / u.s).value
    ev_d = ev.to(u.km / u.s).value
    wt = 1.0 / ev_d**2
    #
    suma = np.sum(wt)
    x_mean = np.sum(x_d * wt) / suma
    y_mean = np.sum(y_d * wt) / suma
    v_mean = np.sum(v_d * wt) / suma
    # relative units
    dx = (x_d - x_mean)
    dy = (y_d - y_mean)
    dv = (v_d - v_mean)
    #
    M = np.array([[np.sum(wt), np.sum(dx * wt), np.sum(dy * wt)], 
        [np.sum(dx * wt), np.sum(dx**2 * wt), np.sum(dx * dy * wt)], 
        [np.sum(dy * wt), np.sum(dx * dy * wt), np.sum(dy**2 * wt)]])
    try:
        covar = np.linalg.inv(M.T)
    except IOError:
        sys.exit('Singular matrix: no solution returned')
    coeffs = np.dot(covar, [[np.sum(dv * wt)], [np.sum(dx * dv * wt)],
                    [np.sum(dy * dv * wt)]])
    #
    errx = np.sqrt(covar[1, 1])
    erry = np.sqrt(covar[2, 2])
    #
    gx = coeffs[1][0]
    gy = coeffs[2][0]
    # best fit mean V_lsr and predicted velocities
    vc = coeffs[0][0] + v_mean
    vp = coeffs[0][0] + gx * dx + gy * dy
    grad = np.sqrt(gx**2 + gy**2)
    posang = np.rad2deg(np.arctan2(gx, gy))# + np.pi/2.)
    # chi^2_red
    red_chisq = np.sum((dv - vp)**2 * wt) / (len(dv) - 3.)

    vc_err = 0.
    # grad_err = np.sqrt((gx * errx)**2 + (gy * erry)**2) / grad
    # paerr = np.rad2deg(np.sqrt((gx / (gx**2 + gy**2))**2 * erry**2 +
    #                  (gy / (gx**2 + gy**2))**2 * errx**2))
    # Associated error including covariance
    grad_err = np.sqrt((gx * errx)**2 + (gy * erry)**2 
                     + 2 * gx * gy * covar[2, 1]) / grad
    
    paerr = np.rad2deg(np.sqrt((gx / (gx**2 + gy**2))**2 * erry**2 +    
                     (gy / (gx**2 + gy**2))**2 * errx**2 
                     - 2 * gx * gy / (gx**2 + gy**2)**2 * covar[2, 1]))
    chisq = red_chisq
    vp += v_mean

    results = {'grad': grad, 'grad_err': grad_err, 
               'posang': posang, 'paerr': paerr, 
               'chisq': chisq, 'vc': vc, 'vc_err': vc_err,
               'v_pred': vp}

    return results


def vfit_local(x, y, v, ev, distance=300.0*u.pc, width=45*u.arcsec, nmin=7.0):
    """
    .. py:function:: vfit_local(x, y, v, ev, [distance=300.0*u.pc],
    [width=45*u.arcsec], [nmin=7])

    It calculates the velocity gradient to a group of velocity 
    measuments. It uses a box of a requested width to determine which 
    pixels to use in the local calculation. The calculation is performed
    at every pixel, and it uses the vfit() function.



    nmin is the number of pixels, a value of 7 is a good value in the case
        of a Nyquist sampled dataset.
    :param float x: The array with the offsets in x-axis (usually delta RA) 
    in angle units.
    :param float y: The array with the offsets in y-axis (usually delta Dec) 
    in angle units.
    :param float v: The array with the velocity measurements (usually Vlsr) 
    in velocity units.
    :param float ev: The array with the uncertainty in the velocity 
    measurements (usually Vlsr) in velocity units.
    :float distance: The distance to the region studied with distance 
    units. Default value is 300 pc.
    :float width: The angular size of the box used in the pixel selection for 
    calculation, in units of angle. Default value is 45 arcsec.
    :float nmin: Minimum number of pixels in the requested box to perform the
    gradient fit. Default value is 7.
    :return: structure with results from all the fits.
    """
    
    grad = np.empty(x.shape)
    grad_err = np.empty(x.shape)
    posang = np.empty(x.shape)
    paerr = np.empty(x.shape)
    vc = np.empty(x.shape)
    chisq = np.empty(x.shape)
    #
    grad[:] = np.nan #
    grad_err[:] = np.nan #
    posang[:] = np.nan #
    paerr[:] = np.nan #
    vc[:] = np.nan #
    chisq[:] = np.nan

    for index, (x_i, y_i) in enumerate(zip(x, y)):
        # determine distances between all pixels to the one processing
        # and select all pixels within a radius of width.
        # dxy = np.sqrt((x - x_i)**2 + (y - y_i)**2)
        # gd = (dxy <= width)
        dx = np.abs(x - x_i)
        dy = np.abs(y - y_i)
        gd = (dx <= width * 0.5) & (dy <= width * 0.5)
        # determine the gradient only of there are enough pixels
        #
        if np.sum(gd) > nmin:
            res_grad = vfit(x[gd], y[gd], v[gd], ev[gd], distance=distance)
            grad[index] = res_grad['grad']
            grad_err[index] = res_grad['grad_err']
            posang[index] = res_grad['posang']
            paerr[index] = res_grad['paerr']
            vc[index] = res_grad['vc']
            chisq[index] = res_grad['chisq']

    # after all the calculations, we store the results and return them
    results = {'grad': grad, 'grad_err': grad_err, 
           'posang': posang, 'paerr': paerr, 
           'chisq': chisq, 'vc': vc, }
    return results


def vfit_image(file_vc, file_evc, distance=300.0*u.pc, width=None, 
        n_oversample=7):
    """
    .. py:function:: vfit_image(file_vc, file_evc, [distance=300.0*u.pc],
    [width=45*u.arcsec], [n_oversample=7])

    It loads the velocity and error in the velocity maps and calculates 
    the global velocity gradient and all the local velocity gradients 
    possible.
    It uses a box of a requested width to determine which 
    pixels to use in the local calculation.


    nmin is the number of pixels, a value of 7 is a good value in the case
        of a Nyquist sampled dataset.
    :param str file_vc: The filename to the FITS file with the centroid 
    velocity measurements. If no unit is present in the header, 
    it assumes km/s.
    :param str file_evc: The filename to the FITS file with the uncertainty 
    to the centroid velocity measurements. If no unit is present in the header, 
    it assumes km/s.
    :float distance: The distance to the region studied with distance 
    units. Default value is 300 pc.
    :float width: The angular size of the box used in the pixel selection for 
    calculation, in units of angle. Default value is 45 arcsec.
    :float n_oversample: Minimum number of pixels in the requested box to 
    perform the gradient fit (in a Nyquist sampled map). In the case of 
    oversampled maps (e.g., interferometers), the code takes the oversampling
    into account to correct for this. Default value is 7. 
    :return: structure with results from all the fits.
    """
    if (os.path.isfile(file_vc) == False):
        sys.exit('Velocity file not found')
    if (os.path.isfile(file_evc) == False):
        sys.exit('Uncertainty velocity file not found')
    
    vc, hd = fits.getdata(file_vc, header=True)
    evc = fits.getdata(file_evc)
    if (hd['BUNIT'] == 'm/s') | (hd['BUNIT'] == 'm s-1'):
        unit_velo = 1e-3 * (u.km / u.s)
    else:
        unit_velo = u.km / u.s
    x = np.arange(hd['NAXIS1']) * hd['CDELT1'] * u.deg
    y = np.arange(hd['NAXIS2']) * hd['CDELT2'] * u.deg
    vc *= unit_velo
    evc *= unit_velo
    xv, yv = np.meshgrid(x, y, indexing='xy')
    gd = (np.isfinite(vc * evc) & (evc.value != 0))
    
    if width == None:
        width = 2 * hd['BMAJ'] * u.deg

    beam_pix = np.abs((hd['BMAJ'] * hd['BMIN']) / 
                     (4 * hd['CDELT1'] * hd['CDELT2']))
    result_all = vfit(xv[gd], yv[gd], vc[gd], evc[gd], distance=distance)
    result = vfit_local(xv[gd], yv[gd], vc[gd], evc[gd], distance=distance, 
        width=width, nmin=n_oversample * beam_pix)
    grad = np.empty(vc.shape) + np.nan
    grad_err = np.empty(vc.shape) + np.nan
    grad_PA = np.empty(vc.shape) + np.nan
    grad_PA_err = np.empty(vc.shape) + np.nan
    grad_vc = np.empty(vc.shape) + np.nan
    # v_pred = grad

    grad[gd] = result['grad']
    grad_err[gd] = result['grad_err']
    grad_PA[gd] = result['posang']
    grad_PA_err[gd] = result['paerr']
    grad_vc[gd] = result['vc']
    #
    # Store global fit results in header
    # 
    print(result_all)
    hd['VG'] = result_all['grad']
    hd['VGerr'] = result_all['grad_err']
    hd['VG_PA'] = result_all['posang']
    hd['VG_PAerr'] = result_all['paerr']
    hd['VG_vc'] = result_all['vc']
    results_image = {'grad': grad, 'grad_err': grad_err, 
           'posang': grad_PA, 'paerr': grad_PA_err, 
           'vc': vc, 'header':hd}
    return results_image
