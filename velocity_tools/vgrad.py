import numpy as np
import sys

import astropy.units as u
from astropy.io import fits

def vfit(x, y, v, ev, distance=300.0*u.pc):
"""
    It calculates the velocity gradient to a group of velocity 
    measuments.

    INPUTS:
        x:          array of RA off-set positions with units.
        y:          array of Dec off-set positions with units.
        v:          Velocity at the position with units.
        ev:         Uncertainty in velocity with units.
    KEYWORD PARAMETERS:
        Width:      FWHM of the weight function, in degrees.
        distance:   Distance to the region in pc.
    OUTPUTS:
        Vc:         Mean centroid velocity in km/s
        Vc_err:     Uncertainty of the mean centroid velocity in km/s
        Grad:       Velocity gradient in km/s/pc
        Grad_Err:   Uncertainty associated to the velocity gradient (km/s/pc)
        PosAng:     Position angle of the fitted gradient, in degrees
        PAErr:      Uncertainty of the position angle (degrees)
        ChiSq:      Chi^2
"""
    # TODO:
    # -Check for no NaNs in the input arrays
    # -Check input arrays have same length
    # -Check distance is positive
    #

    npts = x.shape
    wt = 1.0 / ev**2
    # dtor = np.pi/180.*1000
    if (distance >= 0):
        # conv = distance * dtor #X and Y are in degrees => pc
        # comb = [resolution.to(u.rad).value, 
        #         x.to(u.rad).value, y.to(u.rad).value] * u.rad
        # comb_d = comb.to('', equivalencies=u.dimensionless_angles()) * distance
        # resolution_i, x_i, y_i = comb_d
        x_d = x.to('', equivalencies=u.dimensionless_angles()) * distance
        y_d = y.to('', equivalencies=u.dimensionless_angles()) * distance
        # resolution *= conv
        # X = X * conv
        # Y = Y * conv
    else:
        sys.exit('distance is negative!!')
        return

    # Obtain total weight, and average (x,y,v) to create new variables (dx,dy,dv)
    # which provide a lower uncertainty in the fit.
    #
    suma = np.sum(wt)
    x_mean = np.sum(x_d * wt) / suma
    y_mean = np.sum(y_d * wt) / suma
    v_mean = np.sum(v * wt) / suma
    # relative units
    dx = (x_d - x_mean)
    dy = (y_d - y_mean)
    dv = (v - v_mean)
    #
    M = [[np.sum(wt), np.sum(dx * wt), np.sum(dy * wt)], 
        [np.sum(dx * wt), np.sum(dx**2 * wt), np.sum(dx * dy * wt)], 
        [np.sum(dy * wt), np.sum(dx * dy * wt), np.sum(dy**2 * wt)]]
    #print M
    try:
        covar = np.linalg.inv(M)
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
    vc = coeffs[0] + v_mean
    vp = coeffs[0] + coeffs[1] * dx + coeffs[2] * dy
    grad = np.sqrt(coeffs[1]**2 + coeffs[2]**2)
    posang = np.rad2deg(np.arctan2(gy, -gx))# * 180 / pi
    # chi^2_red
    red_chisq = np.sum((dv - vp)**2 * wt) / (len(dv) - 3.)

    # if boot > 0:
    #     # Calculate uncertainty by randomly changing the velocities 
    #     # and recalculate the fit
    #     grad_list = np.zeros(boot + 1)
    #     posang_list = np.zeros(boot + 1)
    #     vc_list = np.zeros(boot + 1)
    #     for i in range(boot):
    #         v_noise = np.random.normal(0., 1., npts) * ev
    #     # Call again this function
    #         result_i = VFit_OverAll(Xold,Yold,V+v_noise,Sigma_v,
    #                                 distance,resolution/conv,boot=0)
    #         grad_list[i] = results_i['grad'] # grad_aux
    #         posang_list[i] = results_i['posang'] # posang_aux
    #         vc_list[i] = results_i['vc']
    #     vc_err = np.std(vc_list)
    #     grad_err = np.std(grad_list)
    #     paerr = np.std(posang_list)
    # else:
    vc_err = 0.
    grad_err = np.sqrt((gx * errx)**2 + (gy * erry)**2) / grad
    grad_err = np.sqrt((gx * errx)**2 + (gy * erry)**2 
                     + 2 * gx * gy * covar[2, 1]) / grad
    paerr = np.rad2deg(np.sqrt((gx / (gx**2 + gy**2))**2 * erry**2 +
                     (gy / (gx**2 + gy**2))**2 * errx**2))
    paerr = np.rad2deg(np.sqrt((gx / (gx**2 + gy**2))**2 * erry**2 +
                     (gy / (gx**2 + gy**2))**2 * errx**2 
                     - 2 * gx * gy / (gx**2 + gy**2)**2 * covar[2, 1]))
    chisq = red_chisq
    vp += v_mean
    # re
    results = {'grad': grad, 'grad_err': grad_err, 
               'posang': posang, 'paerr': paerr, 
               'chisq': chisq, 'vc': vc, 'vc_err': vc_err,
               'v_pred': vp}

    return results
    #return grad[0], grad_err[0], posang, paerr, chisq, vc, vc_err


def vfit_local(x, y, v, ev, distance=300.0*u.pc, 
    width=45*u.arcsec, nmin=7):
    """
    Function to calculate local gradients
    nmin is the number of pixels, a value of 7 is a good value in the case
        of a Nyquist sampled dataset.
    """
    
    grad = numpy.empty(x.shape)
    grad_err = numpy.empty(x.shape)
    posang = numpy.empty(x.shape)
    paerr = numpy.empty(x.shape)
    vc = numpy.empty(x.shape)
    vp = numpy.empty(x.shape)
    chisq = numpy.empty(x.shape)
    #
    grad[:] = np.nan * u.km / u.s / u.pc
    grad_err[:] = np.nan * u.km / u.s / u.pc
    posang[:] = np.nan * u.deg
    paerr[:] = np.nan * u.deg
    vc[:] = np.nan * u.km / u.s
    vp[:] = np.nan * u.km / u.s
    chisq[:] = np.nan

    for index, (x_i, y_i) in enumerate(zip(x, y)):
        dxy = np.sqrt((x - x_i)**2 + (y - y_i)**2)
        gd = (dxy <= width)
        if np.sum(gd) > nmin:
            res_grad = vfit(x[gd], y[gd], v[gd], ev[gd], distance=distance)
            grad[index] = res_grad['grad']
            grad_err[index] = res_grad['grad_err']
            posang[index] = res_grad['posang']
            paerr[index] = res_grad['paerr']
            vc[index] = res_grad['vc']
            vp[index] = res_grad['vp']
            chisq[index] = res_grad['chisq']

    results = {'grad': grad, 'grad_err': grad_err, 
           'posang': posang, 'paerr': paerr, 
           'chisq': chisq, 'vc': vc, 
           'v_pred': vp}

    return results


def vfit_image(file_vc, file_evc, distance=300.0*u.pc, 
                width=None, n_oversample=7):
    vc, hd = fits.getdata(file_vc, header=True)
    evc = fits.getdata(file_evc)
    x = np.arange(hd['NAXIS1']) * hd['CDELT1'] * u.deg
    y = np.arange(hd['NAXIS2']) * hd['CDELT2'] * u.deg
    xv, yv = np.meshgrid(x, y, indexing='xy')
    gd = np.finite(vc * evc)
    if width == None:
        width = 2 * hd['BMAJ'] * u.deg

    beam_pix = np.abs((hd['BMAJ'] * hd['BMIN']) / 
                     (4 * hd['CDELT1'] * hd['CDELT2']))
    result_all = vfit(xv[gd], yv[gd], v[gd], ev[gd], distance=distance, 
        width=width, nmin=n_oversample * beam_pix)
    result = vfit_local(xv[gd], yv[gd], v[gd], ev[gd], distance=distance, 
        width=width, nmin=n_oversample * beam_pix)

    hd['vgrad'] = result_all['grad']
    hd['vgrad_err'] = result_all['grad_err']
    hd['vgrad_PA'] = result_all['posang']
    hd['vgrad_PA_err'] = result_all['paerr']
    hd['vgrad_vc'] = result_all['vc']

