.. _doc_vgrad:

Velocity gradient calculation
=============================

This module calculates the velocity gradient of a given velocity field. The velocity field is assumed to be a solid rotation. 
The velocity gradient is calculated by fitting a plane to the velocities provided, 
this is the same as used in `Goodman et al. (1993) <https://ui.adsabs.harvard.edu/abs/1993ApJ...406..528G>`_, 
`Caselli et al. (2002) <https://ui.adsabs.harvard.edu/abs/2002ApJ...572..238C>`_, and 
`Pineda et al. (2010) <https://ui.adsabs.harvard.edu/abs/2010ApJ...712L.116P>`_. 
The velocity gradient is given in km/s/pc.

The module allows for the determination of the global velocity gradient in the images. 
This called as follows::

    import numpy as np
    import astropy.units as u
    from velocity_tools import vgrad
    from astropy.io import fits
    file_vc = 'test_vel.fits'
    file_evc = 'test_vel_err.fits'
    vc = fits.getdata(file_vc)
    evc = fits.getdata(file_evc)

    ny_pix, nx_pix = vc.shape
    x_t = -1* np.arange(nx_pix) * u.arcsec
    y_t = np.arange(ny_pix) * (u.arcsec)
    xv, yv = np.meshgrid(x_t, y_t, indexing='xy')
    result_all = vgrad.vfit(xv, yv, vc, evc, distance=10.0*u.pc)

where the the result_all is a structure with results from the fit, including uncertainties. 
The gradient is in units of km/s/pc, while the position angle is in degrees measured 
East-from-North. Velocities are reported in km/s.::
    
    result_all.keys()
    result_all = {'grad': grad, 'grad_err': grad_err, 
               'posang': posang, 'paerr': paerr, 
               'chisq': chisq, 'vc': vc, 'vc_err': vc_err,
               'v_pred': vp}

The module also allows for the determination of the local velocity gradient in an image.
This is called as follows::

    result = vfit_local(xv, yv, vc, evc, distance=10.0*u.pc, 
                        width=45*u.arcsec, nmin=7.0)

where the result is a structure with results from the fit, including uncertainties, calculated around 
a box of with ``width`` and centered in each pixel. 
The variable ``nmin`` is the minimum number of pixels to consider in the fit, which 
should be used to weed out calculations with not enough pixels to be reliable.
The gradient is in units of km/s/pc, while the position angle is in degrees measured
East-from-North. Velocities are reported in km/s.::

    result.keys()
    result = {'grad': grad, 'grad_err': grad_err, 
               'posang': posang, 'paerr': paerr, 
               'chisq': chisq, 'vc': vc, 'vc_err': vc,
               'v_pred': vp}

It is very common to report both the global and local velocity gradients in a paper.
The global velocity gradient is a measure of the overall velocity gradient in the image, while the local velocity gradient is a measure of the velocity gradient at each pixel.
Therefore, the function ``vift_image`` is provided as a convinient function to 
obtain the results for a whole field providing only the centroid velocity and uncertainties files.::

    results_image = vfit_image(file_vc, file_evc, distance=300.0*u.pc, 
                    width=None, n_oversample=2.0)

where the results_image is a structure with the results of the fit for the whole field.
The variable ``width`` is the width of the box to calculate the local velocity gradient (default=2x Beam Size),
while the variable ``n_oversample`` is the number of beams required as the minimum value in the gradient calculation.
The gradient is in units of km/s/pc, while the position angle is in degrees measured
East-from-North. Velocities are reported in km/s.
The results are stored in a dictionary containing the different result images of the local gradient calculations, 
while the global gradient calculation is stored in the FITS header also returned 
in the result structure.::

    
    results_image.hd['VG']
    results_image.hd['VGerr']
    results_image.hd['VG_PA']
    results_image.hd['VG_PAerr']
    results_image.hd['VG_vc']
    results_image = {'grad': grad, 'grad_err': grad_err, 
        'posang': grad_PA, 'paerr': grad_PA_err, 
        'vc': vc, 'header':hd}
