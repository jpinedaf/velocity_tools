import numpy as np
from astropy.io import fits
import astropy.units as u
from astropy import wcs
from astropy.coordinates import SkyCoord

class result_container:
    """ Function to create class container
    """
    pass


def convolve_Vlsr(V_lsr, header):
    """ 
    It Convolves a pure theoretical Vlsr map with a requested beam.
    The beam is setup using the FITS header of the expected observation.
    The convolution would mimick (at least at first order) the effect of a 
    finite beam size in the observations.
    The header is also used to convert the beam into pixel units before 
    convolving the map

    param :
    Vlsr : image with Velocity map to be convolved. It handles astropy.units.
    header : FITS header with beam and pixel size information

    TODO:
    -pass a Beam structure instead of the FITS header, to make it more flexible
    """
    from astropy.convolution import convolve
    from radio_beam import Beam
    pixscale = np.abs(header['cdelt1'])*u.Unit(header['cunit1'])
    my_beam = Beam.from_fits_header(header)
    my_beam_kernel = my_beam.as_kernel(pixscale)
    return convolve(V_lsr, my_beam_kernel, boundary='fill', fill_value=np.nan)


def generate_Vlsr( radius, angle, inclination=42.*u.deg,
    R_out=300.*u.au, Mstar = 2.2*u.Msun, Vc= 5.2*u.km/u.s):
    """
    Keplerian velocity field, for a star of mass=Mstar, inclination angle with 
    respect of the sky of inclination.
    The position in the disk is described in polar coordinates, radius and angle.
    The central velocity of the star is Vlsr, 
    and the maximum outer disk radius is Rout
    It makes full use of astropy.units.
    
    param :
    radius : radius in distance units (e.g. u.au or u.pc)
    angle : position angle with respect to major axis

    Mstar : central stellar mass (e.g. u.Msun)
    Vlsr : velocity of the star (u.km/u.s)
    inclination : with units (e.g u.deg or u.rad)
    R_out : Maximum radius of the disk. Position outside this radius are blanked
    """
    Kep_velo = 29.78 * np.sqrt((Mstar/u.Msun / (radius/u.au)).decompose()) * u.km/u.s
    Kep_velo *= np.sin(inclination) * np.cos(angle)
    Kep_velo += Vc
    Kep_velo[radius > R_out] = np.nan
    return Kep_velo


def generate_offsets(header, ra0, dec0, frame='fk5',
                     pa_angle=142.*u.deg, inclination=42.*u.deg):
    """
    Main geometry: major axis in deprojected lon variable, while lat is the minor axis

    :param header: FITS header of the file to be used
    :param ra0: RA of reference point for offset calculation (with units)
    :param dec0: Dec of reference point for offset calculation (with units)
    :param frame: coordinate frame for reference point (default='fk5')
    :param pa_angle: PA angle in deg
    :param inclination: inclination angle in deg
    :return: a structure with the radius and position angle in the deprojected
       coordinates, and the x and y offsets also in deprojected coordinates
    """
    #
    center = SkyCoord(ra0, dec0, frame=frame)
    # Load WCS 
    w = wcs.WCS(header)
    # Create xy array and then coordinates
    xx, yy = np.meshgrid(np.arange(header['naxis1']),
                         np.arange(header['naxis2']))
    world = w.wcs_pix2world(xx.flatten(), yy.flatten(), 0)
    radec = SkyCoord(world[0]*u.deg, world[1]*u.deg, frame=frame)
    radec_off = radec.transform_to(center.skyoffset_frame())
    #
    # Ra = Lon, Dec = Lat
    #
    lon = radec_off[:].lon
    lat = radec_off[:].lat
    lon.shape = xx.shape
    lat.shape = yy.shape
    # Rotate the axes
    c, s = np.cos(pa_angle), np.sin(pa_angle)
    lat_pa = c*lat + s*lon
    lon_pa = -s*lat + c*lon
    # Deprojection 
    # Major axis in in Lon direction
    lon_pa /= np.cos(inclination)
    # deprojected radius
    dep_angle = np.sqrt(lat_pa**2 + lon_pa**2)
    # deprojected angle
    angle_pa = np.arctan2(lon_pa, lat_pa)
    # Store results on class
    results = result_container()
    results.r = dep_angle
    results.theta = angle_pa
    results.lat = lat_pa
    results.lon = lon_pa
    return results


def mask_velocity(cube, v_map, v_width=1.0*u.km/u.s):
    """
    Returns a mask with the pixels in the channel mask within v_width of the
     expected velocity map (v_map). Values of 1.0s and 0.0s

    :param v_width:
    :param cube: SpectralCube cube object to work on.
    :param v_map: Centroid velocity map in velocity units.
    :return: Mask with shape the same as of the input cube,
      1=in the mask, 0=outside of the mask.
    """
    cube2 = cube.with_spectral_unit(u.km/u.s, velocity_convention='radio')
    vaxis = cube2.spectral_axis
    # Load Keplerian velocity model and give proper units
    vmask = np.zeros(cube2.shape)
    for ii in np.arange(0, vaxis.size):
        mask_i = np.abs(v_map - vaxis[ii]) < v_width
        vmask[ii, :, :] = mask_i
    return vmask.astype(np.float)


def vfit_grad(X, Y, V, V_err, nmin=7):
    """
    Function to fit a single gradient to a velocity field.
    It assumes solid body rotation, and it uses the velocity uncertainty.

    X : Off-Set. With appropriate units (e.g. deg)
    Y : Off-Set. With appropriate units (e.g. deg)
    V : Radial velocity. With appropriate units (e.g. km/s)
    V_err : Uncertainty in radial velocity. With appropriate units (e.g. km/s)

    param :
    nmin : Minimum number of pixels required to carry out the fit. Default is 7, 
    which is appropriate for single dish data that is Nyquist sampled

    OUTPUTS:
    Vc:       Mean centroid velocity in km/s
    Vc_err:   Uncertainty of the mean centroid velocity in km/s
    Grad:     Velocity gradient in km/s/pc
    Grad_Err: Uncertainty associated to the velocity gradient (km/s/pc)
    PosAng:   Position angle of the fitted gradient, in degrees
    PAErr:    Uncertainty of the position angle (degrees)
    ChiSq:    Chi^2
    """

    # Xold = copy(X)
    # Yold = copy(Y)
    npts = X.shape
    
    if npts < nmin:
        results = result_container()
        results.Grad = np.nan
        results.Grad_err = np.nan
        results.GradPA = np.nan
        results.GradPA_err = np.nan
        results.Vc = np.nan
        results.Vc_err = np.nan
        return results #
    wt = 1/(V_err**2)
# Obtain total weight, and average (x,y,v) to create new variables (dx,dy,dv)
# which provide a lower uncertainty in the fit.
#
    sumWt = np.sum(wt)
    x_mean = np.sum(X*wt)/sumWt
    y_mean = np.sum(Y*wt)/sumWt
    v_mean = np.sum(V*wt)/sumWt
    dx = (X-x_mean)  # [mask]  # remove mean value from inputs
    dy = (Y-y_mean)  # [mask]  # to reduce fit uncertainties
    dv = (V-v_mean)  # [mask]  #
    M = [[np.sum(wt),   np.sum(dx*wt),    np.sum(dy*wt)], 
        [np.sum(dx*wt), np.sum(dx**2*wt), np.sum(dx*dy*wt)], 
        [np.sum(dy*wt), np.sum(dx*dy*wt), np.sum(dy**2*wt)]]
    #
    from scipy import linalg
    try:
        covar = linalg.inv(M)
    except IOError:
        import sys
        sys.exit('Singular matrix: no solution returned')
    coeffs = np.dot(covar,[[np.sum(dv*wt)], [np.sum(dx*dv*wt)],[np.sum(dy*dv*wt)]])
    #
    errx = np.sqrt(covar[1, 1])
    erry = np.sqrt(covar[2, 2])
    #
    gx = coeffs[1][0]
    gy = coeffs[2][0]
    #
    vc = coeffs[0]+v_mean
    vp = coeffs[0]+coeffs[1]*dx+coeffs[2]*dy
    grad = np.sqrt(coeffs[1]**2+coeffs[2]**2)
    posang = np.arctan2(gy, -gx)*180/pi
    #
    # red_chisq = np.sum( (dv-vp)**2*wt)/(np.len(dv)-3.)

    vc_err = 0.
    # grad_err = np.sqrt((gx*errx)**2+(gy*erry)**2)/grad
    grad_err = np.sqrt((gx*errx)**2+(gy*erry)**2+2*gx*gy*covar[2,1])/grad
    # paerr = 180/pi*sqrt((gx/(gx**2+gy**2))**2*erry**2 +
    #                      (gy/(gx**2+gy**2))**2*errx**2)
    paerr = 180/pi*sqrt((gx/(gx**2+gy**2))**2*erry**2 +
                         (gy/(gx**2+gy**2))**2*errx**2 - 2*gx*gy/(gx**2+gy**2)**2*covar[2,1])
    # chisq = red_chisq
    vp += v_mean
    #
    results = result_container()
    results.Grad = grad
    results.Grad_err = grad_err
    results.GradPA = posang
    results.GradPA_err = paerr
    results.Vc = vc
    results.Vc_err = vc_err
    return results


def average_profile( x, y, dx, dy=None, log=False, oversample=1.):
    """ 
    Averaging function to create a radial profile.

    It returns the average of y in the range [x, x+dx], where the bin size are 
    constant in linear space. 
    If log=True, then the bin sizes are constant in log space, and the average is 
    calculated in the range [log(x), log(x)+dx].

    oversample = 1. 
            This is the number to correct for the correlated number of pixesl within a beam. 
            Default is 1.
            Normal usage should be (beamsize/pixelsize)**2
    
    returns: xbin, ybin, dxbin, dybin
    xbin : the middle of the bin
    ybin : the mean of the y values in the given bin
    dxbin : the width of the bin
    dybin : the uncertainty on the mean value calculation. This is the standard deviation of the points.

    """
    if log == False:
        # Linear space
        xx = x
    else:
        # log space
        xx = np.log(x)

    xmin = np.min(xx)
    xmax = np.max(xx)
    n_bin = int(np.ceil((xmax-xmin)/dx))
    xbin = np.zeros(n_bin)
    dxbin = np.zeros(n_bin)
    ybin = np.zeros(n_bin)
    dybin = np.zeros(n_bin)
    for i in range(n_bin):
        idx = np.where((xx > xmin+dx*i) & (xx < xmin+dx*(i+1)))
        xbin[i] = xmin + dx*(i+0.5)
        #
        if dy is None:
            ybin[i] = np.average(y[idx])
            dybin[i] = np.std(y[idx]) / np.sqrt(y[idx].size / oversample)
        else:
            ybin[i] = np.average(y[idx], weights=1./dy[idx]**2)
            dybin[i] = 1. / np.sqrt(np.sum(1. / dy[idx]**2))
        dxbin[i] = dx*0.5
    if log == False:
        return xbin, ybin, dxbin, dybin
    else:
        return np.exp(xbin), ybin, dxbin, dybin
