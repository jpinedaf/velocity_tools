import numpy as np
import astropy.units as u
from astropy import wcs
from astropy.coordinates import SkyCoord
import warnings


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
    """
    warnings.warn("convolve_Vlsr has been moved to keplerian_field. Please use keplerian_field.convolve_Vlsr instead. It will be removed in the next release.", DeprecationWarning)
    from . import keplerian_field
    return keplerian_field.convolve_Vlsr(V_lsr, header, dilation=False)


def generate_Vlsr(radius, angle, inclination=42.*u.deg,
                  R_out=300.*u.au, Mstar=2.2*u.Msun, Vc=5.2*u.km/u.s) -> None:
    warnings.warn("generate_Vlsr has been deprecated. Please use keplerian_field.generate_Vlsr instead. It will be removed in the next release.", DeprecationWarning)
    return None

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
    lat.shape = xx.shape
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
    results.pre_lat = lat
    results.pre_lon = lon
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


def average_profile(x, y, dx, dy=None, log=False, oversample=1.):
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
