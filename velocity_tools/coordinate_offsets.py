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
    Convolve pure theoretical Vlsr map with beam in header
    This would mimick (first order) the effect of beam size in the observations
    """
    from astropy.convolution import convolve
    from radio_beam import Beam
    pixscale=np.abs(header['cdelt1'])*u.Unit(header['cunit1'])
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
    Kep_velo = 29.78 * np.sqrt( (Mstar/u.Msun / (radius/u.au)).decompose()) * u.km/u.s
    Kep_velo *= np.sin(inclination) * np.cos(angle)
    Kep_velo += Vc
    Kep_velo[radius > R_out] = np.nan
    return Kep_velo

def generate_offsets( header, ra0, dec0, 
    PA_Angle=142.*u.deg, inclination=42.*u.deg):
    """
    Main geometry: major axis in deprojected lon variable, while lat is the minor axis
    """
    #
    center = SkyCoord(ra0, dec0, frame='fk5')
    # Load WCS 
    w = wcs.WCS(header)
    # Create xy array and then coordinates
    x=np.arange(header['naxis1'])
    y=np.arange(header['naxis2'])
    #
    # epsilon will be determined as the pixel size
    # epsilon= (np.abs(header['cdelt1'])*u.deg).to('', 
    #     equivalencies=u.dimensionless_angles())*distance.to(u.au)
    xx, yy = np.meshgrid(x, y)
    world = w.wcs_pix2world(xx.flatten(), yy.flatten(), 0)
    radec = SkyCoord(world[0]*u.deg, world[1]*u.deg, frame='fk5')
    radec_off = radec.transform_to(center.skyoffset_frame())
    #
    # Ra = Lon, Dec = Lat
    #
    lon=radec_off[:].lon
    lat=radec_off[:].lat
    lon.shape=xx.shape
    lat.shape=yy.shape
    # Rotate the axes
    c, s = np.cos(PA_Angle), np.sin(PA_Angle)
    lat_PA =  c*lat + s*lon
    lon_PA = -s*lat + c*lon
    # Deprojection 
    # Major axis in in Lon direction
    lon_PA /= np.cos(inclination)
    # deprojected radius
    dep_angle=np.sqrt( lat_PA**2 + lon_PA**2)
    # deprojected angle
    angle_PA = np.arctan2(lon_PA, lat_PA)
    # Store results on class
    results = result_container()
    results.r = dep_angle
    results.theta= angle_PA
    results.lat= lat_PA
    results.lon= lon_PA
    return results