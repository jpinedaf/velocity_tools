import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from .coordinate_offsets import generate_offsets
from radio_beam import Beam
from astropy.convolution import convolve
from skimage.morphology import erosion

class result_container:
    """ Function to create class container
    """
    pass

def convolve_Vlsr(V_lsr, header, dilation=False):
    """ Convolve Vlsr map with beam in header
    If dilation is True, then mask all convolved data beyond
    one beam.
    """
    pixscale=np.abs(header['cdelt1'])
    if header['cunit1'] == 'arcsec':
        pixscale *= u.arcsec
    else:
        pixscale *= u.deg
    my_beam = Beam.from_fits_header(header)
    my_beam_kernel = my_beam.as_kernel(pixscale)
    if dilation == False:
        return convolve(V_lsr, my_beam_kernel, boundary='fill', fill_value=np.nan)
    else:
        v_convolved = convolve(V_lsr, my_beam_kernel, boundary='fill', fill_value=np.nan)
        nimage = np.ceil(np.abs(header['BMAJ'] / header['CDELT1'])).astype(int)
        if nimage % 2 == 0:
            nimage += 1
        # creates box with zeros with the size of the major axis, 
        # and place a 1 at the center
        center = (nimage - 1) // 2
        box_beam = np.zeros((nimage, nimage))
        box_beam[center, center] = 1
        # convolve box with beam and then threshold it to half-power
        beam_shape = convolve(box_beam, my_beam_kernel)
        beam_footprint = ((beam_shape / beam_shape.max()) > 0.5) + 0
        bad = erosion(np.isnan(V_lsr), beam_footprint)
        v_convolved[bad] = np.nan
        return v_convolved 



def generate_Vlsr(header, ra0, dec0, frame='fk5',
    PA_Angle=142.*u.deg, inclination=42.*u.deg, distance=110.02*u.pc,
    R_out=300.*u.au, Mstar=2.2*u.Msun, Vc=5.2*u.km/u.s, do_plot=False):
    """
    Main geometry: major axis in deprojected lon variable, while lat is the minor axis

    """
    #
    results = generate_offsets(header, ra0, dec0, frame=frame,
                     pa_angle=PA_Angle, inclination=inclination)
    # center = SkyCoord(ra0, dec0, frame='fk5')
    # Load WCS 
    # w = wcs.WCS(header)
    # Create xy array and then coordinates
    # x = np.arange(header['naxis1'])
    # y = np.arange(header['naxis2'])
    #
    # epsilon will be determined as the pixel size
    epsilon = (np.abs(header['cdelt1'])*u.deg).to('',
                equivalencies=u.dimensionless_angles())*distance.to(u.au)
    # xx, yy = np.meshgrid(x, y)
    # world = w.wcs_pix2world(xx.flatten(), yy.flatten(), 0)

    # radec = SkyCoord(world[0]*u.deg, world[1]*u.deg, frame='fk5')
    # radec_off = radec.transform_to(center.skyoffset_frame())
    #
    # Ra = Lon, Dec = Lat
    #
    # lon = radec_off[:].lon
    # lat = radec_off[:].lat
    # lon.shape = xx.shape
    # lat.shape = xx.shape
    # This plots the offset in x and y
    # This is all good now
    # c, s = np.cos(PA_Angle), np.sin(PA_Angle)
    # Rotate the axes
    # Ra = Lon, Dec = Lat
    # lat_PA = c*lat + s*lon
    # lon_PA = -s*lat + c*lon
    # Deprojection 
    # Major axis in in RA direction
    # lon_PA /= np.cos(inclination)
    # dep_angle = np.sqrt(lat_PA**2 + lon_PA**2)
    # angle_PA = np.arctan2(lon_PA, lat_PA)
    #
    # dep_angle
    angle_PA = results.theta
    dep_radius = np.clip(results.r.to('',
                         equivalencies=u.dimensionless_angles())*distance.to(u.au),
                         epsilon.to(u.au), None)
    Kep_velo = 29.78 * np.sqrt(Mstar/u.Msun / (dep_radius/u.au)) \
               * np.abs(np.sin(inclination)) * np.cos(angle_PA)*u.km/u.s
    Kep_velo += Vc
    Kep_velo[dep_radius > R_out] = np.nan
    #
    if do_plot:
        # Plot
        axes_extent = [(results.pre_lon.to(u.arcsec).value).max(),
                       (results.pre_lon.to(u.arcsec).value).min(),
                       (results.pre_lat.to(u.arcsec).value).min(),
                       (results.pre_lat.to(u.arcsec).value).max()]
        plt.ion()
        plt.figure(figsize=(12, 6))
        plt.subplot(1, 2, 1)
        plt.imshow(dep_radius.value, origin='lower', interpolation='none',
                   extent=axes_extent)
        dx_line = (results.pre_lat.to(u.arcsec).value).max() * 0.2
        dv_disk = (np.nanmax(np.abs(Kep_velo - Vc))).value
        plt.contour(dep_radius.value,
                    np.linspace(np.nanmin(dep_radius.value),
                                np.nanmax(dep_radius.value), num=10),
                    colors='k', extent=axes_extent)
        plt.plot([0.0, dx_line*np.sin(PA_Angle)],
                 [0.0, dx_line*np.cos(PA_Angle)], color='gray')
        plt.plot([0.0, dx_line*np.sin(PA_Angle+90.*u.deg)],
                 [0.0, dx_line*np.cos(PA_Angle+90.*u.deg)], color='red')
        plt.title('Deprojected radius, $r_{PA,i}$')
        #
        plt.subplot(1, 2, 2)
        plt.imshow(Kep_velo.value, origin='lower', cmap='RdYlBu_r',
                   interpolation='none', extent=axes_extent)
        plt.contour(Kep_velo.value,
                    np.linspace(Vc.value - dv_disk,
                                Vc.value + dv_disk,
                                #np.nanmin(Kep_velo.value),
                                #np.nanmax(Kep_velo.value),
                                num=10),
                    colors='k', extent=axes_extent)
        plt.title('Projected $V_{Kep}$')
    # Store results on class
    results = result_container()
    results.r = dep_radius
    # results.theta = angle_PA
    results.v = Kep_velo
    # results.lat = lat_PA
    # results.lon = lon_PA
    return results
