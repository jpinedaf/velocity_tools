import keplerian_field
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.io import fits

def example_Vlsr(file_in='fits_files/test_file.fits'):
    """ How to use this thing
    """
    # Central coordinate
    ra0=15*(11+(33 + 25.321/60.)/60.)*u.deg
    dec0=-1*(70+(11 + 41.209/60.)/60.)*u.deg
    # ra0=15*(3+(30+ 25.6/60.)/60.)*u.deg
    # dec0= (32+(10+ 25.6/60.)/60.)*u.deg
    header=fits.getheader(file_in)
    # dep_radius, dep_angle, Kep_velo, dmaj, dmin = 
    # results=generate_Vlsr( header, ra0, dec0,
    #     file_out='fits_files/test_Vc.fits', PA_Angle=90.*u.deg, #142.*u.deg,
    results=keplerian_field.generate_Vlsr( header, ra0, dec0,    
        file_out='fits_files/test_Vc.fits', PA_Angle=142.*u.deg,
        # Inclination angle in Degrees
        inclination=42.*u.deg,
        # Distance in pc
        distance=110.02*u.pc,
        # Outer radius in au
        R_out=300.*u.au,
        # Stellar Mass in Solar Masses
        Mstar = 2.2*u.Msun,
        # V_lsr
        Vc= 5.2*u.km/u.s, do_plot=False)

    plt.figure(figsize=(12,6))
    plt.subplot(1, 3, 1)
    plt.imshow(results.r.value, origin='lowest', interpolation='none')
    plt.title('Deprojected radius, $r_{PA,i}$')
    plt.subplot(1, 3, 2)
    plt.imshow(results.theta.value, origin='lowest', interpolation='none')
    plt.title('Deprojected Angle, $theta$')
    plt.subplot(1, 3, 3)
    plt.imshow(results.v.value, origin='lowest', cmap='RdYlBu_r', interpolation='none')
    plt.title('Projected $V_{Kep}$')

    header2=header
    header2['BUNIT']='km/s'
    fits.writeto( 'fits_files/HD100546_Vc.fits', results.v.value, header2, overwrite=True)

    header2=header
    header2['BUNIT']='au'
    fits.writeto( 'fits_files/HD100546_Radius.fits', results.r.value, header2, overwrite=True)

    header2=header
    header2['BUNIT']='deg'
    fits.writeto( 'fits_files/HD100546_theta.fits', (results.theta.to(u.deg)).value, header2, overwrite=True)

    # header2=header
    # header2['BUNIT']='rad'
    # fits.writeto( 'fits_files/HD100546_theta.fits', results.theta.value, header2, overwrite=True)
    header2=header
    header2['BUNIT']='au'
    fits.writeto( 'fits_files/HD100546_lat.fits', results.lat.value, header2, overwrite=True)

    header2=header
    header2['BUNIT']='au'
    fits.writeto( 'fits_files/HD100546_lon.fits', results.lon.value, header2, overwrite=True)

from spectral_cube import SpectralCube
cube = SpectralCube.read('adv.fits')