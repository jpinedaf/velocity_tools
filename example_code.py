import keplerian_field
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.io import fits
import numpy as np

#ra1 = ra0=15*(11+(33 + 25.318652/60.)/60.)*u.deg
#dec1=-1*(70+(11 + 41.23173/60.)/60.)*u.deg
def example_Vlsr(file_in='fits_files/test_file.fits'):
    """ How to use this thing
    """
    # Central coordinate
    ra0=15*(11+(33 + 25.321/60.)/60.)*u.deg
    dec0=-1*(70+(11 + 41.209/60.)/60.)*u.deg
    header=fits.getheader(file_in)
    results=keplerian_field.generate_Vlsr( header, ra0, dec0,    
        file_out='fits_files/test_Vc.fits', PA_Angle=(142.+180)*u.deg,
        # Inclination angle in Degrees
        inclination=42.*u.deg,
        # Distance in pc
        distance=110.02*u.pc,
        # Outer radius in au
        R_out=3.2*110.02*u.au,   #R_out=300.*u.au,
        # Stellar Mass in Solar Masses
        Mstar = 2.2*u.Msun,
        # V_lsr
        Vc= 5.64*u.km/u.s, do_plot=False)
        # Vc= 5.2*u.km/u.s, do_plot=False)
    conv_Vlsr = keplerian_field.convolve_Vlsr( results.v, header)
    plt.figure(figsize=(12,6))
    plt.subplot(1, 3, 1)
    plt.imshow(results.r.value, origin='lowest', interpolation='none')
    plt.title('Deprojected radius, $r_{PA,i}$')
    plt.subplot(1, 3, 2)
    plt.imshow(results.theta.value, origin='lowest', interpolation='none')
    plt.title('Deprojected Angle, $theta$')
    plt.subplot(1, 3, 3)
    # plt.imshow(results.v.value, origin='lowest', cmap='RdYlBu_r', interpolation='none')
    plt.imshow(conv_Vlsr, origin='lowest', cmap='RdYlBu_r', interpolation='none')
    plt.title('Projected $V_{Kep}$')

    header2=header
    header2['BUNIT']='km/s'
    fits.writeto( 'fits_files/HD100546_Vc.fits', results.v.value, header2, overwrite=True)
    header2=header
    header2['BUNIT']='km/s'
    fits.writeto( 'fits_files/HD100546_Vc_conv.fits', conv_Vlsr, header2, overwrite=True)

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

if False:
    # Load data and convert cube into spectral axis of km/s
    example_Vlsr(file_in='fits_files/HD100546_masked_Vc.fits')
    from spectral_cube import SpectralCube
    cube = SpectralCube.read('fits_files/HD100546_12CO_mscale_cube_3D.fits')
    cube2=cube.with_spectral_unit(u.km/u.s, velocity_convention='radio')
    vaxis=cube2.spectral_axis
    # Load keplerian velocity model and give proper units
    vkep=fits.getdata('fits_files/HD100546_Vc_conv.fits')*u.km/u.s
    vmask=np.zeros( cube2.shape)
    v_width=(1.0*u.km/u.s)
    for ii in np.arange(0,vaxis.size):
        mask_i=np.abs(vkep-vaxis[ii])<v_width
        vmask[ii,:,:]=mask_i
    header_v=fits.getheader('fits_files/HD100546_12CO_mscale_cube_3D.fits')
    file_mask_out='fits_files/test_mask_1kms.fits'
    fits.writeto(file_mask_out,vmask.astype(np.float), header_v, overwrite=True)

    vmask2=np.zeros( cube2.shape)
    v_width2=(2.0*u.km/u.s)
    for ii in np.arange(0,vaxis.size):
        mask_i2=np.abs(vkep-vaxis[ii])<v_width2
        vmask2[ii,:,:]=mask_i2
    header_v=fits.getheader('fits_files/HD100546_12CO_mscale_cube_3D.fits')
    file_mask_out='fits_files/test_mask_2kms.fits'
    fits.writeto(file_mask_out,vmask2.astype(np.float), header_v, overwrite=True)

    file_mask_out='fits_files/test_mask_1kms.fits'
    my_mask=fits.getdata(file_mask_out)
    cube2_mask = cube2.with_mask(my_mask.astype(np.bool))
    m0 = cube2_mask.moment(order=0)
    m1 = cube2_mask.moment(order=1)
    m0.write('fits_files/test_mom0_masked_1kms.fits')
    m1.write('fits_files/test_mom1_masked_1kms.fits')

    file_mask_out='fits_files/test_mask_2kms.fits'
    my_mask=fits.getdata(file_mask_out)
    cube2_mask = cube2.with_mask(my_mask.astype(np.bool))
    m0 = cube2_mask.moment(order=0)
    m1 = cube2_mask.moment(order=1)
    m0.write('fits_files/test_mom0_masked_2kms.fits')
    m1.write('fits_files/test_mom1_masked_2kms.fits')
        # try:
        #     vmask[mask_i]=1
        # except:
