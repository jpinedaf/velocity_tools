import pytest
import astropy.units as u
import numpy as np
from astropy.io import fits
from velocity_tools.keplerian_field import convolve_Vlsr
import velocity_tools.keplerian_field


@pytest.fixture
def generate_dummy_file(tmp_path) -> str:
    """ """
    file_in = tmp_path / "test_keplerian.fits"
    ra0 = 15 * (3 + (30 + 25.6 / 60.0) / 60.0)
    dec0 = 32 + (10 + 25.6 / 60.0) / 60.0
    # create dummy FITS file
    dx = 1500
    data = np.zeros((2 * dx + 1, 2 * dx + 1))
    # hdu container
    hdu = fits.PrimaryHDU(data)
    pixel = 0.01 * (u.arcsec).to(u.deg)
    #
    hdu.header["CTYPE1"] = "RA---TAN"
    hdu.header["CTYPE2"] = "DEC--TAN"
    hdu.header["CUNIT1"] = "deg"
    hdu.header["CUNIT2"] = "deg"
    hdu.header["CDELT1"] = -1 * pixel
    hdu.header["CDELT2"] = pixel
    hdu.header["CRVAL1"] = ra0
    hdu.header["CRVAL2"] = dec0
    hdu.header["CRPIX1"] = dx * 1.0
    hdu.header["CRPIX2"] = dx * 1.0
    hdu.header["RADESYS"] = "FK5"
    hdu.header["EQUINOX"] = 2000.0
    hdu.header["BUNIT"] = "km/s"
    hdu.header["BMAJ"] = 4 * pixel
    hdu.header["BMIN"] = 4 * pixel
    hdu.header["BPA"] = 0.0
    hdu.writeto(file_in, overwrite=True)

    return str(file_in)

def test_convolve_Vlsr(generate_dummy_file) -> None:
    """Test the convolve_Vlsr function."""

    # Load the dummy FITS file
    file_in = generate_dummy_file
    with fits.open(file_in) as hdul:
        header = hdul[0].header
        data = hdul[0].data

    # Create a dummy V_lsr map
    V_lsr = np.random.random(data.shape)

    # Test convolve_Vlsr without dilation
    convolved_data = velocity_tools.keplerian_field.convolve_Vlsr(
        V_lsr, header, dilation=False)
    assert convolved_data.shape == V_lsr.shape
    assert not np.isnan(convolved_data).all()

    # Test convolve_Vlsr with dilation
    convolved_data_dilated = velocity_tools.keplerian_field.convolve_Vlsr(
        V_lsr, header, dilation=True)
    assert convolved_data_dilated.shape == V_lsr.shape
    assert not np.isnan(convolved_data_dilated).all()

def test_generate_Vlsr(generate_dummy_file) -> None:
    # Load the dummy FITS file
    file_in = generate_dummy_file
    with fits.open(file_in) as hdul:
        header = hdul[0].header
        data = hdul[0].data
    ra0 = header['CRVAL1'] * u.deg
    dec0 = header['CRVAL2'] * u.deg
    R_out = 300.*u.au
    # Generate the Vlsr map
    inc0 = 42.*u.deg
    # inc2 = 60.*u.deg
    results = velocity_tools.keplerian_field.generate_Vlsr(
        header, ra0, dec0, PA_Angle=142.*u.deg, inclination=inc0, 
        distance=110.02*u.pc, R_out=R_out, Mstar=2.2*u.Msun, 
        Vc=5.2*u.km/u.s)
    # results2 = velocity_tools.keplerian_field.generate_Vlsr(
    #     header, ra0, dec0, PA_Angle=142.*u.deg, inclination=inc2, 
    #     distance=110.02*u.pc, R_out=R_out, Mstar=2.2*u.Msun, 
    #     Vc=5.2*u.km/u.s, do_plot=False)
    assert results is not None
    assert results.v.shape == data.shape
    # Check that values outside R_out are NaN
    assert np.isnan(results.v[results.r > R_out]).all()


# def test_generate_Vlsr() -> None:
#     rdisk = 300 * u.au
#     radius = np.linspace(1, 10, endpoint=True, num=5) * u.au
#     angle = 0 * u.deg
#     inclination = 0 * u.deg
#     vc0 = 5 * u.km / u.s

#     vel0 = velocity_tools.coordinate_offsets.generate_Vlsr(
#         radius,
#         angle,
#         inclination=inclination,
#         R_out=rdisk,
#         Mstar=1.0 * u.Msun,
#         Vc=0 * u.km / u.s,
#     )
#     vel1 = velocity_tools.coordinate_offsets.generate_Vlsr(
#         radius, 
#         angle, 
#         inclination=inclination, 
#         R_out=rdisk, 
#         Mstar=1.0 * u.Msun, 
#         Vc=vc0
#     )
#     assert ((vel1 - vel0) == vc0).all
