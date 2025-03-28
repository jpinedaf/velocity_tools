import astropy.units as u
import numpy as np
from astropy.io import fits
import velocity_tools.coordinate_offsets as coordinate_offsets

import pytest


def test_generate_offsets() -> None:
    #
    # check header is valid
    data = np.ones(shape=(11, 11))
    hdu = fits.PrimaryHDU(data)
    ra0 = 3.5 * u.deg
    dec0 = 30 * u.deg
    cdelt = 10 * u.arcsec.to(u.deg)
    hdu.header["CRVAL1"] = (ra0.value, "deg")
    hdu.header["CRVAL2"] = (dec0.value, "deg")
    hdu.header["CRPIX1"] = 5.0
    hdu.header["CRPIX2"] = 5.0
    hdu.header["CDELT1"] = -cdelt
    hdu.header["CDELT2"] = cdelt
    hdu.header["CTYPE1"] = "RA---TAN"
    hdu.header["CTYPE2"] = "DEC--TAN"

    results = coordinate_offsets.generate_offsets(
        hdu.header, ra0, dec0, frame="fk5", pa_angle=0 * u.deg, inclination=0 * u.deg
    )
    # Deprojection
    # Major axis in in Lon direction
    dxi = np.arange(-5, 6) * -cdelt * u.deg
    dyi = np.arange(-5, 6) * cdelt * u.deg

    dx, dy = np.meshgrid(dxi, dyi)

    # ad  more comparisons
    assert results.lat.shape == (11, 11)
    assert (dx == results.lat).all
    assert (dy == results.lon).all
    assert (np.hypot(dx, dy) == results.r).all


def test_average_profile() -> None:
    from velocity_tools.coordinate_offsets import average_profile

    # Generate some data
    def function_profile(x):
        return 10 + x**-2
    n_sample = 1000
    x_axis = 1 + 10 * np.random.rand(n_sample)
    y_axis = function_profile(x_axis) + np.random.rand(n_sample) * 0.1

    # Calculate the average profile
    xbin, ybin, dxbin, dybin = coordinate_offsets.average_profile(x_axis, y_axis, dx=1.0)
    assert len(xbin) == len(ybin)
    assert len(xbin) == len(dxbin)
    assert len(xbin) == len(dybin)
    assert pytest.approx(ybin, rel=0.05) == function_profile(xbin)