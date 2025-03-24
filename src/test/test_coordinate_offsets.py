import astropy.units as u
import numpy as np
from astropy.io import fits
import velocity_tools.coordinate_offsets

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

    results = velocity_tools.coordinate_offsets.generate_offsets(
        hdu.header, ra0, dec0, frame="fk5", pa_angle=0 * u.deg, inclination=0 * u.deg
    )
    # Deprojection
    # Major axis in in Lon direction
    dxi = np.arange(-5, 6) * -cdelt * u.deg
    dyi = np.arange(-5, 6) * cdelt * u.deg

    dx, dy = np.meshgrid(dxi, dyi)

    # results.r = dep_angle
    # results.theta = angle_pa
    # results.lat = lat_pa
    # results.lon = lon_pa
    # ad  more comparisons
    assert results.lat.shape == (11, 11)
    assert (dx == results.lat).all
    assert (dy == results.lon).all
    assert (np.hypot(dx, dy) == results.r).all


def test_generate_Vlsr() -> None:
    rdisk = 300 * u.au
    radius = np.linspace(1, 10, endpoint=True, num=5) * u.au
    angle = 0 * u.deg
    inclination = 0 * u.deg
    vc0 = 5 * u.km / u.s

    vel0 = velocity_tools.coordinate_offsets.generate_Vlsr(
        radius,
        angle,
        inclination=inclination,
        R_out=rdisk,
        Mstar=1.0 * u.Msun,
        Vc=0 * u.km / u.s,
    )
    vel1 = velocity_tools.coordinate_offsets.generate_Vlsr(
        radius, angle, inclination=inclination, R_out=rdisk, Mstar=1.0 * u.Msun, Vc=vc0
    )
    assert ((vel1 - vel0) == vc0).all
