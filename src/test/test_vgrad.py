import os
import numpy as np

from astropy.io import fits
from astropy import units as u
from astropy.io import fits
import velocity_tools.vgrad
import pytest

my_distance = 10.0 * u.pc
v0 = 8.5 * u.km / u.s
g0 = 2.0 * u.km / u.s / u.deg
gx = -g0
gy = g0


@pytest.fixture
def velocity_file(tmp_path) -> list:
    """
    make gradient and then
    """
    half_size = 20
    n_pix = 2 * half_size + 1

    x_t = -1 * np.arange(n_pix) * u.arcsec
    y_t = np.arange(n_pix) * (u.arcsec)
    xv, yv = np.meshgrid(x_t, y_t, indexing="xy")
    ev0 = 0.1 * u.km / u.s
    v_t = (xv * gx + yv * gy).to(u.km / u.s) + v0  # km/s
    ev_t = np.empty(v_t.shape) * u.km / u.s
    ev_t[:] = ev0  # km/s
    # Create temporary FITS files for testing
    filenames = []
    pix_size = 1 * u.arcsec.to(u.deg)

    hdu = fits.PrimaryHDU(v_t.value)
    hdu.header["CTYPE1"] = "RA---TAN"
    hdu.header["CTYPE2"] = "DEC--TAN"
    hdu.header["CUNIT1"] = "deg"
    hdu.header["CUNIT2"] = "deg"
    hdu.header["CRPIX1"] = half_size
    hdu.header["CRPIX2"] = half_size
    # somewhere in Perseus
    hdu.header["CRVAL1"] = (15 * 3.5, "deg")
    hdu.header["CRVAL2"] = (31 + 1.0 / 3.0, "deg")
    hdu.header["CDELT1"] = (-1 * pix_size, "deg")
    hdu.header["CDELT2"] = (pix_size, "deg")
    hdu.header["BUNIT"] = ("km / s", "Vlsr")
    # like a single dish
    hdu.header["BMAJ"] = 2 * hdu.header["CDELT2"]
    hdu.header["BMIN"] = hdu.header["BMAJ"]
    hdu.header["BPA"] = 45.0
    filename = tmp_path / "test_vc.fits"
    hdu.writeto(filename)
    filenames.append(str(filename))

    filename = tmp_path / "test_evc.fits"
    # data2 = data * 0.0 + 0.01
    hdu2 = fits.PrimaryHDU(ev_t.value, header=hdu.header)
    hdu2.writeto(filename)
    filenames.append(str(filename))

    return filenames


def test_vfit(velocity_file) -> None:

    vc, hd = fits.getdata(velocity_file[0], header=True)
    evc = fits.getdata(velocity_file[1], header=False)
    x_t = np.arange(hd["NAXIS1"]) * hd["CDELT1"] * (u.deg)
    y_t = np.arange(hd["NAXIS2"]) * hd["CDELT2"] * (u.deg)
    xv, yv = np.meshgrid(x_t, y_t, indexing="xy")
    result_all = velocity_tools.vgrad.vfit(
        xv, yv, vc * u.km / u.s, evc * u.km / u.s, distance=my_distance
    )
    conv_grad = (1 * u.pc.to(u.au) / 10.0 * u.arcsec).to(u.deg) / u.pc
    input_grad = conv_grad * np.sqrt(gx**2 + gy**2)
    input_PA = np.rad2deg(np.arctan2(gx.value, gy.value))

    assert pytest.approx(result_all["grad"]) == input_grad.value
    assert pytest.approx(result_all["posang"]) == input_PA
    assert pytest.approx(result_all["vc"], abs=0.05) == v0.value


def test_vfit_local(velocity_file) -> None:

    vc, hd = fits.getdata(velocity_file[0], header=True)
    evc = fits.getdata(velocity_file[1], header=False)
    x_t = np.arange(hd["NAXIS1"]) * hd["CDELT1"] * (u.deg)
    y_t = np.arange(hd["NAXIS2"]) * hd["CDELT2"] * (u.deg)
    xv, yv = np.meshgrid(x_t, y_t, indexing="xy")
    result_all = velocity_tools.vgrad.vfit_local(
        xv, yv, vc * u.km / u.s, evc * u.km / u.s, distance=my_distance, nmin=5
    )
    conv_grad = (1 * u.pc.to(u.au) / 10.0 * u.arcsec).to(u.deg) / u.pc
    input_grad = conv_grad * np.sqrt(gx**2 + gy**2)
    input_PA = np.rad2deg(np.arctan2(gx.value, gy.value))

    assert pytest.approx(result_all["grad"]) == input_grad.value
    assert pytest.approx(result_all["posang"]) == input_PA
    assert pytest.approx(result_all["vc"], abs=0.05) == v0.value


def test_vfit_image(velocity_file) -> None:
    result_all = velocity_tools.vgrad.vfit_image(
        velocity_file[0], velocity_file[1], distance=my_distance, n_oversample=2
    )
    conv_grad = (1 * u.pc.to(u.au) / 10.0 * u.arcsec).to(u.deg) / u.pc
    input_grad = conv_grad * np.sqrt(gx**2 + gy**2)
    input_PA = np.rad2deg(np.arctan2(gx.value, gy.value))

    assert pytest.approx(result_all["grad"]) == input_grad.value
    assert pytest.approx(result_all["posang"]) == input_PA
    assert pytest.approx(result_all["vc"].value, abs=0.05) == v0.value

    assert pytest.approx(result_all["header"]["VG"]) == input_grad.value
    assert pytest.approx(result_all["header"]["VG_PA"]) == input_PA
    assert pytest.approx(result_all["header"]["VG_vc"], abs=0.05) == v0.value

    with pytest.raises(
        FileNotFoundError, match="File not found: nonexistent_file.fits"
    ):
        # Test with a non-existent file
        velocity_tools.vgrad.vfit_image(
            "nonexistent_file.fits",
            velocity_file[1],
            distance=my_distance,
            n_oversample=2,
        )
    with pytest.raises(
        FileNotFoundError, match="File not found: nonexistent_file.fits"
    ):
        # Test with a non-existent file
        velocity_tools.vgrad.vfit_image(
            velocity_file[0],
            "nonexistent_file.fits",
            distance=my_distance,
            n_oversample=2,
        )
