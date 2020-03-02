import numpy as np
import astropy.units as u
from scipy import optimize
from astropy.constants import G


#
# Implementation of stream lines using the prescription from
# Mendoza et al. (2009)  doi:10.1111/j.1365-2966.2008.14210.x
#

def v_k(radius, mass=0.5 * u.Msun):
    """
    Velocity term that is repeated in all velocity component.
    It corresponds to v_k in Mendoza+(2009)
    """
    return np.sqrt(G * mass / radius).to(u.km / u.s)


def v_r(r, theta, mass=0.5 * u.Msun, theta0=30 * u.deg):
    """
    Radial velocity component
    """
    return v_k(r, mass=mass) * np.sqrt(1 + np.cos(theta) / np.cos(theta0))


def v_phi(r, theta, mass=0.5*u.Msun, theta0=30 * u.deg):
    geom = np.sqrt(1 - np.cos(theta) / np.cos(theta0)) * \
           np.sin(theta0) / np.sin(theta)
    return v_k(r, mass=mass) * geom


def v_theta(r, theta, mass=0.5 * u.Msun, theta0=30 * u.deg):
    """
    Function to calculate the velocity in the theta direction
    :param r: Radius of the streamline in units of r_u
    :param theta:
    :param mass:
    :param theta0:
    :return:
    """
    geom = (np.cos(theta0) - np.cos(theta)) * \
           np.sqrt((np.cos(theta0) + np.cos(theta)) / \
                   (np.cos(theta0) * np.sin(theta) ** 2))
    return v_k(r, mass=mass) * geom


def r_cent(mass=0.5 * u.Msun, omega=1e-14 / u.s, r0=1e4 * u.au):
    """
    Centrifugal radius or disk radius  in the Ulrich (1976)'s model.
    r_u in Mendoza's nomenclature.

    :param mass: Central mass for the protostar
    :param omega: Angular speed at the r0 radius
    :param r0: Initial radius of the streamline
    :return:
    """
    return (r0 ** 4 * omega ** 2 / (G * mass)).to(u.au)


def theta_root(theta, r_to_rc=0.1, theta0=np.radians(30), ecc=1.,
              orb_ang=90*u.deg):
    """
    function to determine theta numerically by finding the root of a function
    This is equation (9) in Mendoza+(2009)

    :param theta: angle in radians
    :param r_to_rc: radius in units of the centrifugal radius
    :param theta0: Initial angle of the streamline
    :param ecc: eccentricity of the orbit (equation 6)
    :param orb_ang: angle in the orbital motion (equation 7)
    :return: returns the difference between the radius and the predicted one,
           a value of 0 corresponds to a proper streamline
    """
    xi = np.arccos(np.cos(theta) / np.cos(theta0)) + orb_ang
    geom = np.sin(theta0)**2 / (1 - ecc * np.cos(xi))
    return r_to_rc - geom


def get_dphi(theta, theta0=np.radians(30)):
    """
    Gets the difference in Phi.

    :param theta: radians
    :param theta0: radians

    :return difference in Phi angle, in radians
    """
    return np.arccos(np.tan(theta0) / np.tan(theta))


def stream_line(r, mass=0.5 * u.Msun, r0=1e4 * u.au, theta0=30 * u.deg,
                omega=1e-14 / u.s, v_r0=0 * u.km / u.s):
    """
    It calculates the stream line following Mendoza et al. (2009)
    It takes the radial velocity and rotation at the streamline
    initial radius and it describes the entire trajectory.

    :param r:
    :param mass:
    :param r0:
    :param theta0:
    :param phi0:
    :param omega:
    :param v_r0: Initial radial velocity
    :return: theta
    """
    # convert theta0 into radians
    rad_theta0 = theta0.to(u.rad).value
    rc = r_cent(mass=mass, omega=omega, r0=r0)
    theta = np.zeros_like(r.value) + np.nan
    # Initial guess at largest radius is tetha0 + epsilon for numerical reasons
    theta_i = rad_theta0 + 1e-3
    # mu and nu are dimensionless
    mu = (rc / r0).decompose().value
    nu = (v_r0 * np.sqrt(rc / (G * mass))).decompose().value
    epsilon = nu**2 + mu**2 * np.sin(theta0)**2 - 2 * mu
    ecc = np.sqrt(1 + epsilon * np.sin(theta0)**2)
    orb_ang = np.arccos((1 - mu * np.sin(theta0)**2) / ecc)
    for ind in np.arange(len(r)):
        r_i = (r[ind] / rc).decompose().value
        if r_i > 1:
            result = optimize.root(theta_root, theta_i,
                                   args=(r_i, rad_theta0, ecc, orb_ang))
            theta_i = result.x
            theta[ind] = theta_i
    return theta * u.rad


def stream_line_vel(r, theta, mass=0.5 * u.Msun, r0=1e4 * u.au, theta0=30 * u.deg,
                omega=1e-14 / u.s, v_r0=0 * u.km / u.s):
    """
    It calculates the stream line following Mendoza et al. (2009)
    It takes the radial velocity and rotation at the streamline
    initial radius and it describes the entire trajectory.

    :param r:
    :param mass:
    :param r0:
    :param theta0:
    :param phi0:
    :param omega:
    :param v_r0: Initial radial velocity
    :return: theta
    """
    # convert theta0 into radians
    rad_theta0 = theta0.to(u.rad).value
    rc = r_cent(mass=mass, omega=omega, r0=r0)
    theta = np.zeros_like(r.value) + np.nan
    # Initial guess at largest radius is tetha0 + epsilon for numerical reasons
    theta_i = rad_theta0 + 1e-3
    # mu and nu are dimensionless
    mu = (rc / r0).decompose().value
    nu = (v_r0 * np.sqrt(rc / (G * mass))).decompose().value
    epsilon = nu**2 + mu**2 * np.sin(theta0)**2 - 2 * mu
    ecc = np.sqrt(1 + epsilon * np.sin(theta0)**2)
    orb_ang = np.arccos((1 - mu * np.sin(theta0)**2) / ecc)
    for ind in np.arange(len(r)):
        r_i = (r[ind] / rc).decompose().value
        if r_i > 1:
            result = optimize.root(theta_root, theta_i,
                                   args=(r_i, rad_theta0, ecc, orb_ang))
            theta_i = result.x
            theta[ind] = theta_i
    return theta * u.rad

def rotate_xyz(x, y, z, inc=30 * u.deg, pa=30 * u.deg):
    """
    Rotate on inclination and PA
    x-axis and y-axis are on the plane on the sky,
    z-axis is the 

    Rotation around x is inclination angle
    Rotation around y is PA angle

    Using example matrices as described in:
    https://en.wikipedia.org/wiki/3D_projection

    :param x: cartesian x-coordinate, in the direction of decreasing RA
    :param y: cartesian y-coordinate, in the direction away of the observer
    :param z: cartesian z-coordinate, in the direction of increasing Dec.
    :param inc: Inclination angle. 0=no change
    :param pa: Change the PA angle. Measured from North due East.
    :return: new x, y, and z-coordinates as observed on the sky, with the
    same units as the input ones.

    """
    xyz = np.column_stack((x, y, z))
    rot_inc = np.array([[1, 0, 0],
                        [0, np.cos(inc), np.sin(inc)],
                        [0, -np.sin(inc), np.cos(inc)]])
    rot_pa = np.array([[np.cos(pa), 0, -np.sin(pa)],
                       [0, 1, 0],
                       [np.sin(pa), 0, np.cos(pa)]])
    xyz_new = rot_pa.dot(rot_inc.dot(xyz.T))
    return xyz_new[0], xyz_new[1], xyz_new[2]


def xyz_stream(mass=0.5 * u.Msun, r0=1e4 * u.au, theta0=30 * u.deg,
               phi0=15 * u.deg, omega=1e-14 / u.s, v_r0=0 * u.km / u.s,
               inc=0 * u.deg, pa=0 * u.deg):
    """
    it gets xyz coordinates for a stream line and it is also rotated 
    in PA and inclination along the line of sight.
    This is a wrapper around stream_line() and rotate_xyz()
    :param mass:
    :param r0:
    :param theta0:
    :param phi0:
    :param omega:
    :param v_r0:
    :param inc:
    :param pa:
    :return:
    """
    rc = r_cent(mass=mass, omega=omega, r0=r0)
    r = np.arange(r0.to(u.au).value, rc.to(u.au).value * 0.99999, step=-10) * u.au
    theta = stream_line(r, mass=mass, r0=r0, theta0=theta0,
                        omega=omega, v_r0=v_r0)
    d_phi = get_dphi(theta, theta0=theta0)
    # Convert from spherical into cartesian coordinates
    z = r * np.cos(theta)
    y = r * np.sin(theta) * np.sin(phi0 + d_phi)
    x = r * np.cos(theta) * np.cos(phi0 + d_phi)
    return rotate_xyz(x, y, z, inc=inc, pa=pa)
