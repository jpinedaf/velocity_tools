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


def r_cent(mass=0.5 * u.Msun, omega=1e-14 / u.s, r0=1e4 * u.au):
    """
    Centrifugal radius or disk radius in the Ulrich (1976)'s model.
    r_u in Mendoza's nomenclature.

    :param mass: Central mass for the protostar
    :param omega: Angular speed at the r0 radius
    :param r0: Initial radius of the streamline
    :return:
    """
    return (r0 ** 4 * omega ** 2 / (G * mass)).to(u.au)


def theta_abs(theta, r_to_rc=0.1, theta0=np.radians(30), ecc=1.,
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
    cos_ratio = np.cos(theta) / np.cos(theta0)
    if cos_ratio > 1.:
        print('theta0={0}, theta_try={1} --> bad arccos calculation'.format(theta0, theta))
        return np.nan
    xi = np.arccos(cos_ratio) + orb_ang.to(u.rad).value
    geom = np.sin(theta0)**2 / (1 - ecc * np.cos(xi))
    return np.abs(r_to_rc - geom)


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
    print("rc={0}".format(rc))
    theta = np.zeros_like(r.value) + np.nan
    # mu and nu are dimensionless
    mu = (rc / r0).decompose().value
    nu = (v_r0 * np.sqrt(rc / (G * mass))).decompose().value
    # epsilon is the dimensionless energy
    epsilon = nu**2 + mu**2 * np.sin(theta0)**2 - 2 * mu
    ecc = np.sqrt(1 + epsilon * np.sin(theta0)**2)
    orb_ang = np.arccos((1 - mu * np.sin(theta0)**2) / ecc)
    # the first element in the streamline is the starting point
    theta[0] = rad_theta0
    # Initial guess at largest radius is theta0 +- initguess towards the midplane
    deltar = np.amin(np.abs(np.roll(r,1) - r))
    print(deltar)
    # we use a constant of 6e-5 for an epsilon of 0.01 km/s
    # this result will be in radians
    tol = (6.e-5 * deltar * omega / (v_r0+ 0.1 * u.km/u.s)).decompose().value
    # the initial guess will be 10 times the tolerance for now, in testing
    initguess = 10 * tol
    print('tolerance ', tol)
    if rad_theta0 < np.radians(90):
        theta_i = rad_theta0 + initguess
        theta_bracket = [(rad_theta0, np.pi/2.)]
    else:
        theta_i = rad_theta0 - initguess
        theta_bracket = [(np.pi/2., rad_theta0)]
    for ind in np.arange(1, len(r)):
        r_i = (r[ind] / rc).decompose().value
        if r_i > 0.5:
            # print('initial guess of theta_i = {0}'.format(theta_i))
            # result = optimize.minimize(theta_abs, theta_i,
            #                            bounds=theta_bracket,
            #                            args=(r_i, rad_theta0, ecc, orb_ang))
            # By default, when minimize receives bounds and no constrains,
            # it uses the L-BFGS-B method:
            # ftol is the tolerance in the function evaluation
            # "The iteration stops when (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= ftol"
            # gtol corresponds to the parameter pgtol in fmin_l_bfgs_b
            # "The iteration will stop when max{|proj g_i | i = 1, ..., n} <= gtol"
            # eps corresponds to the absolute step size used for numerical approximation of the jacobian via forward differences.
            options_dict = {'gtol': tol/10., 'eps': tol, 'ftol': tol}
            result = optimize.minimize(theta_abs, theta_i,
                                       bounds=theta_bracket,
                                       args=(r_i, rad_theta0, ecc, orb_ang),
                                       options=options_dict)
            theta_i = result.x
            # These prints are to diagnose if the minimization is converging
            # print(ind, result.success)
            # print(result.message, result.status, result.nit)
            theta[ind] = theta_i
    return theta * u.rad


def stream_line_vel(r, theta, mass=0.5*u.Msun, r0=1e4*u.au, theta0=30*u.deg,
                omega=1e-14/u.s, v_r0=0*u.km/u.s):
    """
    It calculates the velocity along the stream line following Mendoza+(2009)
    It takes the radial velocity and rotation at the streamline
    initial radius and it describes the entire trajectory.

    :param theta:
    :param r:
    :param mass:
    :param r0:
    :param theta0:
    :param phi0:
    :param omega:
    :param v_r0: Initial radial velocity
    :return: v_r, v_theta, v_phi in units of km/s
    """
    rc = r_cent(mass=mass, omega=omega, r0=r0)
    r_to_rc = (r / rc).decompose().value
    v_k0 = v_k(rc, mass=mass)
    # mu and nu are dimensionless
    mu = (rc / r0).decompose().value
    nu = (v_r0 * np.sqrt(rc / (G * mass))).decompose().value
    epsilon = nu**2 + mu**2 * np.sin(theta0)**2 - 2 * mu
    ecc = np.sqrt(1 + epsilon*np.sin(theta0)**2)
    orb_ang = np.arccos((1 - mu * np.sin(theta0)**2) / ecc)
    cos_ratio = np.cos(theta) / np.cos(theta0)
    xi = np.arccos(cos_ratio) + orb_ang.to(u.rad)#.value
    #
    v_r_all = -ecc * np.sin(theta0) * np.sin(xi) / r_to_rc /(1 - ecc*np.cos(xi))
    v_theta_all = np.sin(theta0) / np.sin(theta) / r_to_rc \
                  * np.sqrt(np.cos(theta0)**2 - np.cos(theta)**2)
    v_phi_all = np.sin(theta0)**2 / np.sin(theta) / r_to_rc

    return v_r_all * v_k0, v_theta_all * v_k0, v_phi_all * v_k0


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


def xyz_stream(mass=0.5*u.Msun, r0=1e4*u.au, theta0=30*u.deg,
               phi0=15*u.deg, omega=1e-14/u.s, v_r0=0*u.km/u.s,
               inc=0*u.deg, pa=0*u.deg, rmin=None, deltar=1*u.au):
    """
    it gets xyz coordinates and velocities for a stream line.
    They are also rotated in PA and inclination along the line of sight.
    This is a wrapper around stream_line() and rotate_xyz()

    Spherical into cartesian transformation is done for position and velocity
    using:
    https://en.wikipedia.org/wiki/Vector_fields_in_cylindrical_and_spherical_coordinates

    :param mass: Central mass
    :param r0: Initial radius of streamline
    :param theta0: Initial polar angle of streamline
    :param phi0: Initial azimuthal angle of streamline
    :param omega: Angular rotation. (defined positive)
    :param v_r0: Initial radial velocity of the streamline
    :param inc: inclination with respect of line-of-sight, inc=0 is an edge-on-disk
    :param pa: Position angle of the rotation axis, measured due East from North. This is usually estimated from the outflow PA, or the disk PA-90deg.
    :param rmin: smallest radius for calculation
    :param deltar: spacing between two consecutive radii in the sampling of the streamer, in au
    :return:
    """
    rc = r_cent(mass=mass, omega=omega, r0=r0)
    if rc > r0:
        print('Centrifugal radius is larger than start of streamline')
    r = np.arange(r0.to(u.au).value, rc.to(u.au).value*0.5, step=-1*deltar.value) * u.au
    theta = stream_line(r, mass=mass, r0=r0, theta0=theta0,
                        omega=omega, v_r0=v_r0)
    d_phi = get_dphi(theta, theta0=theta0)
    phi = phi0 + d_phi
    #
    v_r, v_theta, v_phi = stream_line_vel(r, theta, mass=mass, r0=r0,
                                          theta0=theta0, omega=omega, v_r0=v_r0)
    v_x = v_r * np.sin(theta) * np.cos(phi) \
          + v_theta * np.cos(theta) * np.cos(phi) \
          - v_phi * np.sin(phi)
    v_y = v_r * np.sin(theta) * np.sin(phi) \
          + v_theta * np.cos(theta) * np.sin(phi) \
          + v_phi * np.cos(phi)
    v_z = v_r * np.cos(theta) \
          - v_theta * np.sin(theta)
    # Convert from spherical into cartesian coordinates
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    if rmin is not None:
        gd_rmin = (r > rmin)
        if gd_rmin.sum() > 0:
            return rotate_xyz(x[gd_rmin], y[gd_rmin], z[gd_rmin], inc=inc, pa=pa),\
                rotate_xyz(v_x[gd_rmin], v_y[gd_rmin], v_z[gd_rmin], inc=inc, pa=pa)
        else:
            return [np.nan], [np.nan], [np.nan], [np.nan], [np.nan], [np.nan]
    else:
        return rotate_xyz(x, y, z, inc=inc, pa=pa), \
               rotate_xyz(v_x, v_y, v_z, inc=inc, pa=pa)
