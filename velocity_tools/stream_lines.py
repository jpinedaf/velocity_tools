import numpy as np
from astropy.io import fits
import astropy.units as u
from scipy import optimize
from astropy.constants import G

def v_M( r, M=0.5*u.Msun):
    """ Velocity term that is repeated in all velocity component
    """
    return np.sqrt(G * M/r).to(u.km/u.s)


def v_r( r, theta, M=0.5*u.Msun, theta0=30*u.deg):
    """Radial velocity component
    """
    return v_M( r, M=0.5*u.Msun) * np.sqrt(1 + np.cos(theta) / np.cos(theta0))


def v_phi(r, theta, M=0.5*u.Msun, theta0=30*u.deg):
    geom = np.sqrt(1 - np.cos(theta) / np.cos(theta0)) *\
           np.sin(theta0)/np.sin(theta)
    return v_M( r, M=0.5*u.Msun) * geom


def v_theta(r, theta, M=0.5*u.Msun, theta0=30*u.deg):
    geom = (np.cos(theta0) - np.cos(theta)) *\
           np.sqrt( (np.cos(theta0) + np.cos(theta)) /\
            (np.cos(theta0) * np.sin(theta)**2))
    return v_M( r, M=0.5*u.Msun) * geom

def R_cent(M=0.5*u.Msun, Omega=1e-14/u.s, r0=1e4*u.au):
    return (r0**4 * Omega**2/ (G*M)).to(u.au)

def get_theta(theta, r_to_Rc=0.1, theta0=np.radians(30)):
    geom = np.sin(theta0)**2 / (1 - (np.cos(theta) / np.cos(theta0)))
    return r_to_Rc - geom

def dphi(theta, theta0=np.radians(30)):
    """
    Gets the difference in Phi.
    theta abd theta0 are in radians
    """
    return np.arctan( np.tan(np.arccos(np.cos(theta) / np.cos(theta0))) / np.sin(theta0))

def stream_line(r, M=0.5*u.Msun, theta0=30*u.deg, 
    Omega=1e-14/u.s, r0=1e4*u.au):
    # convert theta0 into radians
    rad_theta0 = theta0.to(u.rad).value
    Rc = R_cent(M=M, Omega=Omega, r0=r0)
    theta = np.zeros_like(r.value) + np.nan
    # Initial guess at largest radius is tetha0 + epsilon for numerical reasons
    theta_i = rad_theta0 + 1e-3
    for ind in np.arange(len(r)):
        r_i = (r[ind] / Rc).decompose().value
        result = optimize.root( get_theta, theta_i, args=(r_i, rad_theta0))
        theta_i = result.x
        theta[ind] = theta_i
        # print('r_to_Rc={0}   theta0={1}  and sol={2}'.format(r_i, rad_theta0, theta[ind]))
    bad = (r < Rc)
    theta[bad] = np.nan
    return theta * u.rad

def rotate_xyz(x, y, z, inc=30*u.deg, PA=30*u.deg):
    """
    Rotate on inclination and PA
    x-axis and y-axis are on the plane on the sky,
    z-axis is the 

    Rotation around x is inclination angle
    Rotation around y is PA angle

    Using example matrices as decribed in:
    https://en.wikipedia.org/wiki/3D_projection
    """
    xyz = np.column_stack((x, y, z))
    Rot_inc = np.array([[1, 0, 0],
                       [0, np.cos(inc), np.sin(inc)],
                       [0, -np.sin(inc), np.cos(inc)]])
    Rot_PA = np.array([[np.cos(PA), 0, -np.sin(PA)],
                       [0, 1, 0],
                       [np.sin(PA), 0, np.cos(PA)]])
    xyz_new = Rot_PA.dot(Rot_inc.dot(xyz.T))
    return xyz_new[0], xyz_new[1], xyz_new[2]