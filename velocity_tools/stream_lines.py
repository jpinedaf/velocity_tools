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
    geom = np.sqrt(1 - np.cos(theta) / np.cos(theta0)) * 
           np.sin(theta0)/np.sin(theta)
    return v_M( r, M=0.5*u.Msun) * geom


def v_theta(r, theta, M=0.5*u.Msun, theta0=30*u.deg):
    geom = (np.cos(theta0) - np.cos(theta)) * 
           np.sqrt( (np.cos(theta0) + np.cos(theta)) / 
            (np.cos(theta0) * np.sin(theta)**2))
    return v_M( r, M=0.5*u.Msun) * geom

def get_theta(theta, r_to_Rc=0.1, theta0=30*u.deg):
    return r_to_Rc * (1 - np.cos(theta) / np.cos(theta0)) - np.sin(theta0)**2

def stream_line(r, M=0.5*u.Msun, theta0=30*u.deg, 
    Omega=1e-21/u.s, r0=1e4*u.au):
    Rc = (r0**4 * Omega**2/ (G*M)).to(u.au)
    theta = np.zeros_like(r)
    theta_i = theta0
    for ind in np.range(len(r)):
        r_i = (r[i]/Rc).decompose().value
        theta_i = optimize.root( get_theta, theta0=theta_i)
        theta[i] = theta_i
    return theta