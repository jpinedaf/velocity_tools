%run stream_lines.py
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from astropy.constants import G

plt.ion()

r = np.linspace(1e3, 20, num=50) * u.au
theta0_1 = 30*u.deg
theta0_2 = 210*u.deg
Omega0 = 4*(1.7 * u.km/u.s/u.pc).decompose()
r0 = 1e3*u.au
Mstar = 0.5*u.Msun

Omega1 = (np.sqrt(G * Mstar * 150*u.au) / r0**2).decompose()

theta1 = stream_line(r, M=Mstar, theta0=theta0_1, 
    Omega=Omega1, r0=r0)

theta2 = stream_line(r, M=Mstar, theta0=theta0_2, 
    Omega=Omega1, r0=r0)

x1 = r * np.cos(theta1)
x2 = r * np.cos(theta2)
y1 = r * np.sin(theta1)
y2 = r * np.sin(theta2)

plt.close()
plt.plot(x1, y1, color='r', marker='o')
plt.plot(r * np.cos(theta0_1), r * np.sin(theta0_1), color='0.9')
plt.plot(x2, y2, color='b', marker='o')
plt.plot(r * np.cos(theta0_2), r * np.sin(theta0_2), color='0.9')
plt.xlim(-r.max().value, r.max().value)
plt.ylim(-r.max().value, r.max().value)