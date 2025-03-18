%run stream_lines.py
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from astropy.constants import G

plt.ion()

theta0_1 = 30*u.deg
theta0_2 = 70*u.deg

Omega0 = 4*(1.7 * u.km/u.s/u.pc).decompose()
r0 = 1e4*u.au
Mstar = 0.5*u.Msun
Omega1 = (np.sqrt(G * Mstar * 300*u.au) / r0**2).decompose()

Rc = R_cent(M=Mstar, Omega=Omega1, r0=r0)
r = np.arange(r0.value, Rc.value, step=-10) * u.au
angles = np.linspace(0, 2 * np.pi, 100)


theta1 = stream_line(r, M=Mstar, theta0=theta0_1, 
    Omega=Omega1, r0=r0)
# theta1_p = stream_line(r, M=Mstar, theta0=theta0_1+10*u.deg,
#     Omega=Omega1, r0=r0)
# theta1_m = stream_line(r, M=Mstar, theta0=theta0_1-10*u.deg, 
#     Omega=Omega1, r0=r0)
theta2 = stream_line(r, M=Mstar, theta0=theta0_2, 
    Omega=Omega1, r0=r0)

dphi2 = dphi(theta2, theta0=theta0_2)
dphi1 = dphi(theta1, theta0=theta0_1)

# theta2_p = stream_line(r, M=Mstar, theta0=theta0_2+10*u.deg,
#     Omega=Omega1, r0=r0)
# theta2_m = stream_line(r, M=Mstar, theta0=theta0_2-10*u.deg, 
#     Omega=Omega1, r0=r0)

phi1 = 15 * u.deg
dphi_f = 15 * u.deg

z1 = r * np.cos(theta1)
y1 = r * np.sin(theta1) * np.sin(phi1 + dphi1)
x1 = r * np.cos(theta1) * np.cos(phi1 + dphi1)

y1_p = r * np.sin(theta1) * np.sin(phi1 + dphi1 + dphi_f)
x1_p = r * np.cos(theta1) * np.cos(phi1 + dphi1 + dphi_f)

y1_m = r * np.sin(theta1) * np.sin(phi1 + dphi1 - dphi_f)
x1_m = r * np.cos(theta1) * np.cos(phi1 + dphi1 - dphi_f)

z2 = r * np.cos(theta2)
y2 = r * np.sin(theta2) * np.sin(phi1 + dphi2)
x2 = r * np.cos(theta2) * np.cos(phi1 + dphi2)

y2_p = r * np.sin(theta2) * np.sin(phi1 + dphi2 + dphi_f)
x2_p = r * np.cos(theta2) * np.cos(phi1 + dphi2 + dphi_f)

y2_m = r * np.sin(theta2) * np.sin(phi1 + dphi2 - dphi_f)
x2_m = r * np.cos(theta2) * np.cos(phi1 + dphi2 - dphi_f)


# x1_m = r * np.cos(theta1_m)
# x1_p = r * np.cos(theta1_p)
# x2 = r * np.cos(theta2)
# x2_m = r * np.cos(theta2_m)
# x2_p = r * np.cos(theta2_p)
# xl_1 = r * np.cos(theta1[-1])
xl_2 = r * np.cos(theta2[-1])
#
# y1 = r * np.sin(theta1)
# y1_m = r * np.sin(theta1_m)
# y1_p = r * np.sin(theta1_p)
# y2 = r * np.sin(theta2)
# y2_m = r * np.sin(theta2_m)
# y2_p = r * np.sin(theta2_p)
# yl_1 = r * np.sin(theta1[-1])
# yl_2 = r * np.sin(theta2[-1])

plt.close()
plt.plot(x1, y1, color='r', marker='o', markersize=3)
plt.plot(x1_m, y1_m, color='r', marker='o', markersize=3)
plt.plot(x1_p, y1_p, color='r', marker='o', markersize=3)
# plt.plot(xl_1, yl_1, color='0.9')
plt.plot(x2, y2, color='b', marker='o', markersize=3)
plt.plot(x2_m, y2_m, color='b', marker='o', markersize=3)
plt.plot(x2_p, y2_p, color='b', marker='o', markersize=3)
# plt.plot(xl_2, yl_2, color='0.9')
plt.plot(Rc * np.cos(angles), Rc * np.sin(angles), color='r')
plt.axis('equal')
plt.xlim(-r.max().value, r.max().value)
plt.ylim(-r.max().value, r.max().value)
plt.tight_layout()
