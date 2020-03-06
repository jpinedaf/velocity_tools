%run stream_lines.py
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt

# Initial test: change in PA for vector = showin x- and y-axes
x_x = 1.
y_x = 0.
z_x = 0.
x_y = 0.
y_y = 1.
z_y = 0.

x_x_new, y_x_new, z_x_new = rotate_xyz(x_x, y_x, z_x, inc=0*u.deg, PA=45*u.deg)
x_y_new, y_y_new, z_y_new = rotate_xyz(x_y, y_y, z_y, inc=0*u.deg, PA=45*u.deg)

# Plot x- and y-axes only
plt.ion()
plt.close()
fig1, ax1 = plt.subplots(figsize=(7,7))
ax1.plot([0, x_x], [0, y_x], color='red')
ax1.plot([0, x_y], [0, y_y], color='red')
ax1.plot([0, x_x_new], [0, y_x_new], color='blue')
ax1.plot([0, x_y_new], [0, y_y_new], color='blue')
ax1.axis('equal')
ax1.set_xlabel('x')
ax1.set_ylabel('y')


# test: change in inclination for vector = showin x- and y-axes
x_x = 1.
y_x = 0.
z_x = 0.
x_y = 0.
y_y = 1.
z_y = 0.

ang_c = np.linspace(0, 2*np.pi, 100)
x_circ = np.cos(ang_c)
y_circ = np.sin(ang_c)
z_circ = np.cos(ang_c) * 0.0
x_circ_new, y_circ_new, z_circ_new = rotate_xyz(x_circ, y_circ, z_circ, inc=60*u.deg, PA=30*u.deg)
# x_y_new, y_y_new, z_y_new = rotate_xyz(x_y, y_y, z_y, inc=60*u.deg, PA=0*u.deg)

# Plot x- and y-axes only
fig2, ax2 = plt.subplots(figsize=(14,7), ncols=2)
ax2[0].plot([0, x_x], [0, y_x], color='red')
ax2[0].plot([0, x_y], [0, y_y], color='red')
ax2[0].plot( x_circ, y_circ, color='red')
ax2[0].plot( x_circ_new, y_circ_new, color='blue', ls=':')
ax2[0].axis('equal')
ax2[0].set_xlabel('x')
ax2[0].set_ylabel('y')

ax2[1].plot([0, z_x], [0, y_x], color='red')
ax2[1].plot([0, x_y], [0, y_y], color='red')
ax2[1].plot( z_circ, y_circ, color='red')
ax2[1].plot( z_circ_new, y_circ_new, color='blue', ls=':')
ax2[1].axis('equal')
ax2[1].set_xlabel('z')
ax2[1].set_ylabel('y')
