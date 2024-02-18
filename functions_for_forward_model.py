# Functions for the forward model of the Mogi source
# Adapted from:
# Scott Henderson 8/31/2012
# https://github.com/scottyhq/cov9/blob/master/mogi.py

import numpy as np


def cart2pol(x1, x2):
    theta = np.arctan2(x2, x1)
    r = np.hypot(x2, x1)
    return theta, r


def pol2cart(theta, r):
    x1 = r * np.cos(theta)
    x2 = r * np.sin(theta)
    return x1, x2

def forward(x, y, z, mogi, mogi_strength, nu, dx, dy, dz, grid_degree):
    # input
    # mogi: 1x3 array [degree] for xy and [m] for z
    pi = np.pi
    # mean radius of the earth [m]
    earth_rad = 6371000
    # earth circumference
    earth_cir = 2 * pi * earth_rad

    # 1) Calculation of coordinates differences between the source and the grid
    # x and y in [degree] and z in [m]

    x = x - mogi[0]
    y = y - mogi[1]
    z = z - mogi[2]


    if (grid_degree == 1):
        x = x * (earth_cir / 360)  # [m]
        y = y * (earth_cir / 360)  # [m]
        dx = dx * (earth_cir / 360)  # [m]
        dy = dy * (earth_cir / 360)  # [m]
        #dz = dz * (earth_cir / 360)  # [m]
    # 2) Calculation of the RADIALE  (cf. in spherical coordinates) distance between the source and the grid
    # We use cylindric coordinates: convert to cylindric coordinates
    th, rho = cart2pol(x, y)
    R = np.hypot(z, rho)

    # 3) Mogi displacement calculation
    # Calculation of the strength if we have dV as input
    V = dx*dy*dz
    # strength proportional to Cm*dP
    dp_initial = 1e6 #Pa 10e6
    est_comp = 1e-10 #Pa-1
    #est_comp = 1e-9
    mogi_strength = mogi_strength*dp_initial*est_comp
    
    C = mogi_strength * ((1 - nu) / np.pi) * V
    ur = C * rho / R ** 3
    uz = C * mogi[2] / R ** 3

    # 4) Convert in cartesian coordinates
    ux, uy = pol2cart(th, ur)
    u = np.array([ux, uy, uz])

    return u 



