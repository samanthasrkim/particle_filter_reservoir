# PF_main_multicomponent.py
"""
Code to compute the maximal value of weights using the reservoir model

Ne = ensemble size
Nx = state space dimension
Ny = observation space dimension

"""

# Library
import numpy as np
from numpy import *
import math

import matplotlib.mlab as mlab
from random import *
import statistics as stat

from matplotlib import cm

import pylab as plt
import numpy.matlib
import scipy.integrate as integrate

from numpy.random import randn
from numpy.linalg import inv

from functions_for_forward_model import forward

# Set the experiment
n_run = 1000
Ne = 100
nu = 0.32
print('Ensemble size Ne=', Ne)

# Set the dimension
dim = np.array([10, 30, 100])
Ny = dim

# Prepare the grid interval for the reservoir
interval = 0.05
dx = interval
dy = interval
dz = 240
depth = 2800

# in degree 
grid_degree = 1
x_center = 6.77
y_center = 53.31
area_size = 0.5
x0_area = x_center - area_size/2
y0_area = y_center - area_size/2
xmax_area = x_center + area_size/2
ymax_area = y_center + area_size/2
step = 0.005
n_size = int((xmax_area - x0_area)/step)

x = np.linspace(x0_area, xmax_area, n_size)
y = np.linspace(y0_area, ymax_area, n_size)
X_grid, Y_grid = np.meshgrid(x, y)
# Prior distribution
sigma_en = 0.1 # Strength
var_en = sigma_en**2
prior_mean = 0

# Observational noise
sigma_eps =  1e-3 # In meter 
var_eps = sigma_eps**2
eps_mean = 0

# more realistic values
sigma_lik = 1e-3 # in meter
mean_lik = 0

p_max = np.zeros((3, n_run))

## Load the position of the nuclei of strain given the dimension
for dim_index in range(3):
    Ny = dim[dim_index]
    print('Dimension:', Ny)
    if Ny == 10:
        mogi = loadtxt('inputs_model/mogi10dx00222.txt')
        x_truth = np.random.normal(prior_mean, sigma_en, Ny)
    elif Ny == 30:
        mogi = loadtxt('inputs_model/mogi30dx00222.txt')
        x_truth = np.random.normal(prior_mean, sigma_en, Ny)
    elif Ny == 100:
        mogi = loadtxt('inputs_model/mogi100dx00222.txt')
        x_truth = np.random.normal(prior_mean, sigma_en, Ny)

    # Initialize array
    data = np.zeros(Ny)
    x_p = np.zeros((Ne, Ny))

    hgrid_true = np.zeros((n_size, n_size))
    
    xdis_true = np.zeros((n_size, n_size))
    ydis_true = np.zeros((n_size, n_size))
    Z_grid = np.zeros((n_size, n_size))
      
    # Set the position of the data
    data = np.zeros((Ny, 3))
    for i in range(Ny):
        data[i, 0] = mogi[i, 0]
        data[i, 1] = mogi[i, 1]
    mogi[:, 2]= depth
    
    # Compute the true deformation at observation point
    z_true_def = np.zeros((1, Ny))
    X_obs = data[:, 0]
    Y_obs = data[:, 1]
    Z_obs = data[:, 2]*0
    
    for i in range(Ny):
        u = forward(X_obs, Y_obs, Z_obs, mogi[i, :], x_truth[i], nu, dx, dy, dz, grid_degree)
        z_true_def += u[2]

    for simu in range(n_run):
        #print('run =', simu)
        for i in range(Ny):
            # Prior
            x_p[:, i] = np.random.normal(prior_mean, sigma_en, Ne)
        # Observations
        data[:, 2] = z_true_def + np.random.normal(eps_mean, sigma_eps, Ny)
    
        """ 
        Filtering part
        """
        # initilization of the weight array
        p_w = np.ones(Ne, dtype=np.float64)
        # dlpsi is squared difference between data and model
        dlpsi = np.zeros((Ne, Ny))#, dtype=np.float64)
        dlpsi_sq = np.zeros((Ne, Ny))#, dtype=np.float64)
        z_p_def = np.zeros((Ny, Ne))     

        # Compute the deformation at observation point given the values of the particles
        for i in range(Ne):
            for j_source in range(Ny):
                u = forward(X_obs, Y_obs, Z_obs, mogi[j_source, :], x_p[i, j_source], nu, dx, dy, dz, grid_degree)
                uz = u[2]
                z_p_def[:, i] = np.add(z_p_def[:, i], uz)

        z_p_transpose  = z_p_def.T

        for i in range(Ne):
        # Calculation of the diff data/model for each particles for each observation points
            dlpsi[i, :] = (data[:, 2] - z_p_transpose[i,:])**2
        
        # Calculation of the observational probability density function (likelihood)
        # temporairement
        for n in range(Ne):
            dlpsi[n, :] = (1/sqrt(2*pi*sigma_lik))*np.exp(-0.5*dlpsi[n, :]/(sigma_lik**2))

        for i in range(Ne):
            for j in range(Ny):
                p_w[i] = p_w[i] * dlpsi[i,j]
        
        # Normalization
        p_w = p_w / sum(p_w)

        # store the value of the p_w max before to run the code again
        p_max[dim_index, simu] = max(p_w)

# Array with values of wmax
savetxt('output/wmax_multicomponent.txt', p_max)
