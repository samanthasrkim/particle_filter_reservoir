# PF_main_onecomponent.py
"""
Particle filter code for the one-component model

Ne = ensemble size
Nx = state space dimension
Ny = observation space dimension

"""

# Library
import numpy as np
from numpy import *
import math

import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from random import *
import statistics as stat

from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

import pylab as plt
import numpy.matlib
import scipy.integrate as integrate
import scipy as sp
from scipy.interpolate import make_interp_spline
from scipy.interpolate import interp1d
import seaborn as sns


from numpy.random import randn
from numpy.linalg import inv


plt.close('all')
plt.rcParams["patch.force_edgecolor"] = True
plt.rcParams.update({'font.size': 30, 'legend.fontsize':30})

# number of time we run the experiment for simulations
n_run = 1000
Ne = 100
print('Ensemble size Ne=', Ne)

# Set the dimension
dim = np.array([10, 30, 100])
Ny = dim

# Prior distribution
var_en = 1 
sigma_en = sqrt(var_en)
prior_mean = 0

# Observational noise
var_eps = 1
sigma_eps = sqrt(var_eps)
eps_mean = 0

p_max = np.zeros((3, n_run))
for dim_index in range(3):
    Ny = dim[dim_index]
    print('Dimension Ny = ', Ny)

    # Initialize array
    data = np.zeros(Ny)
    x_p = np.zeros((Ne, Ny))
      
    # Set true value of states
    x_truth = np.random.normal(prior_mean, sigma_en, Ny)

    for simu in range(n_run):
        #print('run =', simu)
        for i in range(Ny):
            # Prior
            x_p[:, i] = np.random.normal(prior_mean, sigma_en, Ne)
            # Observations
            data = x_truth + np.random.normal(eps_mean, sigma_eps, Ny)
    
        # Likelihood variance
        sigma_lik = 1
        mean_lik = 0
    
        """ 
        Filtering part
        """
        # initilization of the weight array
        p_w = np.ones((Ne, 1), dtype=np.float64)       
    
        # dlpsi is squared difference between data and model
        dlpsi = np.zeros((Ne, Ny), dtype=np.float64)

        for i in range(Ne):
        # Calculation of the diff data/model for each particles for each observation points
            dlpsi[i, :] = (data - x_p[i,:].T)**2
        
            # Calculation of the observational probability density function (likelihood)
        for n in range(Ne):
            dlpsi[n, :] = (1/sqrt(2*pi*sigma_lik))*np.exp(-0.5*dlpsi[n, :]/(sigma_lik**2))

        # Calculation of the weights
        # test normalization: p_w normalized is the same with or without the factor in dlpsi
        dlpsi_fact = dlpsi *6
        for i in range(Ne):
            for j in range(Ny):
                p_w[i] = p_w[i] * dlpsi_fact[i,j]
        # Normalization
        p_w = p_w / sum(p_w)
    

        # store the value of the p_w max before to run the code again
        p_max[dim_index, simu] = max(p_w)


# Array with values of wmax
savetxt('output/wmax_onecomponent.txt', p_max)
#print("wmax Nx10", np.mean(p_max[0,:]))
#print("wmax Nx30", np.mean(p_max[1,:]))
#print("wmax Nx50", np.mean(p_max[2,:]))



