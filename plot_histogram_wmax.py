# plot_histogram_wmax.py

"""
Code to plot the histogram of maximal weight from the particle filter main code
"""

# Library
import numpy as np
from numpy import *
import math

import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

import pylab as plt
import numpy.matlib
import seaborn as sns


plt.close('all')
plt.rcParams["patch.force_edgecolor"] = True
plt.rcParams.update({'font.size': 30, 'legend.fontsize':30})

p_max = np.loadtxt('output/wmax_onecomponent.txt')

# Plot the histograms for the one-component model
fig = plt.figure(figsize=(15,10))
sns.distplot(p_max[0,:], hist = True, kde=False, kde_kws = {'linestyle': '-'}, color='silver', label = '$N_x=10$', hist_kws=dict(alpha=0.8))
sns.distplot(p_max[1,:], hist = True, kde=False, kde_kws = {'linestyle': '--'}, color='dimgrey', label = '$N_x=30$', hist_kws=dict(alpha=0.8))
sns.distplot(p_max[2, :], hist = True, kde=False, kde_kws = {'linestyle': '-.'}, color='black', label = '$N_x=100$', hist_kws=dict(alpha=0.8))
plt.grid(True)
plt.xlim(0, 1)
plt.ylim(0, 300)
plt.title('One-component')
plt.legend()
plt.xlabel('max $\mathrm{w_i}$')
plt.ylabel('occurrences')

plt.savefig('outputfig/histogram_wmax_onecomponent.png')

p_max = np.loadtxt('output/wmax_multicomponent.txt')
# Plot the histograms for the multi-component model
fig = plt.figure(figsize=(15,10))
sns.distplot(p_max[0,:], hist = True, kde=False, kde_kws = {'linestyle': '-'}, color='silver', label = '$N_x=10$', hist_kws=dict(alpha=0.8))
sns.distplot(p_max[1,:], hist = True, kde=False, kde_kws = {'linestyle': '--'}, color='dimgrey', label = '$N_x=30$', hist_kws=dict(alpha=0.8))
sns.distplot(p_max[2, :], hist = True, kde=False, kde_kws = {'linestyle': '-.'}, color='black', label = '$N_x=100$', hist_kws=dict(alpha=0.8))
plt.grid(True)
plt.xlim(0, 1)
plt.ylim(0, 300)
plt.title('Multi-component')
plt.legend()
plt.xlabel('max $\mathrm{w_i}$')
plt.ylabel('occurrences')

plt.savefig('outputfig/histogram_wmax_multicomponent.png')


