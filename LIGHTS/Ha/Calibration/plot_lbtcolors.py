# This script is made plot SDSS MAG - GAIA g mag relation
#
#
# System imports
import os
import sys
import pdb
import warnings

import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import matplotlib.ticker as ticker
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import matplotlib.cm as cm

import os, glob, sys, warnings, array, re, math, time, copy
import numpy               as np
import matplotlib.pyplot   as plt
from   astropy.io          import fits, ascii
import gc
import matplotlib
import matplotlib.ticker as ticker
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import matplotlib.cm as cm
import scipy
from scipy.interpolate import interp1d
from astropy.table import Table
from scipy.stats import norm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt





gal = str(sys.argv[1])
grmax = float(sys.argv[2])
mulsigma = float(sys.argv[3])





# READ the colors of the stars
list = np.genfromtxt('build-' + gal + '/output-catalogs/hc_.cat')
g_r = list[:,3]
data = np.arange(0, len(g_r), 1)



# Plot to understand color distribution
fig, ax = plt.subplots(figsize=(16, 12), nrows=1, ncols=1, sharex=False)


ax.scatter(data, g_r ,color = 'red', s=95, marker='X',label= 'LBT stars')

#ax.hlines(xmin=0,xmax=len(g_r), y=0.)
ax.set_xlabel('Id star', fontsize=25)
ax.set_ylabel('g-r', fontsize=25)


ax.legend(fontsize=25,loc=4)
#ax1.legend(fontsize=25,loc=4)

plt.subplots_adjust(hspace=0.0, wspace=0.2)
plt.savefig('build-' + gal + '/plots/lbt_colors.png', dpi=300)
plt.cla()


# Correct range of colors in my Field to save factor
fig1, (ax1,ax2) = plt.subplots(figsize=(16, 12), nrows=1, ncols=2, sharex=False)

# CUT redder colors
good_colors = []
for i in np.arange(0, len(g_r), 1):
    if g_r[i] < grmax:
        good_colors.append(g_r[i])

# the histogram of the data
n, bins, patches = ax1.hist(good_colors, 50, facecolor='green', alpha=0.75, label= gal + ' field stars')
print(n)        # frequency
m = np.max(n)

# save maximum value
a = None
for i in np.arange(0, len(n), 1):
    if n[i] == m:
        a = bins[i]
        print(bins[i])     # values

#save approximative mu and sigma/2 to select the range where computing the factor
mu, sigma = norm.fit(good_colors)

min = a - sigma/mulsigma
max = a + sigma/mulsigma



# Save into a file the min and max val of colors where computing the factor
file = open('build-' + gal + '/output-catalogs/range.txt', "w")
file.write(str(min))
file.write(" ")
file.write(str(max))


# show maximum
ax1.vlines(ymin=0, ymax=m, x=a, alpha=1, lw=5, color = 'black')

# show range
ax1.vlines(ymin=0, ymax=m-1, x=min, alpha=1, lw=5, color = 'red')
ax1.vlines(ymin=0, ymax=m-1, x=max, alpha=1, lw=5, color = 'red')


ax1.set_xlabel('g - r', fontsize=25)
ax1.set_ylabel('n stars', fontsize=25)


# same but all the group on the right
ax2.hist( g_r, 50, label= gal + 'field stars')
ax2.set_xlabel('g - r', fontsize=25)
ax2.set_ylabel('n stars', fontsize=25)


plt.savefig('build-' + gal + '/plots/lbt_colors_hist.png', dpi=300)
plt.cla()
