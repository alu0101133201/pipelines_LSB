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


c = np.genfromtxt('build-' + gal + '/output-catalogs/range.txt')
min = c[0]
max = c[1]

# READ the Halpha fluxes of the stars
list = np.genfromtxt('build-' + gal + '/output-catalogs/halpha.cat')
halpha = list[:,3]



# READ the LBTr fluxes of the stars
list = np.genfromtxt('build-' + gal + '/output-catalogs/rc.cat')
lbt_r = list[:,3]


# READ the LBTr fluxes of the stars
list = np.genfromtxt('build-' + gal + '/output-catalogs/gc.cat')
lbt_g = list[:,3]


# Plot to understand
colors = cm.tab20c(np.linspace(0, 20, 73))
fig, (ax,ax1) = plt.subplots(figsize=(16, 12), nrows=1, ncols=2, sharex=False)

halpha[halpha>3000] = np.nan

ax.scatter(lbt_r, halpha ,color = 'purple', s=95, marker='X', label='stars from DATA')
  
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel(' F r-Sloan', fontsize=25)
ax.set_ylabel('F Ha', fontsize=25)



#p = np.polyfit(lbt_r ,halpha,  1)
#ax.plot(lbt_r , lbt_r *p[0] + p[1], color = 'red',alpha=0.9, label="1deg Polyfit " + str(p[0]))
ax.legend(fontsize=25,loc=4)


ax1.scatter(-2.5 * np.log10(lbt_g/lbt_r), halpha/lbt_r ,color = 'purple', s=95, marker='X', label='stars from DATA')
#ax1.hlines(y=np.mean(halpha/lbt_r), xmin=min, xmax=max)
ax1.set_xlim(-0.25,1.5)
ax1.set_xlabel(' g-r', fontsize=25)
ax1.set_ylabel('F Ha/F r Sloan', fontsize=25)




plt.subplots_adjust(hspace=0.0, wspace=0.2)
plt.savefig('build-' + gal + '/plots/f_from_data.png', dpi=300)
plt.cla()
