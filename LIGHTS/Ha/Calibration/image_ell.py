# This script is made to build a 4x4 figure for profiles
# using labels and cropping images at the same physical scale
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



colors = cm.tab20(np.linspace(0, 1, 23))
# for c in colors , color=c

name = str(sys.argv[1])




# only ellipses
fig1, ax = plt.subplots(figsize=(16, 12), nrows=1, ncols=1, sharex=True)


for k in range(0, 1):
    ax.set_xlim(0, 400)
    #ax[j, k].set_xticks([0, 50 ,100,])
    
    
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(axis='x',which='both', direction="in")
    ax.tick_params(axis='x',which='major', length=7)
    ax.tick_params(axis='x',which='minor', length=4)
    
    ax.tick_params(axis='y',which='both', direction="in")
    ax.tick_params(axis='y',which='major', length=7)
    ax.tick_params(axis='y',which='minor', length=4)

    ax.set_ylim(17.5, 32.5)
    ax.invert_yaxis()
    ax.set_yticks([ 20, 22, 24, 26, 28, 30, 32])
    ax.set_yticklabels(['20', '22', '24', '26', '28', '30', '32'],fontsize=15)
    ax.set_xticks([0, 100, 200, 300, 400])
    ax.set_xticklabels(['0','100', '200','300', '400'],fontsize=15)
        
        
        


# Common things on the PATH
path_output = os.getcwd()
path_1 = os.getcwd()
path_2 = '/build-' + name + '/'
path_images = path_1 + path_2

# OPEN THE ZP
zp = 22.5

# read the data


###################################################################################



    
#R
lbt = np.genfromtxt(path_images + 'LBT/prof_lights_' + name + '_Sloan-r.txt')
lbt_pix = lbt[:,0]
lbt_counts = lbt[:,1]
lbt_err = np.sqrt(lbt[:,6])
lbt_counts[lbt_counts<0] = np.nan
# Compute distance and SB
lbt_as_r=lbt_pix * 0.333
lbt_sb_r=-2.5*np.log10(lbt_counts)+ zp + 5.*np.log10(0.333)
lbt_err_sb_r = np.abs((-2.5 * np.log10(lbt_counts + lbt_err) + 5 * np.log10(0.333) + zp) - lbt_sb_r )



###################################################################################



#G
lbt = np.genfromtxt(path_images + 'LBT/prof_lights_' + name + '_Sloan-g.txt')
lbt_pix = lbt[:,0]
lbt_counts = lbt[:,1]
lbt_counts[lbt_counts<0] = np.nan
lbt_err = np.sqrt(lbt[:,6])
# Compute distance and SB
lbt_as_g=lbt_pix * 0.333
lbt_sb_g=-2.5*np.log10(lbt_counts)+ zp + 5.*np.log10(0.333)
lbt_err_sb_g = np.abs((-2.5 * np.log10(lbt_counts + lbt_err) + 5 * np.log10(0.333) + zp) - lbt_sb_g )



###################################################################################


# Ha
lbt = np.genfromtxt(path_images + 'LBT/prof_lights_' + name + '_Sloan-ha.txt')
lbt_pix = lbt[:,0]
lbt_counts = lbt[:,1]
lbt_counts[lbt_counts<0] = np.nan
lbt_err = np.sqrt(lbt[:,6])
# Compute distance and SB
lbt_as_ha=lbt_pix * 0.333
lbt_sb_ha=-2.5*np.log10(lbt_counts)+ zp + 5.*np.log10(0.333)
lbt_err_sb_ha = np.abs((-2.5 * np.log10(lbt_counts + lbt_err) + 5 * np.log10(0.333) + zp) - lbt_sb_ha )




        
ax.fill_between(lbt_as_r, lbt_sb_r+lbt_err_sb_r ,lbt_sb_r-lbt_err_sb_r, color="orange",label="r Sloan" , alpha=0.9)
    
ax.fill_between(lbt_as_g, lbt_sb_g+lbt_err_sb_g ,lbt_sb_g-lbt_err_sb_g, color="blue",label="g Sloan" , alpha=0.9)

#ax.fill_between(lbt_as_ha, lbt_sb_ha+lbt_err_sb_ha +0.065,lbt_sb_ha-lbt_err_sb_ha+0.065, color="purple",label="Ha INT extcorr" , alpha=0.9)
ax.fill_between(lbt_as_ha , lbt_sb_ha+lbt_err_sb_ha ,lbt_sb_ha-lbt_err_sb_ha, color="red",label="Ha INT" , alpha=0.9)

#ax.fill_between(lbt_as_ha_sub, lbt_sb_ha_sub + lbt_err_sb_ha_sub ,lbt_sb_ha_sub -lbt_err_sb_ha_sub, color="green",label="Ha INT sub" , alpha=0.9)
    
      
ax.set_ylabel('$\mu_{\mathit{}}$ (mag/arcsec$^2$)', fontsize=25)
ax.set_xlabel('$\mathit{R}$ (arcsec)', fontsize=25)
ax.legend(fontsize=25, loc=1)


plt.subplots_adjust(hspace=0.0, wspace=0.2)
plt.savefig('build-'+ name + '/data_ellipses_' + name + '.png')
plt.cla()

plt.close()
