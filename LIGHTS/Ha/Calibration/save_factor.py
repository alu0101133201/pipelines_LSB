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
from astropy.io import fits


# read filters
def read_filter(filter_file):
     data     = ascii.read(filter_file)
     wavefilt = data['col1']
     fluxfilt = data['col2']
     fluxfilt = fluxfilt / np.amax(fluxfilt)
     return wavefilt, fluxfilt






# Convolution to get magnitudes
def compute_mags(wavemod, fluxmod, filter_files):
     cvel      = 2.99792458e18 # Speed of light in Angstron/sec
     dl        = 1E-5          # 10 pc in Mpc, z=0; for absolute magnitudes
     cfact     = 5.0 * np.log10(1.7684E8 * dl) # from lum[erg/s/A] to flux [erg/s/A/cm2]
     nfilters  = len(filter_files)
     outmag    = np.zeros(nfilters)
     wave_eff  = np.zeros(nfilters)
     veg       = np.zeros(nfilters)

    # Computing the magnitude for each filter
     for ii in np.arange(0, nfilters):
         wavefilt, fluxfilt = read_filter(filter_files[ii])
         wave_eff[ii]       = np.sqrt(np.trapz(fluxfilt, x = wavefilt)/np.trapz(fluxfilt/np.power(wavefilt, 2), x = wavefilt))
   # Finding the wavelength limits of the filters
         good = (wavefilt > 0.0)
         wlow = np.amin(wavefilt[good])
         whi  = np.amax(wavefilt[good])
   # Selecting the relevant pixels in the input spectrum
         w        = (wavemod >= wlow) & (wavemod <= whi)
         tmp_wave = wavemod[w]
         tmp_flux = fluxmod[w]
         if np.amin(wavemod) > wlow or np.amax(wavemod) < whi:
             continue
   # Interpolate the filter response to data wavelength
         interp   = interp1d(wavefilt[good], fluxfilt[good])
         response = interp(tmp_wave)
         
   # Calculating the magnitude in the desired system
         vega   = 1.0/np.power(tmp_wave, 2)
         f      = np.trapz(tmp_flux * response, x = tmp_wave)
         vega_f = np.trapz(vega     * response, x = tmp_wave)
         mag    = -2.5 * np.log10(f/vega_f)
         fmag   = mag + cfact
         fmag   = fmag + 2.5*np.log10(cvel)-48.6 # oke & gunn 83
         outmag[ii] = fmag
         veg[ii] = vega_f
     return wave_eff, outmag , veg






# Code starts here
# ------------------------------------------------------------------
gal = str(sys.argv[1])






# read the list of SDSS stars
list = np.genfromtxt('list.txt')
ra = list[:, 0]
dec = list[:, 1]

# Read the filters ( H alpha, r LBT ..etc)
wave_Ha, trans_Ha = read_filter('filters/INT_WFC.Halpha.dat')
wave_LBTr, trans_LBTr = read_filter('filters/LBT_LBCR.sdss-r.dat')
wave_LBTg, trans_LBTg = read_filter('filters/LBT_LBCB.sdss-g_1.dat')
filters = ['filters/INT_WFC.Halpha.dat', 'filters/LBT_LBCR.sdss-r.dat', 'filters/LBT_LBCB.sdss-g_1.dat']




#
mag_ha = []        #  HA
mag_lbt_r = []     #  LBT R
mag_lbt_g = []      #  LBT g
g_r = []            # color g-r
ha_r = []           # color Ha-r
ra_ = []
ha_g = []
dec_ = []

veg_corrha = []        #  HA
veg_corr_r = []     #  LBT R
veg_corr_g = []     #  LBT g

#/ np.mean(flux)
# Read all the Sed
for i in np.arange(0, len(ra), 1):
    # Read coordinates
    RA_str = str(ra[i])
    DE_str = str(dec[i])
    outfile = np.genfromtxt('seds/spectra_' + RA_str + '-' + DE_str + '.fits.txt')
    lamdba = 10**(outfile[:, 1])
    flux = outfile[:, 0]
    # here I give the normalized flux
    wave_eff_ref, mag_ref, veg = compute_mags(lamdba, flux , filters)
    print("Mags", mag_ref)
    print("Fluxes", 10**(-0.4*mag_ref))

    mag_ha.append( mag_ref[0] )
    mag_lbt_r.append( mag_ref[1])
    mag_lbt_g.append( mag_ref[2])
    
    veg_corrha.append(veg[0])
    veg_corr_r.append(veg[1])
    veg_corr_g.append(veg[2])
        
    g_r.append(mag_ref[2] - mag_ref[1])
    ha_r.append(mag_ref[0] - mag_ref[1])
    ha_g.append(mag_ref[0] - mag_ref[2])

    ra_.append(RA_str)
    dec_.append(DE_str)


f = np.array(-2.5 * np.log10(mag_ha) / -2.5 * np.log10(mag_lbt_r))

# Now I want to save for each model (each age)
# the magnitudes and the factor

table_to_save = Table([ ra_, dec_, mag_ha, mag_lbt_r, mag_lbt_g, g_r, ha_r ,f, veg_corrha, veg_corr_r, veg_corr_g],
                    names=[ 'ra', 'dec', 'mag_ha', 'mag_lbt_r', 'mag_lbt_g', 'g_r', 'ha_r','fha/fr', 'veg_corrha', 'veg_corr_r', 'veg_corr_g'])
table_to_save.write('build-' + gal + '/table_mag_reference' + gal + '.txt', format='ascii.commented_header', overwrite=True)



cfact    = -48.6+ 2.5*np.log10(2.99792458e18) + 5.0 * np.log10(1.7684E8 * 1e-5)
# This is in magnitudes but with the flux in hertz della stella
mag_nu = np.reshape(-2.5*np.log10(flux*np.power(lamdba , 2)) + cfact, -1)


# Plot fluxes subgroups

# READ the Fluxes and colors of the stars
ha_s = 10**(-0.4 * np.array(mag_ha))
r_s = 10**(-0.4 * np.array(mag_lbt_r))


val = np.genfromtxt('build-' + gal + '/output-catalogs/range.txt')
min = val[0]
max = val[1]




# Plot to understand
colors = cm.tab20c(np.linspace(0, 20, 73))
fig, (ax, ax1) = plt.subplots(figsize=(16, 12), nrows=1, ncols=2, sharex=False)


group = np.zeros((len(ha_s), 3))
        
for i in np.arange( 0, len(ha_s), 1):
  if min < g_r[i] < max:
      if ha_s[i]/r_s[i] < 100:
        ax.scatter(r_s[i], ha_s[i] ,color = 'purple', s=95, marker='X')
        group[i,:] = [r_s[i], ha_s[i], g_r[i]]
    
  else:
    ax.scatter(r_s[i], ha_s[i] ,color = 'orange', s=95, marker='X')
  
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel(' F r-Sloan', fontsize=25)
ax.set_ylabel('F Ha', fontsize=25)





for i in np.arange( 0, len(ha_s), 1):
  if min < g_r[i] < max:
    if ha_s[i]/r_s[i] < 1.2:
       ax1.scatter(g_r[i] , ha_s[i] / r_s[i], color = 'purple', s=95, marker='X', label = 'range in my field')
  else:
    ax1.scatter(g_r[i] , ha_s[i] / r_s[i], color = 'orange', s=95, marker='X',label = 'from sdss spectra')
    
p = np.polyfit(group[:,0] ,group[:,1],  1)
print(p)

#ax.plot(group[:,0], group[:,0]*p[0] + p[1], color = 'red',alpha=0.9, label="1deg Polyfit " + str(p[0]))

# Save values into a file
file = open('build-' + gal + '/pol_coefficients.txt', "w")
file.write(str(p[0]))
file.write(" ")
file.write(str(p[1]))


ax1.set_ylabel(' F_Ha/F_r', fontsize=25)
ax1.set_xlabel('g-r', fontsize=25)
ax1.set_xlim(-0.25,1.5)
ax1.set_ylim(0.9,1.2)

ax.legend(fontsize=25,loc=4)

plt.subplots_adjust(hspace=0.0, wspace=0.2)
plt.savefig('build-' + gal + '/plots/f.png', dpi=300)
plt.cla()


# Now save the tables
table_to_save = Table([ group[:,0], group[:,1], group[:,2], np.array(group[:,1]) / np.array(group[:,0])],
                    names=[ 'flux_lbt_r','flux_ha', 'g_r','fha/fr'])
table_to_save.write('build-' + gal + '/output-catalogs/table_subgroup.txt', format='ascii.commented_header', overwrite=True)
