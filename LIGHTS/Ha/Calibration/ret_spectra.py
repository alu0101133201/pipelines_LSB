# Original author:
# Giulia Golini <giulia.golini@gmail.com>
# Contributing author(s)
# Copyright (C) 2020, Giulia Golini.
#
# This Python script is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your
# option) any later version.
#
# This Python script is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details. See <http://www.gnu.org/licenses/>.

# System imports

import os
import pdb
import sys
import warnings


from astroquery.sdss import SDSS
from astropy import coordinates as coords
from astropy.io import fits
import matplotlib.pyplot   as plt
import numpy as np





# read the list

# read the list
list = np.genfromtxt('list.txt')
ra = list[:, 0]
dec = list[:, 1]





for i in np.arange(0, len(ra), 1):
    # Read coordinates
    RA_str = str(ra[i])
    DE_str = str(dec[i])
    outfile = 'spectra/spectra_' + RA_str + '-' + DE_str + '.fits'

    pos = coords.SkyCoord(RA_str + ' ' + DE_str, unit="deg", frame='icrs')
    print(RA_str + ' ' + DE_str)
    
    # Do the actual query
    xid = SDSS.query_region(pos, spectro=True) #, data_release=12)
    sp = SDSS.get_spectra(matches=xid, data_release=12)
    # Save file
    sp[0].writeto(outfile, overwrite=True)
