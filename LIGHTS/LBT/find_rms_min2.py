# Code to find the minimum std value of the sky

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

# 3rd parties


import numpy as np
import matplotlib.pyplot as plt





filter = sys.argv[1]        # band
start = int(sys.argv[2])
end = int(sys.argv[3])
h = str(sys.argv[4])        #n ccd
path_images = str(sys.argv[5])
obj = str(sys.argv[6])

rms = np.full(end+1, np.nan)
# Plot
# Read std values
for ii in range(start, end+1 ,1):
    n = str(ii)
    cat = np.genfromtxt(path_images + '/data-reduction/noise-sky-2-2/' + obj + '_Sloan-' + filter + '_' + n + '_ccd' + h + '.txt')
    rms[ii] = cat[1]

rms_min = np.amin(rms[1:])

file = open(path_images + '/data-reduction/rms_min_val-2_ccd' + h + '.txt', "w")
file.write(str(rms_min))

