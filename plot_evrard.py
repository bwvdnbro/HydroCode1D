#! /usr/bin/python

################################################################################
# This file is part of HydroCode1D
# Copyright (C) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
#               2018 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
#
# HydroCode1D is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# HydroCode1D is distributed in the hope that it will be useful,
# but WITOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with HydroCode1D. If not, see <http://www.gnu.org/licenses/>.
################################################################################

##
# @file plot_evrard.py
#
# @brief Script to plot the output of the Evrard test.
#
# @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
##

import numpy as np
import matplotlib
matplotlib.use("Agg")
import pylab as pl
import scipy.special.lambertw as lambertw
import glob
import cycler

pl.rcParams["text.usetex"] = True
pl.rcParams["figure.figsize"] = (6, 8)
pl.rcParams["font.size"] = 14
pl.rcParams["axes.labelsize"] = 18

def plot(f, ax):
  ifile = open(f, 'r')
  timeline = ifile.readline()
  time = float(timeline.split()[2])
  ifile.close()
  data = np.loadtxt(f)
  ax[0].loglog(data[:,0], data[:,1], "-")
  ax[1].semilogx(data[:,0], data[:,2], "-",
                 label = "$t = {0:.1f}~{{\\rm{{}}s}}$".format(time))

fig, ax = pl.subplots(2, 1, sharex = True)

files = sorted(glob.glob("snapshot_????.txt"))
ncolor = len(files)
cm = pl.get_cmap("gist_rainbow")
for axcell in ax:
  axcell.set_prop_cycle(
    cycler.cycler("color", [cm(1. * i / ncolor) for i in range(ncolor)]))
for f in files:
  plot(f, ax)

ax[0].set_ylabel("$\\rho{}$ (kg m$^{-3}$)")
ax[1].set_ylabel("$v$ (m s$^{-1}$)")

ax[1].legend(loc = "best", ncol = 2)

ax[1].set_xlabel("$r$ (m)")

pl.tight_layout()
pl.savefig("evrard.png")
pl.close()
