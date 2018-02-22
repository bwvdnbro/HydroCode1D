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
# @file plot_bondi.py
#
# @brief Script to plot the output of the Bondi test.
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

## Bondi

k_in_si = 1.38064852e-23 # m^2 kg s^-1 K^-1
mH_in_si = 1.674e-27 # kg

# input parameters
# physical
mass_point_mass = 1.
T = 1.e-4
bondi_rho = 1.
# practical
r_min = 0.1
r_max = 1.

# derived parameters
cs2 = T * k_in_si / mH_in_si
bondi_r = 0.5 * mass_point_mass / cs2
cs = np.sqrt(cs2)

def bondi(r):
  global cs, cs2, bondi_r, bondi_rho
  u = bondi_r / r
  omega = -u**4 * np.exp(3. - 4. * u)
  v = np.where(r < bondi_r,
                 -cs * np.sqrt(-lambertw(omega, -1).real),
                 -cs * np.sqrt(-lambertw(omega, 0).real))
  rho = -bondi_rho * bondi_r**2 * cs / r**2 / v
  P = cs2 * rho
  return rho, v, P

ra = np.linspace(r_min, r_max, 1000)
rhoa, va, Pa = bondi(ra)

def plot(f, ax):
  ifile = open(f, 'r')
  timeline = ifile.readline()
  time = float(timeline.split()[2])
  ifile.close()
  data = np.loadtxt(f)
  ax[0].semilogy(data[:,0], data[:,1], "-",
                 label = "$t = {0:.1f}~{{\\rm{{}}s}}$".format(time))
  ax[1].plot(data[:,0], data[:,2], "-")

fig, ax = pl.subplots(2, 1, sharex = True)

files = sorted(glob.glob("snapshot_????.txt"))
ncolor = len(files)
cm = pl.get_cmap("gist_rainbow")
for axcell in ax:
  axcell.set_prop_cycle(
    cycler.cycler("color", [cm(1. * i / ncolor) for i in range(ncolor)]))
  axcell.axvline(x = bondi_r, linestyle = "--", color = "k", linewidth = 0.8)
for f in files:
  plot(f, ax)

ax[0].plot(ra, rhoa, "k--", label = "Bondi")
ax[1].plot(ra, va, "k--")
ax[0].set_ylabel("$\\rho{}$ (kg m$^{-3}$)")
ax[1].set_ylabel("$v$ (m s$^{-1}$)")

ax[0].legend(loc = "best", ncol = 2)

ax[1].set_xlabel("$r$ (m)")

pl.tight_layout()
pl.savefig("bondi.png")
pl.close()
