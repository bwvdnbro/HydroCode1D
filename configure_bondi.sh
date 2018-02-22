#! /bin/bash

################################################################################
# This file is part of HydroCode1D
# Copyright (C) 2018 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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

if [ $# -lt 1 ]
then
  echo "Usage: ./configure_bondi.sh <PATH TO SOURCE FOLDER>"
  exit 1
fi

cmake -DCMAKE_BUILD_TYPE=Release -Drmin_in_m=0.1 -Drmax_in_m=1. -Dncell=1000 \
      -Dgamma=1.001 -Dmaxtime_in_s=1. -Dnumber_of_snaps=10 \
      -Deos=EOS_ISOTHERMAL -Disothermal_temperature_in_k=1.e-4 \
      -Dboundaries=BOUNDARIES_CUSTOM -Dpotential=POTENTIAL_POINT_MASS \
      -Dg_internal=1. -Dmass_point_mass_in_kg=1. \
      -Dcourant_factor=0.4 -Driemannsolver_type=RIEMANNSOLVER_TYPE_EXACT \
      -Ddimensionality=DIMENSIONALITY_3D -Dhydro_order=HYDRO_ORDER_2 \
      -Dpredefined_test=bondi $1
