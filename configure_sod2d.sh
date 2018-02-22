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
  echo "Usage: ./configure_sod2d.sh <PATH TO SOURCE FOLDER>"
  exit 1
fi

cmake -DCMAKE_BUILD_TYPE=Release -Drmin_in_m=0. -Drmax_in_m=1. -Dncell=1000 \
      -Dgamma=5./3. -Dmaxtime_in_s=0.1 -Dnumber_of_snaps=1 -Deos=EOS_IDEAL \
      -Dboundaries=BOUNDARIES_OPEN -Dpotential=POTENTIAL_NONE \
      -Dcourant_factor=0.4 -Driemannsolver_type=RIEMANNSOLVER_TYPE_EXACT \
      -Ddimensionality=DIMENSIONALITY_2D -Dhydro_order=HYDRO_ORDER_2 \
      -Dpredefined_test=sod $1
