#! /bin/bash

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
