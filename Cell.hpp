/*******************************************************************************
 * This file is part of HydroCode1D
 * Copyright (C) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 *               2018 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 *
 * HydroCode1D is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HydroCode1D is distributed in the hope that it will be useful,
 * but WITOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with HydroCode1D. If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/

/**
 * @file Cell.hpp
 *
 * @brief Single cell of the grid.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef CELL_HPP
#define CELL_HPP

#include <cstdint>

/**
 * @brief Single cell of the grid.
 */
class Cell {
public:
  // primitive (volume dependent) variables

  /*! @brief Density (in internal units of M L^-3). */
  double _rho;
  /*! @brief Fluid velocity (in internal units of L T^-1). */
  double _u;
  /*! @brief Pressure (in internal units of M L^-1 T^-2). */
  double _P;

  /*! @brief Density gradient (in internal units of M L^-4). */
  double _grad_rho;
  /*! @brief Velocity gradient (in internal units of T^-1). */
  double _grad_u;
  /*! @brief Pressure gradient (in internal units of M L^-2 T^-2). */
  double _grad_P;

  // conserved variables

  /*! @brief Mass (in internal units of M). */
  double _m;
  /*! @brief Momentum (in internal units of L M T^-1). */
  double _p;
  /*! @brief Total energy (in internal units of M L^2 T^-2). */
  double _E;
  /*! @brief Mass flux (for gravity; in internal units of M L^2 T). */
  double _mflux;

  // geometrical quantities

  /*! @brief 1D cell volume: radial length of the spherical shell (in internal
   *  units of L). */
  double _V;
  /*! @brief Radial coordinate of the center of the shell (in internal units of
   *  L). */
  double _midpoint;
  /*! @brief Radial coordinate of the lower boundary of the shell (in internal
   *  units of L). */
  double _lowlim;
  /*! @brief Radial coordinate of the upper boundary of the shell (in internal
   *  units of L). */
  double _uplim;
  /*! @brief Real volume of the cell: differs from _V for 2D and 3D spherically
   *  symmetric simulations (in internal units of L^d). */
  double _V_real;
  /*! @brief Volume of the part of the cell enclosed in _lowlim and _midpoint
   *  (in internal units of L^d). */
  double _V_real_half;

  // gravitational quantities

  /*! @brief Gravitational acceleration (in internal units of L T^-2). */
  double _a;

  /*! @brief Gravitational potential (in internal units of L^2 T^-2). */
  double _pot;
};

#endif // CELL_HPP
