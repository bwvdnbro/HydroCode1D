/*******************************************************************************
 * This file is part of HydroCode1D
 * Copyright (C) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 *               2018 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 *
 * HydroCodeSpherical1D is free software: you can redistribute it and/or modify
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
 * @file Boundaries.hpp
 *
 * @brief Boundary conditions.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef BOUNDARIES_HPP
#define BOUNDARIES_HPP

#include "SafeParameters.hpp"

/**
 * @brief Apply boundary conditions after the primitive variable conversion.
 *
 * @param cells Cells to update.
 * @param ncell Number of cells.
 */
#if BOUNDARIES == BOUNDARIES_OPEN
#define boundary_conditions_primitive_variables(cells, ncell)                  \
  /* just mirror the values across the boundary */                             \
  cells[0]._dt = cells[1]._dt;                                                 \
  cells[0]._rho = cells[1]._rho;                                               \
  cells[0]._u = cells[1]._u;                                                   \
  cells[0]._P = cells[1]._P;                                                   \
  cells[ncell + 1]._dt = cells[ncell]._dt;                                     \
  cells[ncell + 1]._rho = cells[ncell]._rho;                                   \
  cells[ncell + 1]._u = cells[ncell]._u;                                       \
  cells[ncell + 1]._P = cells[ncell]._P;
#elif BOUNDARIES == BOUNDARIES_REFLECTIVE
#define boundary_conditions_primitive_variables(cells, ncell)                  \
  /* mirror the variables and reverse the sign of the velocity */              \
  cells[0]._dt = cells[1]._dt;                                                 \
  cells[0]._rho = cells[1]._rho;                                               \
  cells[0]._u = -cells[1]._u;                                                  \
  cells[0]._P = cells[1]._P;                                                   \
  cells[ncell + 1]._dt = cells[ncell]._dt;                                     \
  cells[ncell + 1]._rho = cells[ncell]._rho;                                   \
  cells[ncell + 1]._u = -cells[ncell]._u;                                      \
  cells[ncell + 1]._P = cells[ncell]._P;
#elif BOUNDARIES == BOUNDARIES_PERIODIC
#define boundary_conditions_primitive_variables(cells, ncell)                  \
  cells[0]._dt = cells[ncell]._dt;                                             \
  cells[0]._rho = cells[ncell]._rho;                                           \
  cells[0]._u = cells[ncell]._u;                                               \
  cells[0]._P = cells[ncell]._P;                                               \
  cells[ncell + 1]._dt = cells[1]._dt;                                         \
  cells[ncell + 1]._rho = cells[1]._rho;                                       \
  cells[ncell + 1]._u = cells[1]._u;                                           \
  cells[ncell + 1]._P = cells[1]._P;
#elif BOUNDARIES == BOUNDARIES_CUSTOM
#define boundary_conditions_primitive_variables(cells, ncell)                  \
  get_left_boundary(cells[1]._dt, cells[1]._rho, cells[1]._u, cells[1]._P,     \
                    cells[0]._midpoint, cells[0]._dt, cells[0]._rho,           \
                    cells[0]._u, cells[0]._P);                                 \
  get_right_boundary(cells[ncell]._dt, cells[ncell]._rho, cells[ncell]._u,     \
                     cells[ncell]._P, cells[ncell + 1]._midpoint,              \
                     cells[ncell + 1]._dt, cells[ncell + 1]._rho,              \
                     cells[ncell + 1]._u, cells[ncell + 1]._P);
#endif

/**
 * @brief Apply boundary conditions after the gradient computation.
 *
 * @param cells Cells to update.
 * @param ncell Number of cells.
 */
#if BOUNDARIES == BOUNDARIES_OPEN
#define boundary_conditions_gradients(cells, ncell)                            \
  /* just mirror the gradients across the boundary */                          \
  cells[0]._grad_rho = cells[1]._grad_rho;                                     \
  cells[0]._grad_u = cells[1]._grad_u;                                         \
  cells[0]._grad_P = cells[1]._grad_P;                                         \
  cells[ncell + 1]._grad_rho = cells[ncell]._grad_rho;                         \
  cells[ncell + 1]._grad_u = cells[ncell]._grad_u;                             \
  cells[ncell + 1]._grad_P = cells[ncell]._grad_P;
#elif BOUNDARIES == BOUNDARIES_REFLECTIVE
#define boundary_conditions_gradients(cells, ncell)                            \
  /* reverse the sign of the gradients and do some magic to get an accurate    \
     gradient for the velocities */                                            \
  cells[0]._grad_rho = -cells[1]._grad_rho;                                    \
  cells[0]._grad_u =                                                           \
      -cells[1]._grad_u - cells[1]._grad_u -                                   \
      4. * cells[0]._u / (cells[0]._midpoint - cells[1]._midpoint);            \
  cells[0]._grad_P = -cells[1]._grad_P;                                        \
  cells[ncell + 1]._grad_rho = -cells[ncell]._grad_rho;                        \
  cells[ncell + 1]._grad_u =                                                   \
      -cells[ncell]._grad_u -                                                  \
      4. * cells[ncell]._u /                                                   \
          (cells[ncell]._midpoint - cells[ncell + 1]._midpoint);               \
  cells[ncell + 1]._grad_P = -cells[ncell]._grad_P;
#elif BOUNDARIES == BOUNDARIES_PERIODIC
#define boundary_conditions_gradients(cells, ncell)                            \
  cells[0]._grad_rho = cells[ncell]._grad_rho;                                 \
  cells[0]._grad_u = cells[ncell]._grad_u;                                     \
  cells[0]._grad_P = cells[ncell]._grad_P;                                     \
  cells[ncell + 1]._grad_rho = cells[1]._grad_rho;                             \
  cells[ncell + 1]._grad_u = cells[1]._grad_u;                                 \
  cells[ncell + 1]._grad_P = cells[1]._grad_P;
#elif BOUNDARIES == BOUNDARIES_CUSTOM
#define boundary_conditions_gradients(cells, ncell)                            \
  /* set all gradients to 0 for now */                                         \
  cells[0]._grad_rho = 0.;                                                     \
  cells[0]._grad_u = 0.;                                                       \
  cells[0]._grad_P = 0.;                                                       \
  cells[ncell + 1]._grad_rho = 0.;                                             \
  cells[ncell + 1]._grad_u = 0.;                                               \
  cells[ncell + 1]._grad_P = 0.;
#endif

#endif // BOUNDARIES_HPP
