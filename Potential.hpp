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
 * @file Potential.hpp
 *
 * @brief External potential.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef POTENTIAL_HPP
#define POTENTIAL_HPP

#if POTENTIAL == POTENTIAL_PM_SELF_GRAVITY
#include "FFTWGravitySolver.hpp"
#endif

/**
 * @brief Initialize the gravity solver.
 *
 * @param ncell Number of cells in the grid.
 */
#if POTENTIAL == POTENTIAL_PM_SELF_GRAVITY
#define gravity_init(ncell)                                                    \
  double *rho_grav = new double[ncell];                                        \
  double *a_grav = new double[ncell];                                          \
  FFTWGravitySolver pm_solver(ncell, RMAX - RMIN);
#else
#define gravity_init(ncell)
#endif

/**
 * @brief Clean up gravity solver.
 */
#if POTENTIAL == POTENTIAL_PM_SELF_GRAVITY
#define free_gravity()                                                         \
  delete[] rho_grav;                                                           \
  delete[] a_grav;
#else
#define free_gravity()
#endif

/**
 * @brief Add the gravitational acceleration.
 *
 * @param cells Cells to update.
 * @param ncell Number of cells.
 * @param dt_momentum Momentum update time step.
 * @param dt_energy Energy update time step.
 */
#if POTENTIAL == POTENTIAL_POINT_MASS
#define do_gravity(cells, ncell, dt_momentum,                                  \
                   dt_energy) /* add gravitational acceleration */             \
  _Pragma("omp parallel for") for (uint_fast32_t i = 1; i < ncell + 1; ++i) {  \
    const double r = cells[i]._midpoint;                                       \
    const double a = -G_INTERNAL * MASS_POINT_MASS / (r * r);                  \
    cells[i]._a = a;                                                           \
    cells[i]._pot = -G_INTERNAL * MASS_POINT_MASS / r;                         \
    const double m = cells[i]._V * cells[i]._rho;                              \
    cells[i]._p += 0.5 * dt_momentum * cells[i]._a * m;                        \
  }
#elif POTENTIAL == POTENTIAL_SELF_GRAVITY
#define do_gravity(cells, ncell, dt_momentum, dt_energy)                       \
  {                                                                            \
    double Mtot = 0.;                                                          \
    /* Get the acceleration for each cell (needs to be evaluated in a fixed    \
     * order, serially). */                                                    \
    for (uint_fast32_t i = 1; i < ncell + 1; ++i) {                            \
      const double mcell = Mtot + cells[i]._rho * cells[i]._V_real_half;       \
      const double r = cells[i]._midpoint;                                     \
      cells[i]._a = -G_INTERNAL * mcell / (r * r);                             \
      cells[i]._pot = -G_INTERNAL * mcell / r;                                 \
      Mtot += cells[i]._rho * cells[i]._V_real;                                \
    }                                                                          \
    /* Now apply gravity to each cell. */                                      \
    _Pragma("omp parallel for") for (uint_fast32_t i = 1; i < ncell + 1;       \
                                     ++i) {                                    \
      if (cells[i]._m > 0.) {                                                  \
        const double gravfac = 0.5 * cells[i]._a;                              \
        cells[i]._E += gravfac * cells[i]._p * dt_energy;                      \
        cells[i]._p += gravfac * cells[i]._m * dt_momentum;                    \
      }                                                                        \
      assert_condition(cells[i]._E >= 0., "cells[%" PRIiFAST32 "]._E = %g", i, \
                       cells[i]._E);                                           \
    }                                                                          \
  }
#elif POTENTIAL == POTENTIAL_PM_SELF_GRAVITY
#define do_gravity(cells, ncell, dt_momentum, dt_energy)                       \
  {                                                                            \
    for (uint_fast32_t i = 1; i < ncell + 1; ++i) {                            \
      rho_grav[i - 1] = cells[i]._rho;                                         \
    }                                                                          \
    pm_solver.compute_accelerations(rho_grav, a_grav);                         \
    for (uint_fast32_t i = 1; i < ncell + 1; ++i) {                            \
      cells[i]._a = a_grav[i - 1] * G_INTERNAL;                                \
    }                                                                          \
    /* Now apply gravity to each cell. */                                      \
    _Pragma("omp parallel for") for (uint_fast32_t i = 1; i < ncell + 1;       \
                                     ++i) {                                    \
      if (cells[i]._m > 0.) {                                                  \
        const double gravfac = 0.5 * cells[i]._a;                              \
        cells[i]._E += gravfac * cells[i]._p * dt_momentum;                    \
        cells[i]._E += gravfac * cells[i]._mflux * dt_energy;                  \
        cells[i]._p += gravfac * cells[i]._m * dt_momentum;                    \
      }                                                                        \
      assert_condition(cells[i]._E >= 0., "cells[%" PRIiFAST32 "]._E = %g", i, \
                       cells[i]._E);                                           \
    }                                                                          \
  }
#elif POTENTIAL == POTENTIAL_NONE
#define do_gravity(cells, ncell, dt_momentum, dt_energy) ;
#endif

/**
 * @brief Do the gravitational half time step prediction for the given cell.
 *
 * @param cell Cell.
 * @param half_dt Half the particle time step (in internal units of T).
 */
#if POTENTIAL != POTENTIAL_NONE
#define add_gravitational_prediction(cell, half_dt)                            \
  if (cell._rho > 0.) {                                                        \
    cell._u += 0. * half_dt * cell._a;                                         \
  }
#else
#define add_gravitational_prediction(cell, half_dt)
#endif

#endif // POTENTIAL_HPP
