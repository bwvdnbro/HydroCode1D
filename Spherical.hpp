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
 * @file Spherical.hpp
 *
 * @brief Spherical source terms that make the 1D method 3D spherical.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef SPHERICAL_HPP
#define SPHERICAL_HPP

/*! @brief Norm of the correction source term in 2D and 3D radial symmetry
 *  (this is the value \f$\alpha{}\f$ in Toro (2009) divided by 2). */
#if DIMENSIONALITY == DIMENSIONALITY_2D
#define DIMENSIONALITY_ALPHA 0.5
#elif DIMENSIONALITY == DIMENSIONALITY_3D
#define DIMENSIONALITY_ALPHA 1.
#endif

/**
 * @brief Add the spherical source terms.
 *
 * See Toro, 2009, chapter 17.
 * We use a second order Runge-Kutta step and apply an operator splitting method
 * to couple the source term to the hydro step.
 *
 * @param cells Cells to update.
 * @param ncell Number of cells.
 */
#if DIMENSIONALITY == DIMENSIONALITY_1D
#define add_spherical_source_term(cells, ncell, dt)
#else
#define add_spherical_source_term(cells, ncell, dt)                            \
  _Pragma("omp parallel for") for (uint_fast32_t i = 1; i < ncell + 1; ++i) {  \
    if (cells[i]._m > 0.) {                                                    \
      const double r = cells[i]._midpoint;                                     \
      const double rinv = DIMENSIONALITY_ALPHA / r;                            \
      const double Vinv = 1. / cells[i]._V;                                    \
      const double Ui[3] = {cells[i]._m * Vinv, cells[i]._p * Vinv,            \
                            cells[i]._E * Vinv};                               \
      const double Ui0inv = 1. / Ui[0];                                        \
      const double Ui12 = Ui[1] * Ui[1];                                       \
      const double p1 = (GAMMA - 1.) * (Ui[2] - 0.5 * Ui12 * Ui0inv);          \
      const double K1[3] = {-dt * Ui[1] * rinv, -dt * Ui12 * Ui0inv * rinv,    \
                            -dt * Ui[1] * (Ui[2] + p1) * Ui0inv * rinv};       \
      double U[3] = {Ui[0] + K1[0], Ui[1] + K1[1], Ui[2] + K1[2]};             \
      const double U0inv = 1. / U[0];                                          \
      const double U12 = U[1] * U[1];                                          \
      const double p2 = (GAMMA - 1.) * (U[2] - 0.5 * U12 * U0inv);             \
      const double K2[3] = {-dt * U[1] * rinv, -dt * U12 * U0inv * rinv,       \
                            -dt * U[1] * (U[2] + p2) * U0inv * rinv};          \
      U[0] = Ui[0] + 0.5 * (K1[0] + K2[0]);                                    \
      U[1] = Ui[1] + 0.5 * (K1[1] + K2[1]);                                    \
      U[2] = Ui[2] + 0.5 * (K1[2] + K2[2]);                                    \
                                                                               \
      cells[i]._m = U[0] * cells[i]._V;                                        \
      cells[i]._p = U[1] * cells[i]._V;                                        \
      cells[i]._E = U[2] * cells[i]._V;                                        \
                                                                               \
      assert_condition(cells[i]._m >= 0., "cells[%" PRIiFAST32 "]._m = %g", i, \
                       cells[i]._m);                                           \
      assert_condition(cells[i]._E >= 0., "cells[%" PRIiFAST32 "]._E = %g", i, \
                       cells[i]._E);                                           \
    }                                                                          \
  }
#endif

#endif // SPHERICAL_HPP
