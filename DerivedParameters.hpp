/*******************************************************************************
 * This file is part of HydroCodeSpherical1D
 * Copyright (C) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 *
 * HydroCodeSpherical1D is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HydroCodeSpherical1D is distributed in the hope that it will be useful,
 * but WITOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with HydroCodeSpherical1D. If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/

/**
 * @file DerivedParameters.hpp
 *
 * @brief Parameters that are derived from other parameters.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef DERIVEDPARAMETERS_HPP
#define DERIVEDPARAMETERS_HPP

#include "PhysicalConstants.hpp"

/*! @brief Inner boundary radius (in internal units of L). */
#define RMIN (RMIN_IN_M)

/*! @brief Outer boundary radius (in internal units of L). */
#define RMAX (RMAX_IN_M)

/*! @brief Maximum simulation time (in internal units of T). */
#define MAXTIME (MAXTIME_IN_S)

/*! @brief Mass of the central point mass (in internal units of M). */
#define MASS_POINT_MASS (MASS_POINT_MASS_IN_KG)

/*! @brief Size of the simulation "box" (in internal units of L). */
#define BOXSIZE (RMAX - RMIN)

/*! @brief Size of a single "cell" of the simulation (in internal units of
 *  L). */
#define CELLSIZE (BOXSIZE / ncell)

/*! @brief Half the size of a single "cell" of the simulation (in internal units
 *  of L). */
#define HALF_CELLSIZE (0.5 * CELLSIZE)

/*! @brief Isothermal sound speed squared in the neutral region (if
 *  EOS_ISOTHERMAL is selected, in internal units of L T^-1). */
#define ISOTHERMAL_C_SQUARED                                                   \
  (ISOTHERMAL_TEMPERATURE_IN_K * BOLTZMANN_K_IN_SI / HYDROGEN_MASS_IN_SI)

#endif // DERIVEDPARAMETERS_HPP
