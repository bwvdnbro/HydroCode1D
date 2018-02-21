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
 * @file OptionNames.hpp
 *
 * @brief Useful aliases for define options.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef OPTIONNAMES_HPP
#define OPTIONNAMES_HPP

// Possible types of boundary conditions.

/*! @brief Open boundaries: inflow or outflow depending on the local flow
 *  velocity */
#define BOUNDARIES_OPEN 1
/*! @brief Reflective boundaries: outgoing flows are reflected inwards. */
#define BOUNDARIES_REFLECTIVE 2

// Possible types of equation of state.

/*! @brief Ideal gas (adiabatic) equation of state. */
#define EOS_IDEAL 1
/*! @brief Isothermal equation of state. */
#define EOS_ISOTHERMAL 2

// Possible types of external potentials.

/*! @brief No external potential (no gravity). */
#define POTENTIAL_NONE 1
/*! @brief Point mass external potential. */
#define POTENTIAL_POINT_MASS 2

// Possible types of Riemann solver

/*! @brief Exact Riemann solver. */
#define RIEMANNSOLVER_TYPE_EXACT 1
/*! @brief HLLC Riemann solver. */
#define RIEMANNSOLVER_TYPE_HLLC 2

// Possible types of dimensionality

/*! @brief 1D solver. */
#define DIMENSIONALITY_1D 1
/*! @brief 3D spherically symmetric solver. */
#define DIMENSIONALITY_3D 2

#endif // OPTIONNAMES_HPP
