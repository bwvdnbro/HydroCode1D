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
/*! @brief Periodic boundaries: outgoing flows enter the box through the other
 *  side (note that this only works for 1D). */
#define BOUNDARIES_PERIODIC 3
/*! @brief Spherical open boundaries: reflective left boundary and open right
 *  boundary. */
#define BOUNDARIES_SPHERICAL 4
/*! @brief Custom boundaries: boundary conditions are provided by the user in
 *  (UserInput.hpp). */
#define BOUNDARIES_CUSTOM 5

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
/*! @brief Spherical gravity. */
#define POTENTIAL_SELF_GRAVITY 3
/*! @brief 1D particle mesh self-gravity. */
#define POTENTIAL_PM_SELF_GRAVITY 4

// Possible types of Riemann solver

/*! @brief Exact Riemann solver. */
#define RIEMANNSOLVER_TYPE_EXACT 1
/*! @brief HLLC Riemann solver. */
#define RIEMANNSOLVER_TYPE_HLLC 2

// Possible types of dimensionality

/*! @brief 1D solver. */
#define DIMENSIONALITY_1D 1
/*! @brief 2D spherically symmetric solver. */
#define DIMENSIONALITY_2D 2
/*! @brief 3D spherically symmetric solver. */
#define DIMENSIONALITY_3D 3

// Possible hydro integration scheme orders

/*! @brief First order. */
#define HYDRO_ORDER_1 1
/*! @brief Second order. */
#define HYDRO_ORDER_2 2

// Possible types of time stepping

/*! @brief Normal time stepping. */
#define TIMELINE_NORMAL 1
/*! @brief Co-moving time integration. */
#define TIMELINE_COMOVING 2

#endif // OPTIONNAMES_HPP
