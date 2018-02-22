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
 * @file SafeParameters.hpp
 *
 * @brief File that should be included by all files that need parameters.
 * Contains sanity checks on the chosen parameter values. Parameter values
 * should be set in Parameters.hpp.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef SAFEPARAMETERS_HPP
#define SAFEPARAMETERS_HPP

// first include the option names
#include "OptionNames.hpp"

// now include the actual parameters
#include "Parameters.hpp"

// check the parameters with options

/*! @brief Macro that prints the actual value of a macro. */
#define value_of_macro_as_string(macro) #macro
/*! @brief Macro that prints a macro name and its actual value. */
#define value_of_macro(macro) #macro ": " value_of_macro_as_string(macro)

// check boundary conditions
#ifndef BOUNDARIES
#error "No boundary conditions selected!"
#else
#if BOUNDARIES != BOUNDARIES_OPEN && BOUNDARIES != BOUNDARIES_REFLECTIVE &&    \
    BOUNDARIES != BOUNDARIES_PERIODIC && BOUNDARIES != BOUNDARIES_CUSTOM
#pragma message(value_of_macro(BOUNDARIES))
#error "Invalid boundary conditions selected!"
#endif
#endif

// check the equation of state
#ifndef EOS
#error "No equation of state selected!"
#else
#if EOS != EOS_IDEAL && EOS != EOS_ISOTHERMAL
#pragma message(value_of_macro(EOS))
#error "Invalid equation of state selected!"
#endif
#endif

// check the potential
#ifndef POTENTIAL
#error "No external potential selected!"
#else
#if POTENTIAL != POTENTIAL_NONE && POTENTIAL != POTENTIAL_POINT_MASS
#pragma message(value_of_macro(POTENTIAL))
#error "Invalid potential selected!"
#endif
#endif

// check Riemann solver type
#ifndef RIEMANNSOLVER_TYPE
#error "No Riemann solver selected!"
#else
#if RIEMANNSOLVER_TYPE != RIEMANNSOLVER_TYPE_EXACT &&                          \
    RIEMANNSOLVER_TYPE != RIEMANNSOLVER_TYPE_HLLC
#pragma message(value_of_macro(RIEMANNSOLVER_TYPE))
#error "Invalid Riemann solver selected!"
#endif
#endif

// check dimensionality
#ifndef DIMENSIONALITY
#error "No dimensionality selected!"
#else
#if DIMENSIONALITY != DIMENSIONALITY_1D &&                                     \
    DIMENSIONALITY != DIMENSIONALITY_2D && DIMENSIONALITY != DIMENSIONALITY_3D
#pragma message(value_of_macro(DIMENSIONALITY))
#error "Invalid dimensionality selected!"
#endif
#endif

// check that we are not setting periodic boundaries in 2D and 3D
#if BOUNDARIES == BOUNDARIES_PERIODIC &&                                       \
    (DIMENSIONALITY == DIMENSIONALITY_2D ||                                    \
     DIMENSIONALITY == DIMENSIONALITY_3D)
#error "Periodic boundary conditions only work with DIMENSIONALITY_1D!"
#endif

// check hydro integration scheme order
#ifndef HYDRO_ORDER
#error "No hydro integration scheme order selected!"
#else
#if HYDRO_ORDER != HYDRO_ORDER_1 && HYDRO_ORDER != HYDRO_ORDER_2
#pragma message(value_of_macro(HYDRO_ORDER))
#error "Invalid hydro integration scheme order selected!"
#endif
#endif

// include derived parameters
#include "DerivedParameters.hpp"

#endif // SAFEPARAMETERS_HPP
