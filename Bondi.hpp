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
 * @file Bondi.hpp
 *
 * @brief All Bondi problem specific code.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef BONDI_HPP
#define BONDI_HPP

#include "Cell.hpp"           // Cell classe
#include "LambertW.hpp"       // Lambert W function implementation
#include "SafeParameters.hpp" // Safe way to include Parameters.hpp

#include <cmath>

/*! @brief Bondi density: density at the neutral Bondi radius (in internal units
 *  of M L^-3). */
#define BONDI_DENSITY (BONDI_DENSITY_IN_SI / UNIT_DENSITY_IN_SI)

/*! @brief Neutral Bondi radius (in internal units of L). */
#define RBONDI (0.5 * G * MASS_POINT_MASS / ISOTHERMAL_C_SQUARED)

// equation of state functionality for EOS_BONDI
#if EOS == EOS_BONDI

/**
 * @brief Initialize ionisation variables.
 *
 * We need to initialize the smooth transition (if LINEAR_TRANSITION was
 * selected when configuring the code), and the ionising luminosity.
 * We also open the log file in which we will write the ionisation radius as a
 * function of time.
 */
#define ionization_initialize()                                                \
  double rion_old = 0.;                                                        \
                                                                               \
  std::ofstream bondi_rfile("ionization_radius.dat");                          \
                                                                               \
  const double bondi_S =                                                       \
      (transition_width > 0.) ? 3. / (4. * transition_width) : 0.;             \
  const double bondi_A = -32. * bondi_S * bondi_S * bondi_S / 27.;             \
                                                                               \
  /* correction factor to grow the central mass with the accreted material     \
     currently disabled */                                                     \
  /*const double bondi_rmin = RMIN - CELLSIZE;                                 \
  const double bondi_volume_correction_factor = 4. * M_PI / 3. / CELLSIZE *    \
    (RMIN * RMIN * RMIN - bondi_rmin * bondi_rmin * bondi_rmin);*/             \
  const double bondi_volume_correction_factor = 0.;                            \
  std::cout << "Bondi volume correction factor: "                              \
            << bondi_volume_correction_factor << std::endl;                    \
                                                                               \
  /* set the Q value to the value that is needed to ionise out until the       \
     requested ionisation radius */                                            \
  _Pragma("omp parallel for") for (uint_fast32_t i = 1; i < ncell + 1; ++i) {  \
    const double rmin = cells[i]._lowlim;                                      \
    const double rmax = cells[i]._uplim;                                       \
    const double Vshell = (rmax * rmax * rmax - rmin * rmin * rmin) / 3.;      \
    const double Cshell = Vshell * cells[i]._rho * cells[i]._rho;              \
    cells[i]._nfac = Cshell;                                                   \
  }                                                                            \
  double const_bondi_Q = 0.;                                                   \
  for (uint_fast32_t i = 1; i < ncell + 1; ++i) {                              \
    const double rmin = cells[i]._lowlim;                                      \
    const double rmax = cells[i]._uplim;                                       \
    if (rmax < INITIAL_IONIZATION_RADIUS) {                                    \
      const_bondi_Q += cells[i]._nfac;                                         \
    } else if (rmin < INITIAL_IONIZATION_RADIUS) {                             \
      const double ifac =                                                      \
          (INITIAL_IONIZATION_RADIUS * INITIAL_IONIZATION_RADIUS *             \
               INITIAL_IONIZATION_RADIUS -                                     \
           rmin * rmin * rmin) /                                               \
          (rmax * rmax * rmax - rmin * rmin * rmin);                           \
      const double Cshell = cells[i]._nfac * ifac;                             \
      const_bondi_Q += Cshell;                                                 \
    }                                                                          \
  }                                                                            \
  std::cout << "Bondi Q: " << const_bondi_Q << std::endl;                      \
  /* current value of the central mass (only used to increase the luminosity   \
     over time). Currently not really used. */                                 \
  double central_mass = MASS_POINT_MASS;

/**
 * @brief Set the initial value for the pressure of the given cell.
 *
 * We simply use the isothermal equation of state for this, as our scheme does
 * not use the initial pressure anyway (it's only used to initialise the total
 * energy, which is ignored in our integration scheme).
 *
 * @param cell Cell.
 */
#define initial_pressure(cell) cell._P = ISOTHERMAL_C_SQUARED * cell._rho

/**
 * @brief Conversion function called during the primitive variable conversion
 * for the given cell.
 *
 * This is the special equation of state that increase the pressure for ionised
 * gas.
 *
 * @param cell Cell.
 */
#define update_pressure(cell)                                                  \
  const double nfac = cells[i]._nfac;                                          \
  const double ifac = 1. - nfac;                                               \
  cells[i]._P = ISOTHERMAL_C_SQUARED * cells[i]._rho *                         \
                (bondi_pressure_contrast * ifac + nfac);

/**
 * @brief Code to determine the neutral fraction of the cells.
 *
 * This code does three loops over all cells (of which two are done in
 * parallel):
 *  - a first loop to compute the integrated density squared for each cell
 *  - a second loop that figures out when the summed integrated density squared
 *    equals the target luminosity and computes the ionisation radius
 *  - a final loop that sets the neutral fraction once the ionisation radius is
 *    known
 */
#define do_ionization()                                                        \
  /* first loop */                                                             \
  _Pragma("omp parallel for") for (uint_fast32_t i = 1; i < ncell + 1; ++i) {  \
    const double rmin = cells[i]._lowlim;                                      \
    const double rmax = cells[i]._uplim;                                       \
    const double Vshell = (rmax * rmax * rmax - rmin * rmin * rmin) / 3.;      \
    const double Cshell = Vshell * cells[i]._rho * cells[i]._rho;              \
    cells[i]._nfac = Cshell;                                                   \
  }                                                                            \
                                                                               \
  /* second loop */                                                            \
  /* note that the get_bondi_Q_factor part defaults to 1 for now */            \
  double Cion =                                                                \
      const_bondi_Q * get_bondi_Q_factor(central_mass / MASS_POINT_MASS);      \
  double rion = 0.;                                                            \
  for (uint_fast32_t i = 1; i < ncell + 1; ++i) {                              \
    if (Cion > 0.) {                                                           \
      const double rmin = cells[i]._lowlim;                                    \
      const double rmax = cells[i]._uplim;                                     \
      const double ifac = std::min(1., Cion / cells[i]._nfac);                 \
      if (ifac < 1.) {                                                         \
        if (ifac > 0.) {                                                       \
          const double nfac = 1. - ifac;                                       \
          rion = std::cbrt(ifac * rmax * rmax * rmax +                         \
                           nfac * rmin * rmin * rmin);                         \
        }                                                                      \
      } else {                                                                 \
        rion = rmax;                                                           \
      }                                                                        \
      /* subtract this shell's ionization budget from the total */             \
      Cion -= ifac * cells[i]._nfac;                                           \
    }                                                                          \
  }                                                                            \
  /* check if we need to write to the file */                                  \
  if (std::abs(rion - rion_old) > 1.e-2 * std::abs(rion + rion_old)) {         \
    const double curtime =                                                     \
        current_integer_time * time_conversion_factor * UNIT_TIME_IN_SI;       \
    const double ionrad = rion * UNIT_LENGTH_IN_SI;                            \
    Cion = const_bondi_Q * get_bondi_Q_factor(central_mass / MASS_POINT_MASS); \
    bondi_rfile.write(reinterpret_cast<const char *>(&curtime),                \
                      sizeof(double));                                         \
    bondi_rfile.write(reinterpret_cast<const char *>(&ionrad),                 \
                      sizeof(double));                                         \
    bondi_rfile.write(reinterpret_cast<const char *>(&Cion), sizeof(double));  \
    bondi_rfile.flush();                                                       \
    rion_old = rion;                                                           \
  }                                                                            \
  /* fix the ionisation radius. This should be a parameter... */               \
  /*rion = INITIAL_IONIZATION_RADIUS;*/                                        \
                                                                               \
  /* third loop */                                                             \
  const double rion_min = rion - 0.5 * transition_width;                       \
  const double rion_max = rion + 0.5 * transition_width;                       \
  _Pragma("omp parallel for") for (uint_fast32_t i = 1; i < ncell + 1; ++i) {  \
    cells[i]._nfac =                                                           \
        get_neutral_fraction(cells[i]._lowlim, cells[i]._uplim, rion,          \
                             rion_min, rion_max, bondi_S, bondi_A);            \
  }

/**
 * @brief Code to handle the mass flux into the inner mask.
 *
 * We set bondi_volume_correction_factor to zero, so this code does not
 * currently do anything.
 *
 * @param mflux Mass flux into the inner mask.
 */
#define flux_into_inner_mask(mflux)                                            \
  central_mass -= mflux * bondi_volume_correction_factor;

#endif // EOS == EOS_BONDI

// boundary condition functionality
#if BOUNDARIES == BOUNDARIES_BONDI

/**
 * @brief Initialize variables used for the boundary conditions.
 *
 * We need to initialize the outer boundary variables.
 */
#define boundary_conditions_initialize()                                       \
  const double bondi_r_inv_high = RBONDI / cells[ncell + 1]._midpoint;         \
  const double bondi_density_high = bondi_density(bondi_r_inv_high);           \
  const double bondi_velocity_high = bondi_velocity(bondi_r_inv_high);         \
  const double bondi_pressure_high = bondi_pressure(bondi_r_inv_high);         \
                                                                               \
  const double bondi_rmax_inv =                                                \
      RBONDI / (cells[ncell + 1]._midpoint + CELLSIZE);                        \
  const double bondi_density_max = bondi_density(bondi_rmax_inv);              \
  const double bondi_velocity_max = bondi_velocity(bondi_rmax_inv);            \
  const double bondi_pressure_max = bondi_pressure(bondi_rmax_inv);

/**
 * @brief Apply boundary conditions after the primitive variable conversion.
 */
#define boundary_conditions_primitive_variables()                              \
  /* impose the Bondi solution at the boundaries */                            \
  /* lower boundary: outflow */                                                \
  cells[0]._rho = cells[1]._rho;                                               \
  cells[0]._u = cells[1]._u;                                                   \
  cells[0]._P = cells[1]._P;                                                   \
                                                                               \
  /* upper boundary: neutral Bondi solution */                                 \
  cells[ncell + 1]._rho = bondi_density_high;                                  \
  cells[ncell + 1]._u = bondi_velocity_high;                                   \
  cells[ncell + 1]._P = bondi_pressure_high;

/**
 * @brief Apply boundary conditions after the gradient computation.
 */
#define boundary_conditions_gradients()                                        \
  /* lower boundary: outflow */                                                \
  {                                                                            \
    cells[0]._grad_rho = cells[1]._grad_rho;                                   \
    cells[0]._grad_u = cells[1]._grad_u;                                       \
    cells[0]._grad_P = cells[1]._grad_P;                                       \
  }                                                                            \
                                                                               \
  /* upper boundary: need to change this to correct expression */              \
  {                                                                            \
    double rhomin, rhoplu, umin, uplu, Pmin, Pplu, dmin, dplu;                 \
    rhomin = cells[ncell + 1]._rho - cells[ncell]._rho;                        \
    rhoplu = cells[ncell + 1]._rho - bondi_density_max;                        \
    umin = cells[ncell + 1]._u - cells[ncell]._u;                              \
    uplu = cells[ncell + 1]._u - bondi_velocity_max;                           \
    Pmin = cells[ncell + 1]._P - cells[ncell]._P;                              \
    Pplu = cells[ncell + 1]._P - bondi_pressure_max;                           \
    if (std::abs(rhomin) < std::abs(rhoplu)) {                                 \
      cells[ncell + 1]._grad_rho = rhomin / CELLSIZE;                          \
    } else {                                                                   \
      cells[ncell + 1]._grad_rho = rhoplu / CELLSIZE;                          \
    }                                                                          \
    if (std::abs(umin) < std::abs(uplu)) {                                     \
      cells[ncell + 1]._grad_u = umin / CELLSIZE;                              \
    } else {                                                                   \
      cells[ncell + 1]._grad_u = uplu / CELLSIZE;                              \
    }                                                                          \
    if (std::abs(Pmin) < std::abs(Pplu)) {                                     \
      cells[ncell + 1]._grad_P = Pmin / CELLSIZE;                              \
    } else {                                                                   \
      cells[ncell + 1]._grad_P = Pplu / CELLSIZE;                              \
    }                                                                          \
  }

#endif // BOUNDARIES == BOUNDARIES_BONDI

/**
 * @brief Squared neutral Bondi velocity divided by the sound speed squared.
 *
 * @param rinv Inverse radius (in units of RBONDI^-1).
 * @return Bondi velocity squared divided by the sound speed squared.
 */
double u2_over_cs2(double rinv) {
  const double lambertarg = -std::exp(3. + 4. * (std::log(rinv) - rinv));
  if (rinv < 1.) {
    return -LambertW::lambert_w(lambertarg, 0);
  } else {
    return -LambertW::lambert_w(lambertarg, -1);
  }
}

/**
 * @brief Get the value of the neutral Bondi density at the given inverse
 * radius.
 *
 * @param rinv Inverse radius (in units of RBONDI^-1).
 * @return Value of the density (in internal units of M L^-3).
 */
double bondi_density(double rinv) {
  // we need to manually disable the density very close to r = 0 to prevent
  // errors in the Lambert W function
  if (rinv < 150.) {
    return BONDI_DENSITY * std::exp(-0.5 * u2_over_cs2(rinv) + 2. * rinv - 1.5);
  } else {
    return 0.;
  }
}

/**
 * @brief Get the value of the neutral Bondi velocity at the given inverse
 * radius.
 *
 * @param rinv Inverse radius (in internal units of L^-1).
 * @return Value of the fluid velocity (in internal units of L T^-1).
 */
double bondi_velocity(double rinv) {
  return -std::sqrt(ISOTHERMAL_C_SQUARED * u2_over_cs2(rinv));
}

/**
 * @brief Get the value of the neutral Bondi pressure at the given inverse
 * radius.
 *
 * @param rinv Inverse radius (in internal units of L^-1).
 * @return Value of the pressure (in internal units of M L^-1 T^-2).
 */
double bondi_pressure(double rinv) {
  return ISOTHERMAL_C_SQUARED * bondi_density(rinv);
}

/**
 * @brief Get the increase in luminosity corresponding to the given increase in
 * central mass.
 *
 * This expression was fitted to the values in Keto (2003).
 *
 * @param M Central mass (in units of the initial central mass).
 * @return Luminosity increase for the given mass increase.
 */
inline static double get_bondi_Q_factor(const double M) {
  return 7.96185873 * (std::pow(M, 2.47692987) - 1.) + 1.;
}

/**
 * @brief Get the neutral fraction within the given interval using an exact
 * integration of the smooth transition expression.
 *
 * @param A Smooth transition parameter A (in internal units of L^-3).
 * @param S Smooth transition parameter S (steepest allowed slope in the
 * transition; in internal units of L^-1).
 * @param rion Ionisation radius (in internal units of L).
 * @param rmin Lower limit of the integration interval (in internal units of L).
 * @param rmax Upper limit of the integration interval (in internal units of L).
 * @return Integrated neutral fraction within the interval (in internal units of
 * L).
 */
inline static double
get_neutral_fraction_integral(const double A, const double S, const double rion,
                              const double rmin, const double rmax) {
  const double rmaxrel = rmax - rion;
  const double rminrel = rmin - rion;
  const double rdiff = rmax - rmin;
  const double rmaxrel2 = rmaxrel * rmaxrel;
  const double rmaxrel4 = rmaxrel2 * rmaxrel2;
  const double rminrel2 = rminrel * rminrel;
  const double rminrel4 = rminrel2 * rminrel2;
  return A * (rmaxrel4 - rminrel4) + S * (rmaxrel2 - rminrel2) + 0.5 * rdiff;
}

/**
 * @brief Get the neutral fraction for the cell with the given midpoint radius.
 *
 * @param rmin Radius of the lower wall of the cell (in internal units of L).
 * @param rmax Radius of the upper wall of the cell (in internal units of L).
 * @param rion Ionisation radius (in internal units of L).
 * @param rion_min Ionisation radius minus smooth transition width (in internal
 * units of L).
 * @param rion_max Ionisation radius plus smooth transition width (in internal
 * units of L).
 * @param A Smooth transition parameter A (in internal units of L^-3).
 * @param S Smooth transition parameter S (steepest allowed slope in the
 * transition; in internal units of L^-1).
 * @return Neutral fraction within the cell.
 */
inline static double get_neutral_fraction(const double rmin, const double rmax,
                                          const double rion,
                                          const double rion_min,
                                          const double rion_max, const double S,
                                          const double A) {
  // HERE WE NEED TO IMPLEMENT LINEAR_TRANSITION
  // Note: we assume rmax - rmin << rion_max - rion_min
  if (rmax < rion_min) {
    return 0.;
  } else if (rmin < rion_min && rmax >= rion_min) {
    return get_neutral_fraction_integral(A, S, rion, rion_min, rmax) /
           (rmax - rmin);
  } else if (rmin >= rion_min && rmax <= rion_max) {
    return get_neutral_fraction_integral(A, S, rion, rmin, rmax) /
           (rmax - rmin);
  } else if (rmin < rion_max && rmax > rion_max) {
    return (get_neutral_fraction_integral(A, S, rion, rmin, rion_max) +
            (rmax - rion_max)) /
           (rmax - rmin);
  } else {
    return 1.;
  }
  //    if(rmax < rion){
  //      return 0;
  //    } else if(rmin < rion){
  //      return (rmax - rion) / (rmax - rmin);
  //    } else {
  //      return 1.;
  //    }
}

// Bondi initial condition functionality
#if IC == IC_BONDI
/**
 * @brief Initialize the given cells.
 *
 * We start with a constant density and velocity everywhere, the value of which
 * corresponds to the inflow values.
 *
 * @param cells Cells to initialize.
 * @param ncell Number of cells.
 */
#define initialize(cells, ncell)                                               \
  _Pragma("omp parallel for") for (unsigned int i = 1; i < ncell + 1; ++i) {   \
    const double r_inv = RBONDI / RMAX;                                        \
    cells[i]._rho = bondi_density(r_inv);                                      \
    cells[i]._u = bondi_velocity(r_inv);                                       \
    cells[i]._P = ISOTHERMAL_C_SQUARED * BONDI_DENSITY;                        \
    const double r2 = cells[i]._midpoint * cells[i]._midpoint;                 \
    cells[i]._a = -G * MASS_POINT_MASS / r2;                                   \
    cells[i]._nfac = 0.;                                                       \
  }

#endif // IC == IC_BONDI

#endif // BONDI_HPP
