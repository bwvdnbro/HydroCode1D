/*******************************************************************************
 * This file is part of HydroCode1D
 * Copyright (C) 2018 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file Timeline.hpp
 *
 * @brief Time integration routines.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef TIMELINE_HPP
#define TIMELINE_HPP

#include "SafeParameters.hpp"

#if TIMELINE_TYPE == TIMELINE_COMOVING
#include "Cosmology.hpp"
#endif

/**
 * @brief Initialize the time stepping variables.
 */
#if TIMELINE_TYPE == TIMELINE_NORMAL
#define timeline_init()                                                        \
  const double maxtime = MAXTIME;                                              \
  const uint_fast64_t integer_maxtime = 0x8000000000000000;                    \
  const double time_conversion_factor = maxtime / integer_maxtime;             \
  const double initial_dt = (MAXTIME / NUMBER_OF_SNAPS);                       \
  const double max_physical_dt = initial_dt;                                   \
  uint_fast64_t current_integer_time = 0;                                      \
  uint_fast64_t snaptime = integer_maxtime / NUMBER_OF_SNAPS;                  \
  uint_fast64_t current_integer_dt = snaptime;
#elif TIMELINE_TYPE == TIMELINE_COMOVING
#define timeline_init()                                                        \
  Cosmology cosmo(1., 0., 0., 0., 0., 0., 1., 0.00990099, 1.);                 \
  double current_scale_factor = 0.00990099;                                    \
  const uint_fast64_t integer_maxtime = 0x8000000000000000;                    \
  const double maxtime = -std::log(0.00990099);                                \
  uint_fast64_t current_integer_time = 0;                                      \
  uint_fast64_t current_integer_dt = 0;                                        \
  double source_dt, gravity_dt_momentum, gravity_dt_energy, hydro_dt;          \
  uint_fast64_t snaptime = integer_maxtime / NUMBER_OF_SNAPS;                  \
  const double time_conversion_factor_inv =                                    \
      static_cast<double>(integer_maxtime) / maxtime;                          \
  const double time_conversion_factor =                                        \
      maxtime / static_cast<double>(integer_maxtime);                          \
  (void)source_dt, (void)gravity_dt_momentum, (void)gravity_dt_energy;
#endif

/**
 * @brief Set the system time step.
 *
 * @param Maximum allowed physical timestep (in s).
 */
#if TIMELINE_TYPE == TIMELINE_NORMAL
#define timeline_set_timestep(min_physical_dt)                                 \
  {                                                                            \
    const uint_fast64_t min_integer_dt =                                       \
        (min_physical_dt / maxtime) * integer_maxtime;                         \
    current_integer_dt = round_power2_down(min_integer_dt);                    \
    while ((integer_maxtime - current_integer_time) % current_integer_dt >     \
           0) {                                                                \
      current_integer_dt >>= 1;                                                \
    }                                                                          \
  }
#elif TIMELINE_TYPE == TIMELINE_COMOVING
#define timeline_set_timestep(min_physical_dt)                                 \
  {                                                                            \
    const uint_fast64_t min_integer_dt = static_cast<uint_fast64_t>(           \
        cosmo.get_scale_factor_interval(current_scale_factor,                  \
                                        min_physical_dt) *                     \
        time_conversion_factor_inv);                                           \
    current_integer_dt = round_power2_down(min_integer_dt);                    \
    while ((integer_maxtime - current_integer_time) % current_integer_dt >     \
           0) {                                                                \
      current_integer_dt >>= 1;                                                \
    }                                                                          \
    const double a_stop =                                                      \
        std::exp(std::log(current_scale_factor) +                              \
                 current_integer_dt * time_conversion_factor);                 \
    source_dt = cosmo.get_factor_a2inv(current_scale_factor, a_stop);          \
    gravity_dt_momentum = cosmo.get_factor_ainv(current_scale_factor, a_stop); \
    gravity_dt_energy = cosmo.get_factor_a(current_scale_factor, a_stop);      \
    hydro_dt = cosmo.get_factor_a2inv(current_scale_factor, a_stop);           \
  }
#endif

/**
 * @brief Do we want a next integration time step?
 *
 * @return True if we want to compute another system time step.
 */
#if TIMELINE_TYPE == TIMELINE_NORMAL
#define timeline_next_step() (current_integer_time < integer_maxtime)
#elif TIMELINE_TYPE == TIMELINE_COMOVING
#define timeline_next_step() (current_integer_time < integer_maxtime)
#endif

/**
 * @brief Get the physical system time step for source terms.
 *
 * @return Physical system time step for source terms (in s).
 */
#if TIMELINE_TYPE == TIMELINE_NORMAL
#define timeline_get_system_source_dt()                                        \
  (current_integer_dt * time_conversion_factor)
#elif TIMELINE_TYPE == TIMELINE_COMOVING
#define timeline_get_system_source_dt() source_dt
#endif

/**
 * @brief Get the physical system time step for gravity terms involving the
 * momentum.
 *
 * @return Physical system time step for gravity terms (in s).
 */
#if TIMELINE_TYPE == TIMELINE_NORMAL
#define timeline_get_system_gravity_dt_momentum()                              \
  (current_integer_dt * time_conversion_factor)
#elif TIMELINE_TYPE == TIMELINE_COMOVING
#define timeline_get_system_gravity_dt_momentum() gravity_dt_momentum
#endif

/**
 * @brief Get the physical system time step for gravity terms involving the
 * energy.
 *
 * @return Physical system time step for gravity terms (in s).
 */
#if TIMELINE_TYPE == TIMELINE_NORMAL
#define timeline_get_system_gravity_dt_energy()                                \
  (current_integer_dt * time_conversion_factor)
#elif TIMELINE_TYPE == TIMELINE_COMOVING
#define timeline_get_system_gravity_dt_energy() gravity_dt_energy
#endif

/**
 * @brief Get the physical system time step for hydro terms.
 *
 * @return Physical system time step for hydro terms (in s).
 */
#if TIMELINE_TYPE == TIMELINE_NORMAL
#define timeline_get_system_hydro_dt()                                         \
  (current_integer_dt * time_conversion_factor)
#elif TIMELINE_TYPE == TIMELINE_COMOVING
#define timeline_get_system_hydro_dt() hydro_dt
#endif

/**
 * @brief Get the maximum allowed physical time step.
 *
 * @return Maximum allowed physical time step (in s).
 */
#if TIMELINE_TYPE == TIMELINE_NORMAL
#define timeline_get_max_physical_dt() (max_physical_dt)
#elif TIMELINE_TYPE == TIMELINE_COMOVING
#define timeline_get_max_physical_dt() 1.e9
#endif

/**
 * @brief Get the current physical time of the simulation.
 *
 * @return Current physical time of the simulation.
 */
#if TIMELINE_TYPE == TIMELINE_NORMAL
#define timeline_get_current_physical_time()                                   \
  (current_integer_time * maxtime / integer_maxtime)
#elif TIMELINE_TYPE == TIMELINE_COMOVING
#define timeline_get_current_physical_time()                                   \
  cosmo.get_time_since_big_bang(current_scale_factor)
#endif

/**
 * @brief Write a snapshot?
 *
 * @return True if we want to write a snapshot dump at the current time.
 */
#if TIMELINE_TYPE == TIMELINE_NORMAL
#define timeline_do_write_snapshot() (current_integer_time >= isnap * snaptime)
#elif TIMELINE_TYPE == TIMELINE_COMOVING
#define timeline_do_write_snapshot() (current_integer_time >= isnap * snaptime)
#endif

/**
 * @brief Print statistics about the system time stepping.
 */
#if TIMELINE_TYPE == TIMELINE_NORMAL
#define timeline_print_system_stats()                                          \
  std::cout << "time " << t << " of " << maxtime << " (" << pct << " %)"       \
            << std::endl;                                                      \
  std::cout << "\t\t\tSystem time step: "                                      \
            << current_integer_dt * time_conversion_factor << std::endl;
#elif TIMELINE_TYPE == TIMELINE_COMOVING
#define timeline_print_system_stats()                                          \
  std::cout << "redshift " << (1. / current_scale_factor - 1.) << std::endl;   \
  std::cout << "\t\t\tSystem time step: "                                      \
            << current_integer_dt * time_conversion_factor << std::endl;       \
  std::cout << "\t\t\tHydro step: " << hydro_dt << std::endl;                  \
  std::cout << "\t\t\tGravity step (momentum): " << gravity_dt_momentum        \
            << std::endl;                                                      \
  std::cout << "\t\t\tGravity step (energy): " << gravity_dt_energy            \
            << std::endl;
#endif

/**
 * @brief Step forward in time.
 */
#if TIMELINE_TYPE == TIMELINE_NORMAL
#define timeline_do_timestep() current_integer_time += current_integer_dt;
#elif TIMELINE_TYPE == TIMELINE_COMOVING
#define timeline_do_timestep()                                                 \
  current_integer_time += current_integer_dt;                                  \
  current_scale_factor = std::exp(                                             \
      ((1. * current_integer_time) / integer_maxtime - 1.) * maxtime);
#endif

/**
 * @brief Get the fraction of the total simulation time that has already been
 * simulated.
 *
 * @return Already simulated fraction (in %).
 */
#if TIMELINE_TYPE == TIMELINE_NORMAL
#define timeline_get_simulated_fraction()                                      \
  (current_integer_time * 100. / integer_maxtime)
#elif TIMELINE_TYPE == TIMELINE_COMOVING
#define timeline_get_simulated_fraction() 100.
#endif

/**
 * @brief Get the index of the first snapshot that will be written.
 *
 * @return Index of the first snapshot file.
 */
#if TIMELINE_TYPE == TIMELINE_NORMAL
#define timeline_get_initial_snap_index() 0
#elif TIMELINE_TYPE == TIMELINE_COMOVING
#define timeline_get_initial_snap_index() 0
#endif

/**
 * @brief Convert all cell variables to internal variables.
 *
 * This method only changes the cell variables if comoving integration is used.
 */
#if TIMELINE_TYPE == TIMELINE_NORMAL
#define timeline_do_variable_conversion(cells, ncell)
#elif TIMELINE_TYPE == TIMELINE_COMOVING
#define timeline_do_variable_conversion(cells, ncell)                          \
  {                                                                            \
    const double astart = 0.00990099;                                          \
    _Pragma("omp parallel for") for (unsigned int i = 1; i < ncell + 1; ++i) { \
      cells[i]._u *= astart;                                                   \
    }                                                                          \
  }
#endif

/**
 * @brief Get the conversion factor from internal velocity to peculiar velocity.
 *
 * @return Conversion factor.
 */
#if TIMELINE_TYPE == TIMELINE_NORMAL
#define timeline_get_velocity_factor() 1.
#elif TIMELINE_TYPE == TIMELINE_COMOVING
#define timeline_get_velocity_factor() (1. / current_scale_factor)
#endif

/**
 * @brief Get the timestamp to add to the snapshot files.
 *
 * @param timestamp Variable to set.
 */
#if TIMELINE_TYPE == TIMELINE_NORMAL
#define timeline_get_timestamp(timestamp)                                      \
  {                                                                            \
    std::stringstream stampstr;                                                \
    stampstr << "# time: "                                                     \
             << (current_integer_time * maxtime / integer_maxtime) << " s";    \
    timestamp = stampstr.str();                                                \
  }
#elif TIMELINE_TYPE == TIMELINE_COMOVING
#define timeline_get_timestamp(timestamp)                                      \
  {                                                                            \
    std::stringstream stampstr;                                                \
    stampstr << "# redshift: " << (1. / current_scale_factor - 1.);            \
    timestamp = stampstr.str();                                                \
  }
#endif

#endif // TIMELINE_HPP
