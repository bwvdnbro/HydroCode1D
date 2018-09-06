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
  uint_fast64_t current_integer_dt = snaptime;                                 \
  uint_fast64_t isnap = 0;
#elif TIMELINE_TYPE == TIMELINE_COMOVING
#define timeline_init()                                                        \
  const double initial_dt = 0.;                                                \
  uint_fast64_t isnap = 0;                                                     \
  const double maxtime = MAXTIME;                                              \
  const uint_fast64_t integer_maxtime = 0x8000000000000000;                    \
  uint_fast64_t current_integer_time = 0;                                      \
  const double max_physical_dt = initial_dt;                                   \
  uint_fast64_t snaptime = integer_maxtime / NUMBER_OF_SNAPS;                  \
  uint_fast64_t current_integer_dt = snaptime;                                 \
  const double time_conversion_factor = maxtime / integer_maxtime;
#endif

/**
 * @brief Set the timesteps for all cells.
 *
 * @param Maximum allowed physical timestep.
 */
#if TIMELINE_TYPE == TIMELINE_NORMAL
#define set_timesteps(min_physical_dt)                                         \
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
#define set_timesteps(min_physical_dt)
#endif

/**
 * @brief Get the current physical time of the simulation.
 *
 * @return Current physical time of the simulation.
 */
#if TIMELINE_TYPE == TIMELINE_NORMAL
#define current_physical_time()                                                \
  (current_integer_time * maxtime / integer_maxtime)
#elif TIMELINE_TYPE == TIMELINE_COMOVING
#define current_physical_time() 0.
#endif

/**
 * @brief Write a snapshot?
 */
#if TIMELINE_TYPE == TIMELINE_NORMAL
#define do_write_snapshot() (current_integer_time >= isnap * snaptime)
#elif TIMELINE_TYPE == TIMELINE_COMOVING
#define do_write_snapshot() false
#endif

/**
 * @brief Step forward in time.
 */
#if TIMELINE_TYPE == TIMELINE_NORMAL
#define do_timestep() current_integer_time += current_integer_dt;
#elif TIMELINE_TYPE == TIMELINE_COMOVING
#define do_timestep()
#endif

#endif // TIMELINE_HPP
