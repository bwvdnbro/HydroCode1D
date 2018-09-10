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
 * @file HydroCode1D.cpp
 *
 * @brief Main program.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */

// project includes
#include "Boundaries.hpp"        // for boundary conditions
#include "Cell.hpp"              // Cell class
#include "EOS.hpp"               // for equations of state
#include "HLLCRiemannSolver.hpp" // fast HLLC Riemann solver
#include "IC.hpp"                // general initial condition interface
#include "Potential.hpp"         // (external) gravity
#include "RiemannSolver.hpp"     // slow exact Riemann solver
#include "SafeParameters.hpp"    // safe way to include Parameter.hpp
#include "Spherical.hpp"         // spherical source terms
#include "Timeline.hpp"          // time stepping routines
#include "Timer.hpp"             // program timers

// standard libraries
#include <cfloat>
#include <cinttypes>
#include <cmath>
#include <cstdint>
#include <ctime>
#include <fstream>
#include <iostream>
#include <omp.h>
#include <sstream>

/**
 * @brief Custom assertion macro.
 *
 * @param c Condition to assert.
 * @param m Message to append to the error output.
 * @param ... Extra parameters for message.
 */
#define assert_condition(c, m, ...)                                            \
  if (!(c)) {                                                                  \
    fprintf(stderr, "\n%s:%s():%i: Assertion failed: %s (" m ")!\n\n",         \
            __FILE__, __FUNCTION__, __LINE__, #c, ##__VA_ARGS__);              \
    exit(1);                                                                   \
  }

/**
 * @brief Get the current time as a string.
 *
 * @return Current system time as a string with format YYYY:MM:DD:HH:MM:SS.
 */
static std::string get_timestamp() {
  const std::time_t timestamp = std::time(nullptr);
  const std::tm *time = std::localtime(&timestamp);
  std::stringstream timestream;
  timestream << (time->tm_year + 1900) << ":";
  if (time->tm_mon < 9) {
    timestream << "0";
  }
  timestream << (time->tm_mon + 1) << ":";
  if (time->tm_mday < 10) {
    timestream << "0";
  }
  timestream << time->tm_mday << ":";
  if (time->tm_hour < 10) {
    timestream << "0";
  }
  timestream << time->tm_hour << ":";
  if (time->tm_min < 10) {
    timestream << "0";
  }
  timestream << time->tm_min << ":";
  if (time->tm_sec < 10) {
    timestream << "0";
  }
  timestream << time->tm_sec;
  return timestream.str();
}

/**
 * @brief Add code information to an open output file.
 *
 * @param ofile Output file to add to.
 */
void add_code_block_to_file(std::ofstream &ofile) {
  ofile << "# Output by HydroCode1D on " << get_timestamp() << "\n";
  ofile << "# https://github.com/bwvdnbro/HydroCode1D\n";
  ofile << "# Git version: " << git_info << "\n";
  ofile << "# Configuration:\n";
  for (uint_fast32_t i = 0; i < CONFIGURATION_NUMBER; ++i) {
    ofile << "#  " << configuration_keys[i] << ": " << configuration_values[i]
          << "\n";
  }
  ofile << "#\n";
}

/**
 * @brief Write a snapshot with the given index.
 *
 * @param istep Index of the snapshot file.
 * @param cells Cells to write.
 * @param ncell Number of cells.
 */
void write_snapshot(uint_fast64_t istep, double time, const Cell *cells,
                    const unsigned int ncell) {
  std::stringstream filename;
  filename << "snapshot_";
  filename.fill('0');
  filename.width(4);
  filename << istep;
  filename << ".txt";
  std::cout << "Writing snapshot " << filename.str() << std::endl;
  std::ofstream ofile(filename.str().c_str());
  ofile << "# time: " << time << " s\n";
  ofile << "#\n";
  add_code_block_to_file(ofile);
  ofile << "# x (m)\trho (kg m^-3)\tu (m s^-1)\tP (kg m^-1 s^-2)\n";
  for (uint_fast32_t i = 1; i < ncell + 1; ++i) {
    ofile << cells[i]._midpoint << "\t" << cells[i]._rho << "\t" << cells[i]._u
          << "\t" << cells[i]._P << "\n";
  }
  ofile.close();
}

/**
 * @brief Round the given integer down to the nearest power of 2.
 *
 * @param x Integer.
 * @return Nearest lower power of 2.
 */
static inline uint_fast64_t round_power2_down(uint_fast64_t x) {
  --x;
  x |= (x >> 1);
  x |= (x >> 2);
  x |= (x >> 4);
  x |= (x >> 8);
  x |= (x >> 16);
  x |= (x >> 32);
  x >>= 1;
  ++x;
  return x;
}

/**
 * @brief Pair-wise slope limiter.
 *
 * @param phiL Original left state value.
 * @param phiR Original right state value.
 * @param phi_dash Unlimited extrapolated value in the left state.
 * @return Limited extrapolated value in the left state.
 */
static inline double pair_wise_limiter(const double phiL, const double phiR,
                                       const double phi_dash) {

  if (phiL == phiR) {

    return phiL;

  } else {

    const double phibar = phiL + 0.5 * (phiR - phiL);
    const double phidiff = std::abs(phiL - phiR);
    const double delta1 = 0.5 * phidiff;
    const double delta2 = 0.25 * phidiff;

    if (phiL < phiR) {

      const double phi_min = std::min(phiL, phiR);
      double phimin;
      if ((phi_min - delta1) * phi_min > 0.) {
        phimin = phi_min - delta1;
      } else {
        phimin = phi_min * std::abs(phi_min) / (std::abs(phi_min) + delta1);
      }
      return std::max(phimin, std::min(phibar + delta2, phi_dash));

    } else {

      const double phi_max = std::max(phiL, phiR);
      double phiplu;
      if ((phi_max + delta1) * phi_max > 0.) {
        phiplu = phi_max + delta1;
      } else {
        phiplu = phi_max * std::abs(phi_max) / (std::abs(phi_max) + delta1);
      }
      return std::min(phiplu, std::max(phibar - delta2, phi_dash));
    }
  }
}

/**
 * @brief Main simulation program.
 *
 * Usage: ./HydroCode1D [NUMBER OF THREADS]
 *
 * This program takes a single (optional) command line argument: the number of
 * shared memory parallel threads to use. If not set, the number of available
 * threads in the environment variable OMP_NUM_THREADS is used. Note that for
 * low cell numbers, using a lot of threads will slow the computation down
 * rather than speed it up; it is advisable to try different thread numbers to
 * find the optimal value for a given problem, resolution and hardware.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  // time the program
  Timer total_time;
  total_time.start();

  const unsigned int ncell = NCELL;

  // figure out how many threads we are using and tell the user about this
  int num_thread;
#pragma omp parallel
  {
#pragma omp single
    {
      num_thread = omp_get_num_threads();
      std::cout << num_thread << " thread(s) available." << std::endl;
    }
  }

  if (argc > 1) {
    const int num_thread_request = atoi(argv[1]);
    if (num_thread_request <= num_thread) {
      num_thread = num_thread_request;
    } else {
      std::cout << "More threads requested than available.\nWill fall back to "
                   "number of available threads."
                << std::endl;
    }
  }
  omp_set_num_threads(num_thread);
#pragma omp parallel
  {
#pragma omp single
    {
      num_thread = omp_get_num_threads();
      std::cout << "Will run using " << num_thread << " thread(s)."
                << std::endl;
    }
  }

#if EOS == EOS_ISOTHERMAL
  std::cout << "Isothermal sound speed: " << std::sqrt(ISOTHERMAL_C_SQUARED)
            << " m s^-1." << std::endl;
#endif

  // initialize the gravity solver
  gravity_init(ncell);

  // initialize the time stepping
  timeline_init();

  // create the 1D spherical grid
  // we create 2 ghost cells to the left and to the right of the simulation box
  // to handle boundary conditions
  Cell *cells = new Cell[ncell + 2];
#pragma omp parallel for
  for (uint_fast32_t i = 0; i < ncell + 2; ++i) {
    // cell positions (lower limit, center and upper limit) are precomputed for
    // maximal efficiency
    cells[i]._lowlim = RMIN + (i - 1.) * CELLSIZE;
    cells[i]._midpoint = RMIN + (i - 0.5) * CELLSIZE;
    cells[i]._uplim = RMIN + i * CELLSIZE;
    cells[i]._V = CELLSIZE;
#if DIMENSIONALITY == DIMENSIONALITY_1D
    cells[i]._V_real = cells[i]._V;
    cells[i]._V_real_half = 0.5 * cells[i]._V;
#elif DIMENSIONALITY == DIMENSIONALITY_2D
    cells[i]._V_real = M_PI * (cells[i]._uplim * cells[i]._uplim -
                               cells[i]._lowlim * cells[i]._lowlim);
    cells[i]._V_real_half = M_PI * (cells[i]._midpoint * cells[i]._midpoint -
                                    cells[i]._lowlim * cells[i]._lowlim);
#elif DIMENSIONALITY == DIMENSIONALITY_3D
    cells[i]._V_real = 4. * M_PI / 3. *
                       (cells[i]._uplim * cells[i]._uplim * cells[i]._uplim -
                        cells[i]._lowlim * cells[i]._lowlim * cells[i]._lowlim);
    cells[i]._V_real_half =
        4. * M_PI / 3. *
        (cells[i]._midpoint * cells[i]._midpoint * cells[i]._midpoint -
         cells[i]._lowlim * cells[i]._lowlim * cells[i]._lowlim);
#endif
  }

  // set up the initial condition
  // this bit is handled by IC.hpp, and user specific code in UserInput.hpp
  ic_initialize(cells, ncell);

  // convert to comoving variables if applicable
  timeline_do_variable_conversion(cells, ncell);

  // Courant factor for the CFL time step criterion
  // we use a very conservative value
  const double courant_factor = COURANT_FACTOR;

  std::cout << "Courant factor: " << courant_factor << std::endl;

  // convert the input primitive variables into conserved variables, and compute
  // the initial time step
  // we use a global time step, which is the minimum time step among all cells
  double min_physical_dt = timeline_get_max_physical_dt();
#pragma omp parallel for reduction(min : min_physical_dt)
  // convert primitive variables to conserved variables
  for (uint_fast32_t i = 1; i < ncell + 1; ++i) {
    // apply the equation of state to get the initial pressure (if necessary)
    initial_pressure(cells[i]);

    assert_condition(cells[i]._rho >= 0., "cells[%" PRIiFAST32 "]._rho = %g", i,
                     cells[i]._rho);
    assert_condition(cells[i]._P >= 0., "cells[%" PRIiFAST32 "]._P = %g", i,
                     cells[i]._P);

    // use the cell volume to convert primitive into conserved variables
    cells[i]._m = cells[i]._rho * cells[i]._V;
    cells[i]._p = cells[i]._m * cells[i]._u;
    cells[i]._E = cells[i]._P * cells[i]._V / (GAMMA - 1.) +
                  0.5 * cells[i]._u * cells[i]._p;

    assert_condition(cells[i]._m >= 0., "cells[%" PRIiFAST32 "]._m = %g", i,
                     cells[i]._m);
    assert_condition(cells[i]._E >= 0., "cells[%" PRIiFAST32 "]._E = %g", i,
                     cells[i]._E);

    // time step criterion
    // only non-vacuum cells are considered for the time step
    if (cells[i]._rho > 0.) {
      const double cs = std::sqrt(GAMMA * cells[i]._P / cells[i]._rho) +
                        std::abs(cells[i]._u);
      const double dt = courant_factor * cells[i]._V / cs;
      min_physical_dt = std::min(min_physical_dt, dt);
    }
  }

  // set the system time step based on the maximum allowed physical timestep
  timeline_set_timestep(min_physical_dt);

#if RIEMANNSOLVER_TYPE == RIEMANNSOLVER_TYPE_HLLC
  HLLCRiemannSolver solver(GAMMA);
#elif RIEMANNSOLVER_TYPE == RIEMANNSOLVER_TYPE_EXACT
  RiemannSolver solver(GAMMA);
#endif

  std::ofstream efile("energy.log");
  add_code_block_to_file(efile);
  efile << "# t (s)\tEkin (kg m^2 s^-2)\tEpot (kg m^2 s^-2)\tEtherm (kg m^2 "
           "s^-2)\tEtot (kg m^2 s^-2)\n";

  // initialize some variables used to guesstimate the remaining run time
  double last_stat_time = total_time.interval();
  const double stat_interval = 60.; // stat output every minute
  Timer step_time;
  double time_since_last = 0.;
  double time_since_start = 0.;
  unsigned int steps_since_last = 0;
  uint_fast64_t isnap = timeline_get_initial_snap_index();
  // main simulation loop: perform NSTEP steps
  while (timeline_next_step()) {

    // start the step timer
    step_time.start();

    // add the spherical source term. Handled by Spherical.hpp
    add_spherical_source_term(cells, ncell, timeline_get_system_source_dt());

    // do first gravity kick, handled by Potential.hpp
    do_gravity(cells, ncell, timeline_get_system_gravity_dt());

    // update the primitive variables based on the values of the conserved
    // variables and the current cell volume
    // also compute the new time step
    min_physical_dt = timeline_get_max_physical_dt();
    double Ekin_tot = 0., Epot_tot = 0., Etherm_tot = 0., Etot_tot = 0.;
#pragma omp parallel for reduction(min : min_physical_dt)                      \
    reduction(+ : Ekin_tot, Epot_tot, Etherm_tot, Etot_tot)
    for (uint_fast32_t i = 1; i < ncell + 1; ++i) {

      assert_condition(cells[i]._m >= 0., "cells[%" PRIiFAST32 "]._m = %g", i,
                       cells[i]._m);

      cells[i]._rho = cells[i]._m / cells[i]._V;
      if (cells[i]._m > 0.) {
        cells[i]._u = cells[i]._p / cells[i]._m;
      } else {
        cells[i]._u = 0.;
      }
      // the pressure update depends on the equation of state
      // this is handled in EOS.hpp
      update_pressure(cells[i]);
      // only non-vacuum cells are considered for the time step
      if (cells[i]._rho > 0.) {
        const double cs = std::sqrt(GAMMA * cells[i]._P / cells[i]._rho) +
                          std::abs(cells[i]._u);
        const double dt = courant_factor * cells[i]._V / cs;
        min_physical_dt = std::min(min_physical_dt, dt);
      }

      // energy statistics

      // first: gather some variables
      const double m = cells[i]._rho * cells[i]._V_real;

      if (m > 0.) {
        const double v2 = cells[i]._u * cells[i]._u;
        const double u = cells[i]._P / ((GAMMA - 1.) * cells[i]._rho);
        const double phi = cells[i]._pot;

        // now compute the energies
        const double Ekin = 0.5 * m * v2;
        const double Epot = m * phi;
        const double Etherm = m * u;
        const double Etot = Ekin + Epot + Etherm;

        // add to the totals
        Ekin_tot += Ekin;
        Epot_tot += Epot;
        Etherm_tot += Etherm;
        Etot_tot += Etot;
      }
    }

    timeline_set_timestep(min_physical_dt);

    // now output energy statistics
    const double t = timeline_get_current_physical_time();
    efile << t << "\t" << Ekin_tot << "\t" << Epot_tot << "\t" << Etherm_tot
          << "\t" << Etot_tot << "\n";

    // check if we need to output a snapshot
    if (timeline_do_write_snapshot()) {
      // write the actual snapshot
      write_snapshot(isnap, t, cells, ncell);
      ++isnap;
    }

    // check if we need to output run time statistics
    // (this is the only useful output the user gets to see before the
    //  simulation finishes)
    const double total_time_interval = total_time.interval();
    if (total_time_interval - last_stat_time > stat_interval) {
      last_stat_time = total_time_interval;
      // yes: display some statistics and a guesstimate of the remaining run
      // time
      const double pct = timeline_get_simulated_fraction();
      std::cout << get_timestamp() << "\t" << ncell << ": ";
      timeline_print_system_stats();
      const double avg_time_since_last = time_since_last / steps_since_last;
      std::cout << "\t\t\tAverage time per step: " << avg_time_since_last
                << " s" << std::endl;
      time_since_start += time_since_last;
      const double time_to_go = time_since_start * (100. - pct) / pct;
      std::cout << "\t\t\tEstimated time to go: " << time_to_go << " s"
                << std::endl;
      // reset guesstimate counters
      time_since_last = 0.;
      steps_since_last = 0;
    }

    // apply boundary conditions
    // handled by Boundaries.hpp
    boundary_conditions_primitive_variables(cells, ncell);

#if HYDRO_ORDER == HYDRO_ORDER_2
// compute slope limited gradients for the primitive variables in each cell
#pragma omp parallel for
    for (uint_fast32_t i = 1; i < ncell + 1; ++i) {
      const double dx = cells[i + 1]._midpoint - cells[i - 1]._midpoint;
      const double dx_inv = 1. / dx;
      const double half_dx = 0.5 * dx;

      const double gradrho = (cells[i + 1]._rho - cells[i - 1]._rho) * dx_inv;
      const double rhomax = std::max(cells[i - 1]._rho, cells[i + 1]._rho);
      const double rhomin = std::min(cells[i - 1]._rho, cells[i + 1]._rho);
      const double rho_ext_plu = half_dx * gradrho;
      const double rho_ext_min = -half_dx * gradrho;
      const double rhoextmax = std::max(rho_ext_min, rho_ext_plu);
      const double rhoextmin = std::min(rho_ext_min, rho_ext_plu);
      const double alpha_rho =
          (gradrho != 0.)
              ? std::min(1.,
                         0.5 * std::min((rhomax - cells[i]._rho) / rhoextmax,
                                        (rhomin - cells[i]._rho) / rhoextmin))
              : 1.;
      cells[i]._grad_rho = alpha_rho * gradrho;

      const double gradu = (cells[i + 1]._u - cells[i - 1]._u) * dx_inv;
      const double umax = std::max(cells[i - 1]._u, cells[i + 1]._u);
      const double umin = std::min(cells[i - 1]._u, cells[i + 1]._u);
      const double u_ext_plu = half_dx * gradu;
      const double u_ext_min = -half_dx * gradu;
      const double uextmax = std::max(u_ext_min, u_ext_plu);
      const double uextmin = std::min(u_ext_min, u_ext_plu);
      const double alpha_u =
          (gradu != 0.)
              ? std::min(1., 0.5 * std::min((umax - cells[i]._u) / uextmax,
                                            (umin - cells[i]._u) / uextmin))
              : 1.;
      cells[i]._grad_u = alpha_u * gradu;

      const double gradP = (cells[i + 1]._P - cells[i - 1]._P) * dx_inv;
      const double Pmax = std::max(cells[i - 1]._P, cells[i + 1]._P);
      const double Pmin = std::min(cells[i - 1]._P, cells[i + 1]._P);
      const double P_ext_plu = half_dx * gradP;
      const double P_ext_min = -half_dx * gradP;
      const double Pextmax = std::max(P_ext_min, P_ext_plu);
      const double Pextmin = std::min(P_ext_min, P_ext_plu);
      const double alpha_P =
          (gradP != 0.)
              ? std::min(1., 0.5 * std::min((Pmax - cells[i]._P) / Pextmax,
                                            (Pmin - cells[i]._P) / Pextmin))
              : 1.;
      cells[i]._grad_P = alpha_P * gradP;

      assert_condition(cells[i]._grad_rho == cells[i]._grad_rho,
                       "cells[%" PRIiFAST32 "]._grad_rho = NaN", i);
      assert_condition(cells[i]._grad_u == cells[i]._grad_u,
                       "cells[%" PRIiFAST32 "]._grad_u = NaN", i);
      assert_condition(cells[i]._grad_P == cells[i]._grad_P,
                       "cells[%" PRIiFAST32 "]._grad_P = NaN", i);
    }

    // apply boundary conditions for the gradients
    // handled by Boundaries.hpp
    boundary_conditions_gradients(cells, ncell);

// evolve all primitive variables forward in time for half a time step
// using the Euler equations and the spatial gradients within the cells
#pragma omp parallel for
    for (uint_fast32_t i = 0; i < ncell + 2; ++i) {
      const double half_dt = 0.5 * timeline_get_system_hydro_dt();
      const double rho = cells[i]._rho;
      const double u = cells[i]._u;
      const double P = cells[i]._P;
      cells[i]._rho -=
          half_dt * (rho * cells[i]._grad_u + u * cells[i]._grad_rho);
      if (rho > 0.) {
        cells[i]._u -=
            half_dt * (u * cells[i]._grad_u + cells[i]._grad_P / rho);
      }
      cells[i]._P -=
          half_dt * (GAMMA * P * cells[i]._grad_u + u * cells[i]._grad_P);

      if (cells[i]._rho < 0.) {
        cells[i]._rho = rho;
      }
      if (cells[i]._P < 0.) {
        cells[i]._P = P;
      }

      // add gravity prediction. Handled by Potential.hpp.
      add_gravitational_prediction(cells[i],
                                   0.5 * timeline_get_system_gravity_dt());
    }

#elif HYDRO_ORDER == HYDRO_ORDER_1

// set all gradients to zero to disable the second order scheme
#pragma omp parallel for
    for (uint_fast32_t i = 0; i < ncell + 2; ++i) {
      cells[i]._grad_rho = 0.;
      cells[i]._grad_u = 0.;
      cells[i]._grad_P = 0.;
    }

#endif

// do the flux exchange
// to avoid thread concurrency, we compute each flux twice, and only update one
// cell at a time
#pragma omp parallel for
    for (uint_fast32_t i = 1; i < ncell + 1; ++i) {
      const double dt = timeline_get_system_hydro_dt();
      // left flux
      {
        // get the variables in the left and right state
        const double rhoL = cells[i - 1]._rho;
        const double uL = cells[i - 1]._u;
        const double PL = cells[i - 1]._P;
        const double rhoR = cells[i]._rho;
        const double uR = cells[i]._u;
        const double PR = cells[i]._P;

        assert_condition(rhoL == rhoL, "cells[%" PRIiFAST32 "]._rho = NaN",
                         i - 1);
        assert_condition(rhoL >= 0., "cells[%" PRIiFAST32 "]._rho = %g", i - 1,
                         rhoL);
        assert_condition(uL == uL, "cells[%" PRIiFAST32 "]._u = NaN", i - 1);
        assert_condition(PL == PL, "cells[%" PRIiFAST32 "]._P = NaN", i - 1);
        assert_condition(PL >= 0., "cells[%" PRIiFAST32 "]._P = %g", i - 1, PL);

        assert_condition(rhoR == rhoR, "cells[%" PRIiFAST32 "]._rho = NaN", i);
        assert_condition(rhoR >= 0., "cells[%" PRIiFAST32 "]._rho = %g", i,
                         rhoR);
        assert_condition(uR == uR, "cells[%" PRIiFAST32 "]._u = NaN", i);
        assert_condition(PR == PR, "cells[%" PRIiFAST32 "]._P = NaN", i);
        assert_condition(PR >= 0., "cells[%" PRIiFAST32 "]._P = %g", i, PR);

        // do the second order spatial reconstruction
        const double dmin = 0.5 * (cells[i]._midpoint - cells[i - 1]._midpoint);
        const double dplu = -dmin;
        double rhoL_dash, rhoR_dash, uL_dash, uR_dash, PL_dash, PR_dash;
        rhoL_dash = rhoL + dmin * cells[i - 1]._grad_rho;
        uL_dash = uL + dmin * cells[i - 1]._grad_u;
        PL_dash = PL + dmin * cells[i - 1]._grad_P;
        rhoR_dash = rhoR + dplu * cells[i]._grad_rho;
        uR_dash = uR + dplu * cells[i]._grad_u;
        PR_dash = PR + dplu * cells[i]._grad_P;

        if (rhoL_dash < 0.) {
          rhoL_dash = rhoL;
        }
        if (rhoR_dash < 0.) {
          rhoR_dash = rhoR;
        }
        if (PL_dash < 0.) {
          PL_dash = PL;
        }
        if (PR_dash < 0.) {
          PR_dash = PR;
        }

        // pair-wise limiter
        rhoL_dash = pair_wise_limiter(rhoL, rhoR, rhoL_dash);
        rhoR_dash = pair_wise_limiter(rhoR, rhoL, rhoR_dash);
        uL_dash = pair_wise_limiter(uL, uR, uL_dash);
        uR_dash = pair_wise_limiter(uR, uL, uR_dash);
        PL_dash = pair_wise_limiter(PL, PR, PL_dash);
        PR_dash = pair_wise_limiter(PR, PL, PR_dash);

        // solve the Riemann problem at the interface between the two cells
        double mflux, pflux, Eflux;
        solver.solve_for_flux(rhoL_dash, uL_dash, PL_dash, rhoR_dash, uR_dash,
                              PR_dash, mflux, pflux, Eflux);

        cells[i]._m += dt * mflux;
        cells[i]._p += dt * pflux;
        cells[i]._E += dt * Eflux;

        assert_condition(cells[i]._m >= 0., "cells[%" PRIiFAST32 "]._m = %g", i,
                         cells[i]._m);
        assert_condition(cells[i]._E >= 0., "cells[%" PRIiFAST32 "]._E = %g", i,
                         cells[i]._E);

        // make sure mass and energy stay positive
        //        cells[i]._m = std::max(cells[i]._m, 0.);
        //        cells[i]._E = std::max(cells[i]._E, 0.);
      }
      // right flux
      {
        // get the variables in the left and right state
        const double rhoL = cells[i]._rho;
        const double uL = cells[i]._u;
        const double PL = cells[i]._P;
        const double rhoR = cells[i + 1]._rho;
        const double uR = cells[i + 1]._u;
        const double PR = cells[i + 1]._P;

        assert_condition(rhoL == rhoL, "cells[%" PRIiFAST32 "]._rho = NaN", i);
        assert_condition(rhoL >= 0., "cells[%" PRIiFAST32 "]._rho = %g", i,
                         rhoL);
        assert_condition(uL == uL, "cells[%" PRIiFAST32 "]._u = NaN", i);
        assert_condition(PL == PL, "cells[%" PRIiFAST32 "]._P = NaN", i);
        assert_condition(PL >= 0., "cells[%" PRIiFAST32 "]._P = %g", i, PL);

        assert_condition(rhoR == rhoR, "cells[%" PRIiFAST32 "]._rho = NaN",
                         i + 1);
        assert_condition(rhoR >= 0., "cells[%" PRIiFAST32 "]._rho = %g", i + 1,
                         rhoR);
        assert_condition(uR == uR, "cells[%" PRIiFAST32 "]._u = NaN", i + 1);
        assert_condition(PR == PR, "cells[%" PRIiFAST32 "]._P = NaN", i + 1);
        assert_condition(PR >= 0., "cells[%" PRIiFAST32 "]._P = %g", i + 1, PR);

        // do the second order spatial reconstruction
        const double dmin = 0.5 * (cells[i + 1]._midpoint - cells[i]._midpoint);
        const double dplu = -dmin;
        double rhoL_dash, rhoR_dash, uL_dash, uR_dash, PL_dash, PR_dash;
        rhoL_dash = rhoL + dmin * cells[i]._grad_rho;
        uL_dash = uL + dmin * cells[i]._grad_u;
        PL_dash = PL + dmin * cells[i]._grad_P;
        rhoR_dash = rhoR + dplu * cells[i + 1]._grad_rho;
        uR_dash = uR + dplu * cells[i + 1]._grad_u;
        PR_dash = PR + dplu * cells[i + 1]._grad_P;

        if (rhoL_dash < 0.) {
          rhoL_dash = rhoL;
        }
        if (rhoR_dash < 0.) {
          rhoR_dash = rhoR;
        }
        if (PL_dash < 0.) {
          PL_dash = PL;
        }
        if (PR_dash < 0.) {
          PR_dash = PR;
        }

        // pair-wise limiter
        rhoL_dash = pair_wise_limiter(rhoL, rhoR, rhoL_dash);
        rhoR_dash = pair_wise_limiter(rhoR, rhoL, rhoR_dash);
        uL_dash = pair_wise_limiter(uL, uR, uL_dash);
        uR_dash = pair_wise_limiter(uR, uL, uR_dash);
        PL_dash = pair_wise_limiter(PL, PR, PL_dash);
        PR_dash = pair_wise_limiter(PR, PL, PR_dash);

        // solve the Riemann problem at the interface between the two cells
        double mflux, pflux, Eflux;
        solver.solve_for_flux(rhoL_dash, uL_dash, PL_dash, rhoR_dash, uR_dash,
                              PR_dash, mflux, pflux, Eflux);

        cells[i]._m -= dt * mflux;
        cells[i]._p -= dt * pflux;
        cells[i]._E -= dt * Eflux;

        assert_condition(cells[i]._m >= 0., "cells[%" PRIiFAST32 "]._m = %g", i,
                         cells[i]._m);
        assert_condition(cells[i]._E >= 0., "cells[%" PRIiFAST32 "]._E = %g", i,
                         cells[i]._E);

        // make sure mass and energy stay positive
        //        cells[i]._m = std::max(cells[i]._m, 0.);
        //        cells[i]._E = std::max(cells[i]._E, 0.);
      }
    }

    // add the spherical source term
    // handled by Spherical.hpp
    add_spherical_source_term(cells, ncell, timeline_get_system_source_dt());

    // do the second gravity kick
    // handled by Potential.hpp
    do_gravity(cells, ncell, timeline_get_system_gravity_dt());

    // stop the step timer, and update guesstimate counters
    step_time.stop();
    time_since_last += step_time.value();
    ++steps_since_last;
    step_time.reset();

    // update the system time
    timeline_do_timestep();
  }

  // write the final snapshots
  write_snapshot(isnap, timeline_get_current_physical_time(), cells, ncell);

  // clean up the gravity solver
  free_gravity();

  // clean up: free cell memory
  delete[] cells;

  // stop timing the program and display run time information
  total_time.stop();
  std::cout << "Total program time: " << total_time.value() << " s."
            << std::endl;

  // all went well: return with exit code 0
  return 0;
}
