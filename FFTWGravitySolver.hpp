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
 * @file FFTWGravitySolver.hpp
 *
 * @brief 1D FFTW Poisson solver for gravity.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */

#include <cinttypes>
#include <fftw3.h>

/**
 * @brief 1D FFTW Poisson solver for gravity.
 */
class FFTWGravitySolver {
private:
  /*! @brief Array used to perform the fast Fourier transforms. */
  fftw_complex *_fftw_array;

  /*! @brief FFTW plan for the forward transform. */
  fftw_plan _forward_plan;

  /*! @brief FFTW plan for the backward transform. */
  fftw_plan _backward_plan;

  /*! @brief Value of Newton's gravitational constant. */
  const double _G;

  /*! @brief Number of cells. */
  const uint_fast32_t _ncell;

  /*! @brief Length of the box (in m). */
  const double _box_length;

public:
  /**
   * @brief Constructor.
   *
   * @param ncell Number of cells in the grid.
   * @param box_length Length of the box (in m).
   * @param G Value for Newton's gravitational constant (in N m^2 kg^-2).
   */
  inline FFTWGravitySolver(const uint_fast32_t ncell, const double box_length,
                           const double G = 6.674e-11)
      : _G(G), _ncell(ncell), _box_length(box_length) {

    _fftw_array = reinterpret_cast<fftw_complex *>(
        fftw_malloc(sizeof(fftw_complex) * ncell));

    _forward_plan = fftw_plan_dft_1d(ncell, _fftw_array, _fftw_array,
                                     FFTW_FORWARD, FFTW_MEASURE);
    _backward_plan = fftw_plan_dft_1d(ncell, _fftw_array, _fftw_array,
                                      FFTW_BACKWARD, FFTW_MEASURE);
  }

  /**
   * @brief Destructor.
   */
  inline ~FFTWGravitySolver() {
    fftw_destroy_plan(_backward_plan);
    fftw_destroy_plan(_forward_plan);
    fftw_free(_fftw_array);
  }

  /**
   * @brief Compute the accelerations for the given density field.
   *
   * @param rho Input density field (in kg m^-3).
   * @param a Output acceleration (in m s^-2).
   */
  inline void compute_accelerations(const double *rho, double *a) {

    // copy density into FFTW array
    for (uint_fast32_t i = 0; i < _ncell; ++i) {
      _fftw_array[i][0] = rho[i];
      _fftw_array[i][1] = 0.;
    }

    // apply Fourier transform
    fftw_execute(_forward_plan);

    // compute potential in Fourier space
    // the k = 0 mode is discarded
    _fftw_array[0][0] = 0.;
    _fftw_array[0][1] = 0.;
    // positive k modes
    for (uint_fast32_t i = 1; i <= _ncell / 2; ++i) {
      const double k = 2. * M_PI * i / _box_length;
      const double kfac = -4. * M_PI * _G / (k * k);
      _fftw_array[i][0] *= kfac;
      _fftw_array[i][1] *= kfac;
    }
    // negative k modes
    for (uint_fast32_t i = _ncell / 2 + 1; i < _ncell; ++i) {
      const double k = 2. * M_PI * (i - _ncell) / _box_length;
      const double kfac = -4. * M_PI * _G / (k * k);
      _fftw_array[i][0] *= kfac;
      _fftw_array[i][1] *= kfac;
    }

    // inverse Fourier transform
    fftw_execute(_backward_plan);

    // compute accelerations
    // we factor in the fact that applying forward + backward transform results
    // in an additional factor _ncell
    // no idea what happens to the factor 1/2...
    const double dxinv = 1. / _box_length;
    double phi_m, phi_p;
    // first element: apply periodic boundary
    phi_m = _fftw_array[_ncell - 1][0];
    phi_p = _fftw_array[1][0];
    a[0] = -(phi_p - phi_m) * dxinv;
    // other elements (except last)
    for (uint_fast32_t i = 1; i < _ncell - 1; ++i) {
      phi_m = _fftw_array[i - 1][0];
      phi_p = _fftw_array[i + 1][0];
      a[i] = -(phi_p - phi_m) * dxinv;
    }
    // last element
    phi_m = _fftw_array[_ncell - 2][0];
    phi_p = _fftw_array[0][0];
    a[_ncell - 1] = -(phi_p - phi_m) * dxinv;
  }
};
