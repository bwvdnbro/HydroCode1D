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
 * @file Cosmology.hpp
 *
 * @brief Cosmological integration routines.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef COSMOLOGY_HPP
#define COSMOLOGY_HPP

/*! @brief Number of elements in the cosmological interpolation tables. */
#define COSMOLOGY_NTAB 10000

#include "PhysicalConstants.hpp"

#include <cinttypes>
#include <cmath>
#include <gsl/gsl_integration.h>

/**
 * @brief Cosmological integration routines.
 */
class Cosmology {
private:
  /*! @brief Omega matter. */
  const double _Omega_m;

  /*! @brief Omega radiation. */
  const double _Omega_r;

  /*! @brief Omega curvature. */
  const double _Omega_k;

  /*! @brief Omega dark energy. */
  const double _Omega_Lambda;

  /*! @brief Dark energy EoS constant. */
  const double _w_0;

  /*! @brief Dark energy EoS factor. */
  const double _w_a;

  /*! @brief Hubble constant. */
  const double _H0;

  /*! @brief Natural logarithm of the minimum scale factor. */
  const double _log_amin;

  /*! @brief Natural logarithm of the maximum scale factor. */
  const double _log_amax;

  /*! @brief Current scale factor. */
  double _a_current;

  /*! @brief @f$\frac{1}{a}@f$ factor interpolation table. */
  double *_table_ainv;

  /*! @brief @f$\frac{1}{a^2}@f$ factor interpolation table. */
  double *_table_a2inv;

  /*! @brief @f$a@f$ factor interpolation table. */
  double *_table_a;

  /*! @brief Time interpolation table. */
  double *_table_time;

  /*! @brief Offset of the time interpolation table. */
  double _time_offset;

  /**
   * @brief Return the linearly interpolated table value corresponding to the
   * given value.
   *
   * @param table Table to interpolate on.
   * @param x Interpolation point.
   * @param xmin Minimum value in the table.
   * @param xmax Maximum value in the table.
   * @param size Size of the table.
   * @return Value at the interpolation point.
   */
  static inline double interpolate(const double *table, const double x,
                                   const double xmin, const double xmax,
                                   const int_fast32_t size) {

    const double xx = ((x - xmin) / (xmax - xmin)) * size;
    const int_fast32_t i = static_cast<int_fast32_t>(xx);
    const int_fast32_t ii = (i >= size) ? size - 1 : i;

    if (ii <= 1) {
      return table[0] * xx;
    } else {
      return table[ii - 1] + (table[ii] - table[ii - 1]) * (xx - ii);
    }
  }

  /**
   * @brief Get the value of the dimensionless Hubble evolution function for
   * the given scale factor.
   *
   * @param a Scale factor.
   * @return Dimensionless Hubble evolution function.
   */
  inline double E(const double a) const {

    const double ainv = 1. / a;
    const double ainv2 = ainv * ainv;

    const double wtilda = (a - 1.) * _w_a - (1. + _w_0 + _w_a) * std::log(a);

    return std::sqrt(_Omega_m * ainv * ainv2 + _Omega_r * ainv2 * ainv2 +
                     _Omega_k * ainv2 + _Omega_Lambda * std::exp(3. * wtilda));
  }

  /**
   * @brief Get the Hubble constant for the given scale factor.
   *
   * @param a Scale factor.
   * @return Hubble constant (in s^-1).
   */
  inline double H(const double a) const { return _H0 * E(a); }

  /**
   * @brief Get the integrand of the @f$\frac{1}{a}@f$ factor.
   *
   * @param a Scale factor.
   * @param param Void-cast Cosmology pointer.
   * @return Integrand value.
   */
  inline static double integrand_ainv(const double a, void *param) {

    const Cosmology *c = reinterpret_cast<Cosmology *>(param);

    const double H = c->H(a);
    const double ainv = 1. / a;
    return (1. / H) * ainv * ainv;
  }

  /**
   * @brief Get the integrand of the @f$\frac{1}{a^2}@f$ factor.
   *
   * @param a Scale factor.
   * @param param Void-cast Cosmology pointer.
   * @return Integrand value.
   */
  inline static double integrand_a2inv(const double a, void *param) {

    const Cosmology *c = reinterpret_cast<Cosmology *>(param);

    const double H = c->H(a);
    const double ainv = 1. / a;
    return (1. / H) * ainv * ainv * ainv;
  }

  /**
   * @brief Get the integrand of the @f$a@f$ factor.
   *
   * @param a Scale factor.
   * @param param Void-cast Cosmology pointer.
   * @return Integrand value.
   */
  inline static double integrand_a(const double a, void *param) {

    const Cosmology *c = reinterpret_cast<Cosmology *>(param);

    const double H = c->H(a);
    return (1. / H);
  }

  /**
   * @brief Get the time integrand.
   *
   * @param a Scale factor.
   * @param param Void-cast Cosmology pointer.
   * @return Integrand value.
   */
  inline static double integrand_time(const double a, void *param) {

    const Cosmology *c = reinterpret_cast<Cosmology *>(param);

    const double H = c->H(a);
    return (1. / H) * (1. / a);
  }

public:
  /**
   * @brief Constructor.
   *
   * @param Omega_m Omega matter.
   * @param Omega_r Omega radiation.
   * @param Omega_k Omega curvature.
   * @param Omega_Lambda Omega dark energy.
   * @param w_0 Dark energy EoS constant.
   * @param w_a Dark energy EoS factor.
   * @param h Reduced Hubble constant.
   * @param amin Minimum scale factor.
   * @param amax Maximum scale factor.
   */
  inline Cosmology(const double Omega_m, const double Omega_r,
                   const double Omega_k, const double Omega_Lambda,
                   const double w_0, const double w_a, const double h,
                   const double amin, const double amax)
      : _Omega_m(Omega_m), _Omega_r(Omega_r), _Omega_k(Omega_k),
        _Omega_Lambda(Omega_Lambda), _w_0(w_0), _w_a(w_a),
        _H0(h * HUBBLE_IN_SI), _log_amin(std::log(amin)),
        _log_amax(std::log(amax)), _a_current(amin) {

    _table_ainv = new double[COSMOLOGY_NTAB];
    _table_a2inv = new double[COSMOLOGY_NTAB];
    _table_a = new double[COSMOLOGY_NTAB];
    _table_time = new double[COSMOLOGY_NTAB];

    const double delta_a = (_log_amax - _log_amin) / COSMOLOGY_NTAB;
    double *atable = new double[COSMOLOGY_NTAB];
    for (uint_fast32_t i = 0; i < COSMOLOGY_NTAB; ++i) {
      atable[i] = std::exp(_log_amin + delta_a * (i + 1.));
    }

    gsl_integration_workspace *gsl_ws =
        gsl_integration_workspace_alloc(COSMOLOGY_NTAB);

    double result, abserr;

    gsl_function gsl_F = {&integrand_ainv, this};
    for (uint_fast32_t i = 0; i < COSMOLOGY_NTAB; ++i) {
      gsl_integration_qag(&gsl_F, amin, atable[i], 0., 1.e-10, COSMOLOGY_NTAB,
                          GSL_INTEG_GAUSS61, gsl_ws, &result, &abserr);
      _table_ainv[i] = result;
    }

    gsl_F.function = &integrand_a2inv;
    for (uint_fast32_t i = 0; i < COSMOLOGY_NTAB; ++i) {
      gsl_integration_qag(&gsl_F, amin, atable[i], 0., 1.e-10, COSMOLOGY_NTAB,
                          GSL_INTEG_GAUSS61, gsl_ws, &result, &abserr);
      _table_a2inv[i] = result;
    }

    gsl_F.function = &integrand_a;
    for (uint_fast32_t i = 0; i < COSMOLOGY_NTAB; ++i) {
      gsl_integration_qag(&gsl_F, amin, atable[i], 0., 1.e-10, COSMOLOGY_NTAB,
                          GSL_INTEG_GAUSS61, gsl_ws, &result, &abserr);
      _table_a[i] = result;
    }

    gsl_F.function = &integrand_time;
    for (uint_fast32_t i = 0; i < COSMOLOGY_NTAB; ++i) {
      gsl_integration_qag(&gsl_F, amin, atable[i], 0., 1.e-10, COSMOLOGY_NTAB,
                          GSL_INTEG_GAUSS61, gsl_ws, &result, &abserr);
      _table_time[i] = result;
    }
    gsl_integration_qag(&gsl_F, 0., amin, 0., 1.e-10, COSMOLOGY_NTAB,
                        GSL_INTEG_GAUSS61, gsl_ws, &result, &abserr);
    _time_offset = result;

    gsl_integration_workspace_free(gsl_ws);
    delete[] atable;
  }

  /**
   * @brief Destructor.
   */
  inline ~Cosmology() {
    delete[] _table_ainv;
    delete[] _table_a2inv;
    delete[] _table_a;
    delete[] _table_time;
  }

  /**
   * @brief Get the @f$\frac{1}{a}@f$ factor for the given scale factor step.
   *
   * @param astart Scale factor at the start of the step.
   * @param astop Scale factor at the end of the step.
   * @return Integrated factor during the step.
   */
  inline double get_factor_ainv(const double astart, const double astop) const {

    return interpolate(_table_ainv, std::log(astop), _log_amin, _log_amax,
                       COSMOLOGY_NTAB) -
           interpolate(_table_ainv, std::log(astart), _log_amin, _log_amax,
                       COSMOLOGY_NTAB);
  }

  /**
   * @brief Get the @f$\frac{1}{a^2}@f$ factor for the given scale factor step.
   *
   * @param astart Scale factor at the start of the step.
   * @param astop Scale factor at the end of the step.
   * @return Integrated factor during the step.
   */
  inline double get_factor_a2inv(const double astart,
                                 const double astop) const {

    return interpolate(_table_a2inv, std::log(astop), _log_amin, _log_amax,
                       COSMOLOGY_NTAB) -
           interpolate(_table_a2inv, std::log(astart), _log_amin, _log_amax,
                       COSMOLOGY_NTAB);
  }

  /**
   * @brief Get the @f$a@f$ factor for the given scale factor step.
   *
   * @param astart Scale factor at the start of the step.
   * @param astop Scale factor at the end of the step.
   * @return Integrated factor during the step.
   */
  inline double get_factor_a(const double astart, const double astop) const {

    return interpolate(_table_a, std::log(astop), _log_amin, _log_amax,
                       COSMOLOGY_NTAB) -
           interpolate(_table_a, std::log(astart), _log_amin, _log_amax,
                       COSMOLOGY_NTAB);
  }

  /**
   * @brief Get the time interval corresponding to the given scale factor step.
   *
   * @param astart Scale factor at the start of the step.
   * @param astop Scale factor at the end of the step.
   * @return Corresponding time interval (in s).
   */
  inline double get_time_interval(const double astart,
                                  const double astop) const {

    return interpolate(_table_time, std::log(astop), _log_amin, _log_amax,
                       COSMOLOGY_NTAB) -
           interpolate(_table_time, std::log(astart), _log_amin, _log_amax,
                       COSMOLOGY_NTAB);
  }

  /**
   * @brief Get the scale factor interval corresponding to the given time
   * interval at the given scale factor.
   *
   * @param a Current value of the scale factor.
   * @param dt Time interval (in s).
   * @return Corresponding scale factor interval.
   */
  inline double get_scale_factor_interval(const double a,
                                          const double dt) const {

    return H(a) * dt;
  }

  /**
   * @brief Get the time since the Big Bang for the given scale factor.
   *
   * @param a Scale factor.
   * @return Time since Big Bang (in s).
   */
  inline double get_time_since_big_bang(const double a) const {

    return _time_offset + interpolate(_table_time, std::log(a), _log_amin,
                                      _log_amax, COSMOLOGY_NTAB);
  }
};

#endif // COSMOLOGY_HPP
