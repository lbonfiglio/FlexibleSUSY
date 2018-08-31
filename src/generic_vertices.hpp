// ====================================================================
// This file is part of FlexibleSUSY.
//
// FlexibleSUSY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// FlexibleSUSY is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FlexibleSUSY.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

#ifndef GENERIC_VERTICES_H
#define GENERIC_VERTICES_H

#include "error.hpp"
#include "numerics2.hpp"

#include <complex>

namespace flexiblesusy {

namespace cxx_qft {

/**
 * @class SSSVertex
 */
class SSSVertex {
public:
   explicit SSSVertex(std::complex<double> v)
      : val(v) {}

   std::complex<double> value() const { return val; }

   bool isZero() const {
      return (is_zero(val.real()) && is_zero(val.imag()));
   }

private:
   std::complex<double> val;
};

/**
 * @class FFSVertex
 */
class FFSVertex {
public:
   FFSVertex(const std::complex<double>& left,
                const std::complex<double>& right)
      : value(left, right) {}

   std::complex<double> left() const { return value.first; }
   std::complex<double> right() const { return value.second; }

   bool isZero() const {
      return (is_zero(value.first.real()) && is_zero(value.first.imag()) &&
              is_zero(value.second.real()) && is_zero(value.second.imag()));
   }

private:
   std::pair<std::complex<double>, std::complex<double>> value;
};

/**
 * @class SSVVertex
 */
class SSVVertex {
public:
   SSVVertex(std::complex<double> v, int mi, int si )
      : val(v), minuendIndex(mi), subtrahendIndex(si) {}

   std::complex<double> value(int mi, int si) const {
      if( mi == minuendIndex && si == subtrahendIndex ) {
         return val;
      }

      if( mi == subtrahendIndex && si == minuendIndex ) {
         return -val;
      }

      throw std::invalid_argument(
         "MomentumDifferenceVertex: Wrong index combination" );
      return 0.0;
   }

   bool isZero() const {
      return (is_zero(val.real()) && is_zero(val.imag()));
   }

private:
  std::complex<double> val;
  int minuendIndex;
  int subtrahendIndex;
};

/**
 * @class SVVVertex
 */
class SVVVertex {
public:
   explicit SVVVertex(std::complex<double> v)
      : val(v) {}

   std::complex<double> value() const { return val; }

   bool isZero() const {
      return (is_zero(val.real()) && is_zero(val.imag()));
   }

private:
  std::complex<double> val;
};

} // namespace cxx_qft

} // namespace flexiblesusy

#endif
