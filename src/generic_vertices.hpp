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
      return is_zero(val);
   }

private:
   std::complex<double> val;
};

/**
 * @class SSSSVertex
 */
class SSSSVertex {
public:
   explicit SSSSVertex(std::complex<double> v)
      : val(v) {}

   std::complex<double> value() const { return val; }

   bool isZero() const {
      return is_zero(val);
   }

private:
   std::complex<double> val;
};

/**
 * @class FFSVertex
 */
class FFSVertex {
public:
   FFSVertex(std::complex<double> left,
             std::complex<double> right)
      : value(left, right) {}

   std::complex<double> left() const { return value.first; }
   std::complex<double> right() const { return value.second; }

   bool isZero() const {
      return is_zero(value.first) && is_zero(value.second);
   }

private:
   std::pair<std::complex<double>, std::complex<double>> value;
};

/**
 * @class FFVVertex
 */
class FFVVertex {
public:
   FFVVertex(std::complex<double> left,
             std::complex<double> right)
      : value(left, right) {}

   std::complex<double> left() const { return value.first; }
   std::complex<double> right() const { return value.second; }

   bool isZero() const {
      return is_zero(value.first) && is_zero(value.second);
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
         "SSVVertex: Wrong index combination" );
      return 0.0;
   }

   bool isZero() const {
      return is_zero(val);
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
      return is_zero(val);
   }

private:
   std::complex<double> val;
};

/**
 * @class SSVVVertex
 */
class SSVVVertex {
public:
   explicit SSVVVertex(std::complex<double> v)
      : val(v) {}

   std::complex<double> value() const { return val; }

   bool isZero() const {
      return is_zero(val);
   }

private:
   std::complex<double> val;
};

/**
 * @class VVVVertex
 */
class VVVVertex {
public:
   explicit VVVVertex(std::complex<double> v)
      : val(v) {}

   std::complex<double> value() const { return val; }

   bool isZero() const {
      return is_zero(val);
   }

private:
   std::complex<double> val;
};

/**
 * @class VVVVVertex
 */
class VVVVVertex {
public:
   VVVVVertex(std::complex<double> v1, std::complex<double> v2,
              std::complex<double> v3)
      : coeff_12_34(v1), coeff_13_24(v2), coeff_14_23(v3) {}

   bool isZero() const {
      return is_zero(coeff_12_34) && is_zero(coeff_13_24) &&
         is_zero(coeff_14_23);
   }

private:
   std::complex<double> coeff_12_34;
   std::complex<double> coeff_13_24;
   std::complex<double> coeff_14_23;
};

/**
 * @class GGSVertex
 */
class GGSVertex {
public:
   explicit GGSVertex(std::complex<double> v)
      : val(v) {}

   std::complex<double> value() const { return val; }

   bool isZero() const {
      return is_zero(val);
   }

private:
   std::complex<double> val;
};

/**
 * @class GGVVertex
 */
class GGVVertex {
public:
   explicit GGVVertex(std::complex<double> v)
      : val(v) {}

   std::complex<double> value() const { return val; }

   bool isZero() const {
      return is_zero(val);
   }

private:
   std::complex<double> val;
};

} // namespace cxx_qft

} // namespace flexiblesusy

#endif
