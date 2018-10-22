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

#include <array>
#include <complex>

namespace flexiblesusy {

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

   std::complex<double> value(int mi, int si) const;

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
   /**
    * Set coefficient of kinematic vector when
    * vertex is written in the order
    * V[i] V[j] V[k].
    */
   VVVVertex(int i, int j, int k, std::complex<double> c)
      : vec_0_id(i), vec_1_id(j), vec_2_id(k), val(c)
      {}
   ~VVVVertex() = default;
   VVVVertex(const VVVVertex&) = default;
   VVVVertex(VVVVertex&&) = default;
   VVVVertex& operator=(const VVVVertex&) = default;
   VVVVertex& operator=(VVVVertex&&) = default;

   std::complex<double> value(int, int, int) const;

   bool isZero() const {
      return is_zero(val);
   }

private:
   int vec_0_id;
   int vec_1_id;
   int vec_2_id;
   std::complex<double> val;

   int permutation_sign(int, int, int) const;
};

/**
 * @class VVVVVertex
 */
class VVVVVertex {
public:
   /**
    * Set coefficients of elements of kinematic vector
    * when vertex is written in the order
    * V[i] V[j] V[k] V[l].
    */
   VVVVVertex(int i, int j, int k, int l,
              std::complex<double> c1, std::complex<double> c2,
              std::complex<double> c3);
   ~VVVVVertex() = default;
   VVVVVertex(const VVVVVertex&) = default;
   VVVVVertex(VVVVVertex&&) = default;
   VVVVVertex& operator=(const VVVVVertex&) = default;
   VVVVVertex& operator=(VVVVVertex&&) = default;

   std::complex<double> value_1(int, int, int, int) const;
   std::complex<double> value_2(int, int, int, int) const;
   std::complex<double> value_3(int, int, int, int) const;

   bool isZero() const {
      return is_zero(coefficients[0]) && is_zero(coefficients[1]) &&
         is_zero(coefficients[2]);
   }

private:
   std::array<int, 4> indices;
   std::array<std::complex<double>, 3> coefficients;

   bool check_indices_valid(int, int, int, int) const;
   std::array<int, 4> to_positions_in_ref_order(int, int, int, int) const;
   std::complex<double> get_coefficient(int, int) const;
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

} // namespace flexiblesusy

#endif
