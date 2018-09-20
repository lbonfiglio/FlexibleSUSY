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

#include "generic_vertices.hpp"

#include <stdexcept>

namespace flexiblesusy {

namespace cxx_qft {

std::complex<double> SSVVertex::value(int mi, int si) const
{
   if( mi == minuendIndex && si == subtrahendIndex ) {
      return val;
   }

   if( mi == subtrahendIndex && si == minuendIndex ) {
      return -val;
   }

   throw std::invalid_argument(
      "SSVVertex: Wrong index combination" );
}

VVVVVertex::VVVVVertex(int i, int j, int k, int l,
                       std::complex<double> c1, std::complex<double> c2,
                       std::complex<double> c3)
   : indices({i, j, k, l})
   , coefficients({c1, c2, c3})
{
}

std::complex<double> VVVVVertex::value_1(int i, int j, int k, int l) const
{
   if (!check_indices_valid(i, j, k, l)) {
      throw std::invalid_argument("VVVVVertex: invalid indices");
   }

   const auto permuted_indices = to_positions_in_ref_order(i, j, k, l);

   return get_coefficient(permuted_indices[0], permuted_indices[1]);
}

std::complex<double> VVVVVertex::value_2(int i, int j, int k, int l) const
{
   if (!check_indices_valid(i, j, k, l)) {
      throw std::invalid_argument("VVVVVertex: invalid indices");
   }

   const auto permuted_indices = to_positions_in_ref_order(i, j, k, l);

   return get_coefficient(permuted_indices[0], permuted_indices[2]);
}

std::complex<double> VVVVVertex::value_3(int i, int j, int k, int l) const
{
   if (!check_indices_valid(i, j, k, l)) {
      throw std::invalid_argument("VVVVVertex: invalid indices");
   }

   const auto permuted_indices = to_positions_in_ref_order(i, j, k, l);

   return get_coefficient(permuted_indices[0], permuted_indices[3]);
}

bool VVVVVertex::check_indices_valid(int i, int j, int k, int l) const
{
   bool has_vec_0_id = false;
   bool has_vec_1_id = false;
   bool has_vec_2_id = false;
   bool has_vec_3_id = false;

   if (i == indices[0]) {
      has_vec_0_id = true;
   } else if (i == indices[1]) {
      has_vec_1_id = true;
   } else if (i == indices[2]) {
      has_vec_2_id = true;
   } else if (i == indices[3]) {
      has_vec_3_id = true;
   } else {
      return false;
   }

   const auto is_valid_idx =
      [this, &has_vec_0_id, &has_vec_1_id, &has_vec_2_id, &has_vec_3_id]
      (int idx) {
      if (idx == this->indices[0] && !has_vec_0_id) {
         has_vec_0_id = true;
      } else if (idx == this->indices[1] && !has_vec_1_id) {
         has_vec_1_id = true;
      } else if (idx == this->indices[2] && !has_vec_2_id) {
         has_vec_2_id = true;
      } else if (idx == this->indices[3] && !has_vec_3_id) {
         has_vec_3_id = true;
      } else {
         return false;
      }
      return true;
   };

   if (!is_valid_idx(j)) {
      return false;
   }

   if (!is_valid_idx(k)) {
      return false;
   }

   if (!is_valid_idx(l)) {
      return false;
   }

   return true;
}

std::array<int, 4> VVVVVertex::to_positions_in_ref_order(
   int i, int j, int k, int l) const
{
   bool has_vec_0_id = false;
   bool has_vec_1_id = false;
   bool has_vec_2_id = false;
   bool has_vec_3_id = false;

   std::array<int, 4> result{{-1, -1, -1, -1}};
   const auto set_ref_pos =
      [this, &result, &has_vec_0_id, &has_vec_1_id,
       &has_vec_2_id, &has_vec_3_id] (int pos, int val) {
      if (this->indices[0] == val && !has_vec_0_id) {
         result[pos] = 0;
         has_vec_0_id = true;
      } else if (indices[1] == val && !has_vec_1_id) {
         result[pos] = 1;
         has_vec_1_id = true;
      } else if (indices[2] == val && !has_vec_2_id) {
         result[pos] = 2;
         has_vec_2_id = true;
      } else if (indices[3] == val && !has_vec_3_id) {
         result[pos] = 3;
         has_vec_3_id = true;
      }
   };

   set_ref_pos(0, i);
   set_ref_pos(1, j);
   set_ref_pos(2, k);
   set_ref_pos(3, l);

   return result;
}

std::complex<double> VVVVVertex::get_coefficient(int i, int j) const
{
   const auto key = (i << 2) | j;
   switch (key) {
      /* case 0: i == j is invalid */
   case 1: return coefficients[0];
   case 2: return coefficients[1];
   case 3: return coefficients[2];
   case 4: return coefficients[0];
      /* case 5: i == j is invalid */
   case 6: return coefficients[2];
   case 7: return coefficients[1];
   case 8: return coefficients[1];
   case 9: return coefficients[2];
      /* case 10: i == j is invalid */
   case 11: return coefficients[0];
   case 12: return coefficients[2];
   case 13: return coefficients[1];
   case 14: return coefficients[0];
      /* case 15: i == j is invalid */
   default:
      throw std::invalid_argument(
         "VVVVVertex: invalid index combination");
   }
}

} // namespace cxx_qft

} // namespace flexiblesusy
