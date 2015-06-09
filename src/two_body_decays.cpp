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

#include "two_body_decays.hpp"
#include "wrappers.hpp"
#include "numerics2.hpp"

/**
 * @file two_body_decays.cpp
 * @brief contains the implementation of the two-body decay functions
 */

namespace flexiblesusy {

namespace two_body_decays {

double kallen(double x, double y, double z)
{
   return x * x + y * y + z * z - 2 * x * y - 2 * x * z - 2 * y * z;
}

double scalar_to_fermion_fermion(double ms, double mf1, double mf2,
                                 std::complex<double> cpl_left,
                                 std::complex<double> cpl_right)
{
   // return 0 if not kinematically allowed
   if (Abs(ms) < Abs(mf1) + Abs(mf2))
      return 0.;

   const double ms2 = Sqr(ms);
   const double mf12 = Sqr(mf1);
   const double mf22 = Sqr(mf2);

   const double kappa = Sqrt(kallen(ms2, mf12, mf22)); 

   // note - no symmetry factor here at the moment
   const double width = Pi * oneOver16PiSqr * kappa *
      ((AbsSqr(cpl_left) + AbsSqr(cpl_right)) * (ms2 - mf12 - mf22)
       - 4.0 * Re(Conj(cpl_left) * cpl_right) * mf1 * mf2) /
      Power(Abs(ms), 3);

   return width;
}

} // namespace two_body_decays

} // namespace flexiblesusy
