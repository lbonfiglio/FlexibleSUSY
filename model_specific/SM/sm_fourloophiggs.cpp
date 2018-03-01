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

#include "sm_fourloophiggs.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {
namespace sm_fourloophiggs {

/**
 * Standard Model Higgs mass 4-loop contribution, \f$O(\alpha_t
 * \alpha_s^3)\f$.  Taken from arxiv:1508.00912, Eq. (5.5).
 *
 * @note The result contains the 4-loop tadpole diagrams.  It is
 * therefore not 1-particle irreducible (1PI).
 *
 * @param scale renormalization scale
 * @param mt MS-bar top mass
 * @param yt MS-bar Yukawa coupling
 * @param g3 MS-bar strong gauge coupling
 *
 * @return real part of 4-loop correction \f$O(\alpha_t \alpha_s^3)\f$
 */
double delta_mh_4loop_at_as_as_as_sm(
   double scale, double mt, double yt, double g3)
{
   const double yt2 = Sqr(yt);
   const double mt2 = Sqr(mt);
   const double g36 = Power6(g3);
   const double Q2 = Sqr(scale);
   const double LogT = FiniteLog(mt2 / Q2);
   const double LogT2 = Sqr(LogT);
   const double LogT3 = Power3(LogT);
   const double LogT4 = Power4(LogT);

   const double result =
      g36*yt2*mt2 * (23925.974863118638
                     + 6435.327201984095*LogT
                     - 20675.746270332107*LogT2
                     - 3456*LogT3
                     + 5520*LogT4);

   return result * fourLoop;
}

} // namespace sm_fourloophiggs
} // namespace flexiblesusy
