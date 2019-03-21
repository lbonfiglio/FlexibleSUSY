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

// File generated at Thu 21 Mar 2019 14:09:02

#include "CNMSSMEFTHiggs_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Lambdax.
 *
 * @return 1-loop beta function
 */
double CNMSSMEFTHiggs_susy_parameters::calc_beta_Lambdax_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_Lambdax;

   beta_Lambdax = Re(0.2*oneOver16PiSqr*Lambdax*(15*traceYdAdjYd + 5*
      traceYeAdjYe + 15*traceYuAdjYu + 10*AbsSqr(Kappa) + 20*AbsSqr(Lambdax) -
      3*Sqr(g1) - 15*Sqr(g2)));


   return beta_Lambdax;
}

/**
 * Calculates the 2-loop beta function of Lambdax.
 *
 * @return 2-loop beta function
 */
double CNMSSMEFTHiggs_susy_parameters::calc_beta_Lambdax_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_Lambdax;

   beta_Lambdax = Re(-0.02*twoLoop*Lambdax*(450*traceYdAdjYdYdAdjYd + 300*
      traceYdAdjYuYuAdjYd + 150*traceYeAdjYeYeAdjYe + 450*traceYuAdjYuYuAdjYu +
      450*traceYdAdjYd*AbsSqr(Lambdax) + 150*traceYeAdjYe*AbsSqr(Lambdax) + 450
      *traceYuAdjYu*AbsSqr(Lambdax) + 600*AbsSqr(Kappa)*AbsSqr(Lambdax) - 207*
      Quad(g1) - 375*Quad(g2) + 20*traceYdAdjYd*Sqr(g1) - 60*traceYeAdjYe*Sqr(
      g1) - 40*traceYuAdjYu*Sqr(g1) - 60*AbsSqr(Lambdax)*Sqr(g1) - 300*AbsSqr(
      Lambdax)*Sqr(g2) - 90*Sqr(g1)*Sqr(g2) - 800*traceYdAdjYd*Sqr(g3) - 800*
      traceYuAdjYu*Sqr(g3) + 400*Sqr(Conj(Kappa))*Sqr(Kappa) + 500*Sqr(Conj(
      Lambdax))*Sqr(Lambdax)));


   return beta_Lambdax;
}

/**
 * Calculates the 3-loop beta function of Lambdax.
 *
 * @return 3-loop beta function
 */
double CNMSSMEFTHiggs_susy_parameters::calc_beta_Lambdax_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambdax;

   beta_Lambdax = 0;


   return beta_Lambdax;
}

/**
 * Calculates the 4-loop beta function of Lambdax.
 *
 * @return 4-loop beta function
 */
double CNMSSMEFTHiggs_susy_parameters::calc_beta_Lambdax_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambdax;

   beta_Lambdax = 0;


   return beta_Lambdax;
}

/**
 * Calculates the 5-loop beta function of Lambdax.
 *
 * @return 5-loop beta function
 */
double CNMSSMEFTHiggs_susy_parameters::calc_beta_Lambdax_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambdax;

   beta_Lambdax = 0;


   return beta_Lambdax;
}

} // namespace flexiblesusy
