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

// File generated at Tue 16 Jul 2019 14:55:54

#include "MSSM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of g2.
 *
 * @return 1-loop beta function
 */
double MSSM_susy_parameters::calc_beta_g2_1_loop(const Susy_traces& susy_traces) const
{


   double beta_g2;

   beta_g2 = Re(oneOver16PiSqr*Cube(g2));


   return beta_g2;
}

/**
 * Calculates the 2-loop beta function of g2.
 *
 * @return 2-loop beta function
 */
double MSSM_susy_parameters::calc_beta_g2_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_g2;

   beta_g2 = Re(0.2*twoLoop*Cube(g2)*(-30*traceYdAdjYd - 10*traceYeAdjYe - 30*
      traceYuAdjYu + 9*Sqr(g1) + 125*Sqr(g2) + 120*Sqr(g3)));


   return beta_g2;
}

/**
 * Calculates the 3-loop beta function of g2.
 *
 * @return 3-loop beta function
 */
double MSSM_susy_parameters::calc_beta_g2_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)

   const double traceAdjYdYd = TRACE_STRUCT.traceAdjYdYd;
   const double traceAdjYeYe = TRACE_STRUCT.traceAdjYeYe;
   const double traceAdjYuYu = TRACE_STRUCT.traceAdjYuYu;
   const double traceAdjYdYdAdjYdYd = TRACE_STRUCT.traceAdjYdYdAdjYdYd;
   const double traceAdjYeYeAdjYeYe = TRACE_STRUCT.traceAdjYeYeAdjYeYe;
   const double traceAdjYuYuAdjYdYd = TRACE_STRUCT.traceAdjYuYuAdjYdYd;
   const double traceAdjYuYuAdjYuYu = TRACE_STRUCT.traceAdjYuYuAdjYuYu;


   double beta_g2;

   beta_g2 = Re(0.04*threeLoop*Cube(g2)*(600*traceAdjYdYdAdjYdYd + 300*
      traceAdjYdYd*traceAdjYeYe + 200*traceAdjYeYeAdjYeYe + 300*
      traceAdjYuYuAdjYdYd + 600*traceAdjYuYuAdjYuYu - 457*Quad(g1) + 875*Quad(
      g2) + 1100*Quad(g3) - 55*traceAdjYdYd*Sqr(g1) - 105*traceAdjYeYe*Sqr(g1)
      - 145*traceAdjYuYu*Sqr(g1) - 825*traceAdjYdYd*Sqr(g2) - 275*traceAdjYeYe*
      Sqr(g2) - 825*traceAdjYuYu*Sqr(g2) + 45*Sqr(g1)*Sqr(g2) - 800*
      traceAdjYdYd*Sqr(g3) - 800*traceAdjYuYu*Sqr(g3) - 40*Sqr(g1)*Sqr(g3) +
      600*Sqr(g2)*Sqr(g3) + 450*Sqr(traceAdjYdYd) + 50*Sqr(traceAdjYeYe) + 450*
      Sqr(traceAdjYuYu)));


   return beta_g2;
}

/**
 * Calculates the 4-loop beta function of g2.
 *
 * @return 4-loop beta function
 */
double MSSM_susy_parameters::calc_beta_g2_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_g2;

   beta_g2 = 0;


   return beta_g2;
}

/**
 * Calculates the 5-loop beta function of g2.
 *
 * @return 5-loop beta function
 */
double MSSM_susy_parameters::calc_beta_g2_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_g2;

   beta_g2 = 0;


   return beta_g2;
}

} // namespace flexiblesusy
