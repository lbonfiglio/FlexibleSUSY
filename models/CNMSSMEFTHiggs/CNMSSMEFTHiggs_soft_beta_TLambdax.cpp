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

// File generated at Thu 21 Mar 2019 14:09:19

#include "CNMSSMEFTHiggs_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

namespace {

template <typename Derived>
typename Eigen::MatrixBase<Derived>::PlainObject operator+(const Eigen::MatrixBase<Derived>& m, double n)
{
   return m + Eigen::Matrix<double,
                            Eigen::MatrixBase<Derived>::RowsAtCompileTime,
                            Eigen::MatrixBase<Derived>::ColsAtCompileTime>::Identity() * n;
}

template <typename Derived>
typename Eigen::MatrixBase<Derived>::PlainObject operator+(double n, const Eigen::MatrixBase<Derived>& m)
{
   return m + Eigen::Matrix<double,
                            Eigen::MatrixBase<Derived>::RowsAtCompileTime,
                            Eigen::MatrixBase<Derived>::ColsAtCompileTime>::Identity() * n;
}

template <typename Derived>
typename Eigen::MatrixBase<Derived>::PlainObject operator-(const Eigen::MatrixBase<Derived>& m, double n)
{
   return m - Eigen::Matrix<double,
                            Eigen::MatrixBase<Derived>::RowsAtCompileTime,
                            Eigen::MatrixBase<Derived>::ColsAtCompileTime>::Identity() * n;
}

template <typename Derived>
typename Eigen::MatrixBase<Derived>::PlainObject operator-(double n, const Eigen::MatrixBase<Derived>& m)
{
   return - m + Eigen::Matrix<double,
                            Eigen::MatrixBase<Derived>::RowsAtCompileTime,
                            Eigen::MatrixBase<Derived>::ColsAtCompileTime>::Identity() * n;
}

} // anonymous namespace

/**
 * Calculates the 1-loop beta function of TLambdax.
 *
 * @return 1-loop beta function
 */
double CNMSSMEFTHiggs_soft_parameters::calc_beta_TLambdax_1_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;


   double beta_TLambdax;

   beta_TLambdax = Re(0.2*oneOver16PiSqr*(30*traceAdjYdTYd*Lambdax + 10*
      traceAdjYeTYe*Lambdax + 30*traceAdjYuTYu*Lambdax + 6*MassB*Lambdax*Sqr(g1
      ) + 30*MassWB*Lambdax*Sqr(g2) + 20*Conj(Kappa)*Lambdax*TKappa + 15*
      traceYdAdjYd*TLambdax + 5*traceYeAdjYe*TLambdax + 15*traceYuAdjYu*
      TLambdax + 10*AbsSqr(Kappa)*TLambdax + 60*AbsSqr(Lambdax)*TLambdax - 3*
      Sqr(g1)*TLambdax - 15*Sqr(g2)*TLambdax));


   return beta_TLambdax;
}

/**
 * Calculates the 2-loop beta function of TLambdax.
 *
 * @return 2-loop beta function
 */
double CNMSSMEFTHiggs_soft_parameters::calc_beta_TLambdax_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYdTYdAdjYd = TRACE_STRUCT.traceYdAdjYdTYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYdAdjYuTYuAdjYd = TRACE_STRUCT.traceYdAdjYuTYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYeAdjYeTYeAdjYe = TRACE_STRUCT.traceYeAdjYeTYeAdjYe;
   const double traceYuAdjYdTYdAdjYu = TRACE_STRUCT.traceYuAdjYdTYdAdjYu;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;
   const double traceYuAdjYuTYuAdjYu = TRACE_STRUCT.traceYuAdjYuTYuAdjYu;


   double beta_TLambdax;

   beta_TLambdax = Re(0.02*twoLoop*(-1800*traceYdAdjYdTYdAdjYd*Lambdax - 600*
      traceYdAdjYuTYuAdjYd*Lambdax - 600*traceYeAdjYeTYeAdjYe*Lambdax - 600*
      traceYuAdjYdTYdAdjYu*Lambdax - 1800*traceYuAdjYuTYuAdjYu*Lambdax - 828*
      MassB*Lambdax*Quad(g1) - 1500*MassWB*Lambdax*Quad(g2) - 40*traceAdjYdTYd*
      Lambdax*Sqr(g1) + 120*traceAdjYeTYe*Lambdax*Sqr(g1) + 80*traceAdjYuTYu*
      Lambdax*Sqr(g1) + 40*MassB*traceYdAdjYd*Lambdax*Sqr(g1) - 120*MassB*
      traceYeAdjYe*Lambdax*Sqr(g1) - 80*MassB*traceYuAdjYu*Lambdax*Sqr(g1) -
      180*MassB*Lambdax*Sqr(g1)*Sqr(g2) - 180*MassWB*Lambdax*Sqr(g1)*Sqr(g2) +
      1600*traceAdjYdTYd*Lambdax*Sqr(g3) + 1600*traceAdjYuTYu*Lambdax*Sqr(g3) -
      1600*MassG*traceYdAdjYd*Lambdax*Sqr(g3) - 1600*MassG*traceYuAdjYu*Lambdax
      *Sqr(g3) - 900*traceAdjYdTYd*Conj(Lambdax)*Sqr(Lambdax) - 300*
      traceAdjYeTYe*Conj(Lambdax)*Sqr(Lambdax) - 900*traceAdjYuTYu*Conj(Lambdax
      )*Sqr(Lambdax) - 120*MassB*Conj(Lambdax)*Sqr(g1)*Sqr(Lambdax) - 600*
      MassWB*Conj(Lambdax)*Sqr(g2)*Sqr(Lambdax) - 1600*Kappa*Lambdax*Sqr(Conj(
      Kappa))*TKappa - 1200*Conj(Kappa)*Conj(Lambdax)*Sqr(Lambdax)*TKappa - 450
      *traceYdAdjYdYdAdjYd*TLambdax - 300*traceYdAdjYuYuAdjYd*TLambdax - 150*
      traceYeAdjYeYeAdjYe*TLambdax - 450*traceYuAdjYuYuAdjYu*TLambdax - 1350*
      traceYdAdjYd*AbsSqr(Lambdax)*TLambdax - 450*traceYeAdjYe*AbsSqr(Lambdax)*
      TLambdax - 1350*traceYuAdjYu*AbsSqr(Lambdax)*TLambdax - 1800*AbsSqr(Kappa
      )*AbsSqr(Lambdax)*TLambdax + 207*Quad(g1)*TLambdax + 375*Quad(g2)*
      TLambdax - 20*traceYdAdjYd*Sqr(g1)*TLambdax + 60*traceYeAdjYe*Sqr(g1)*
      TLambdax + 40*traceYuAdjYu*Sqr(g1)*TLambdax + 180*AbsSqr(Lambdax)*Sqr(g1)
      *TLambdax + 900*AbsSqr(Lambdax)*Sqr(g2)*TLambdax + 90*Sqr(g1)*Sqr(g2)*
      TLambdax + 800*traceYdAdjYd*Sqr(g3)*TLambdax + 800*traceYuAdjYu*Sqr(g3)*
      TLambdax - 400*Sqr(Conj(Kappa))*Sqr(Kappa)*TLambdax - 2500*Sqr(Conj(
      Lambdax))*Sqr(Lambdax)*TLambdax));


   return beta_TLambdax;
}

/**
 * Calculates the 3-loop beta function of TLambdax.
 *
 * @return 3-loop beta function
 */
double CNMSSMEFTHiggs_soft_parameters::calc_beta_TLambdax_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_TLambdax;

   beta_TLambdax = 0;


   return beta_TLambdax;
}

/**
 * Calculates the 4-loop beta function of TLambdax.
 *
 * @return 4-loop beta function
 */
double CNMSSMEFTHiggs_soft_parameters::calc_beta_TLambdax_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_TLambdax;

   beta_TLambdax = 0;


   return beta_TLambdax;
}

/**
 * Calculates the 5-loop beta function of TLambdax.
 *
 * @return 5-loop beta function
 */
double CNMSSMEFTHiggs_soft_parameters::calc_beta_TLambdax_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_TLambdax;

   beta_TLambdax = 0;


   return beta_TLambdax;
}

} // namespace flexiblesusy
