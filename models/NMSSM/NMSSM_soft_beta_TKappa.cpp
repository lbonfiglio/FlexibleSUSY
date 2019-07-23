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

// File generated at Tue 23 Jul 2019 15:15:32

#include "NMSSM_soft_parameters.hpp"
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
 * Calculates the 1-loop beta function of TKappa.
 *
 * @return 1-loop beta function
 */
double NMSSM_soft_parameters::calc_beta_TKappa_1_loop(const Soft_traces& soft_traces) const
{


   double beta_TKappa;

   beta_TKappa = Re(6*oneOver16PiSqr*(3*AbsSqr(Kappa)*TKappa + AbsSqr(Lambdax)*
      TKappa + 2*Conj(Lambdax)*Kappa*TLambdax));


   return beta_TKappa;
}

/**
 * Calculates the 2-loop beta function of TKappa.
 *
 * @return 2-loop beta function
 */
double NMSSM_soft_parameters::calc_beta_TKappa_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;


   double beta_TKappa;

   beta_TKappa = Re(-1.2*twoLoop*(30*traceAdjYdTYd*AbsSqr(Lambdax)*Kappa + 10*
      traceAdjYeTYe*AbsSqr(Lambdax)*Kappa + 30*traceAdjYuTYu*AbsSqr(Lambdax)*
      Kappa + 6*MassB*AbsSqr(Lambdax)*Kappa*Sqr(g1) + 30*MassWB*AbsSqr(Lambdax)
      *Kappa*Sqr(g2) + 15*traceYdAdjYd*AbsSqr(Lambdax)*TKappa + 5*traceYeAdjYe*
      AbsSqr(Lambdax)*TKappa + 15*traceYuAdjYu*AbsSqr(Lambdax)*TKappa + 60*
      AbsSqr(Kappa)*AbsSqr(Lambdax)*TKappa - 3*AbsSqr(Lambdax)*Sqr(g1)*TKappa -
      15*AbsSqr(Lambdax)*Sqr(g2)*TKappa + 100*Sqr(Conj(Kappa))*Sqr(Kappa)*
      TKappa + 10*Sqr(Conj(Lambdax))*Sqr(Lambdax)*TKappa + 30*traceYdAdjYd*Conj
      (Lambdax)*Kappa*TLambdax + 10*traceYeAdjYe*Conj(Lambdax)*Kappa*TLambdax +
      30*traceYuAdjYu*Conj(Lambdax)*Kappa*TLambdax - 6*Conj(Lambdax)*Kappa*Sqr(
      g1)*TLambdax - 30*Conj(Lambdax)*Kappa*Sqr(g2)*TLambdax + 40*Kappa*Lambdax
      *Sqr(Conj(Lambdax))*TLambdax + 40*Conj(Kappa)*Conj(Lambdax)*Sqr(Kappa)*
      TLambdax));


   return beta_TKappa;
}

/**
 * Calculates the 3-loop beta function of TKappa.
 *
 * @return 3-loop beta function
 */
double NMSSM_soft_parameters::calc_beta_TKappa_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_TKappa;

   beta_TKappa = 0;


   return beta_TKappa;
}

/**
 * Calculates the 4-loop beta function of TKappa.
 *
 * @return 4-loop beta function
 */
double NMSSM_soft_parameters::calc_beta_TKappa_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_TKappa;

   beta_TKappa = 0;


   return beta_TKappa;
}

/**
 * Calculates the 5-loop beta function of TKappa.
 *
 * @return 5-loop beta function
 */
double NMSSM_soft_parameters::calc_beta_TKappa_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_TKappa;

   beta_TKappa = 0;


   return beta_TKappa;
}

} // namespace flexiblesusy
