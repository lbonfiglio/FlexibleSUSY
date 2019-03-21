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

// File generated at Thu 21 Mar 2019 19:35:32

#include "CNMSSM_soft_parameters.hpp"
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
 * Calculates the 1-loop beta function of mHu2.
 *
 * @return 1-loop beta function
 */
double CNMSSM_soft_parameters::calc_beta_mHu2_1_loop(const Soft_traces& soft_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceconjTYuTpTYu = TRACE_STRUCT.traceconjTYuTpTYu;
   const double tracemq2AdjYuYu = TRACE_STRUCT.tracemq2AdjYuYu;
   const double tracemu2YuAdjYu = TRACE_STRUCT.tracemu2YuAdjYu;
   const double Tr11 = TRACE_STRUCT.Tr11;


   double beta_mHu2;

   beta_mHu2 = Re(0.2*oneOver16PiSqr*(3.872983346207417*g1*Tr11 + 30*
      traceconjTYuTpTYu + 30*tracemq2AdjYuYu + 30*tracemu2YuAdjYu + 30*mHu2*
      traceYuAdjYu + 10*mHd2*AbsSqr(Lambdax) + 10*mHu2*AbsSqr(Lambdax) + 10*ms2
      *AbsSqr(Lambdax) + 10*AbsSqr(TLambdax) - 6*AbsSqr(MassB)*Sqr(g1) - 30*
      AbsSqr(MassWB)*Sqr(g2)));


   return beta_mHu2;
}

/**
 * Calculates the 2-loop beta function of mHu2.
 *
 * @return 2-loop beta function
 */
double CNMSSM_soft_parameters::calc_beta_mHu2_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceconjTYdTpYd = TRACE_STRUCT.traceconjTYdTpYd;
   const double traceconjTYdTpTYd = TRACE_STRUCT.traceconjTYdTpTYd;
   const double traceconjTYeTpYe = TRACE_STRUCT.traceconjTYeTpYe;
   const double traceconjTYeTpTYe = TRACE_STRUCT.traceconjTYeTpTYe;
   const double traceconjTYuTpYu = TRACE_STRUCT.traceconjTYuTpYu;
   const double traceconjTYuTpTYu = TRACE_STRUCT.traceconjTYuTpTYu;
   const double tracemd2YdAdjYd = TRACE_STRUCT.tracemd2YdAdjYd;
   const double traceme2YeAdjYe = TRACE_STRUCT.traceme2YeAdjYe;
   const double traceml2AdjYeYe = TRACE_STRUCT.traceml2AdjYeYe;
   const double tracemq2AdjYdYd = TRACE_STRUCT.tracemq2AdjYdYd;
   const double tracemq2AdjYuYu = TRACE_STRUCT.tracemq2AdjYuYu;
   const double tracemu2YuAdjYu = TRACE_STRUCT.tracemu2YuAdjYu;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYdAdjYuTYuAdjTYd = TRACE_STRUCT.traceYdAdjYuTYuAdjTYd;
   const double traceYdAdjTYuTYuAdjYd = TRACE_STRUCT.traceYdAdjTYuTYuAdjYd;
   const double traceYuAdjYdTYdAdjTYu = TRACE_STRUCT.traceYuAdjYdTYdAdjTYu;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;
   const double traceYuAdjYuTYuAdjTYu = TRACE_STRUCT.traceYuAdjYuTYuAdjTYu;
   const double traceYuAdjTYdTYdAdjYu = TRACE_STRUCT.traceYuAdjTYdTYdAdjYu;
   const double traceYuAdjTYuTYuAdjYu = TRACE_STRUCT.traceYuAdjTYuTYuAdjYu;
   const double tracemd2YdAdjYuYuAdjYd = TRACE_STRUCT.tracemd2YdAdjYuYuAdjYd;
   const double tracemq2AdjYdYdAdjYuYu = TRACE_STRUCT.tracemq2AdjYdYdAdjYuYu;
   const double tracemq2AdjYuYuAdjYdYd = TRACE_STRUCT.tracemq2AdjYuYuAdjYdYd;
   const double tracemq2AdjYuYuAdjYuYu = TRACE_STRUCT.tracemq2AdjYuYuAdjYuYu;
   const double tracemu2YuAdjYdYdAdjYu = TRACE_STRUCT.tracemu2YuAdjYdYdAdjYu;
   const double tracemu2YuAdjYuYuAdjYu = TRACE_STRUCT.tracemu2YuAdjYuYuAdjYu;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr22 = TRACE_STRUCT.Tr22;


   double beta_mHu2;

   beta_mHu2 = Re(0.04*twoLoop*(77.45966692414834*g1*Tr31 - 150*
      tracemd2YdAdjYuYuAdjYd - 150*tracemq2AdjYdYdAdjYuYu - 150*
      tracemq2AdjYuYuAdjYdYd - 900*tracemq2AdjYuYuAdjYuYu - 150*
      tracemu2YuAdjYdYdAdjYu - 900*tracemu2YuAdjYuYuAdjYu - 150*
      traceYdAdjTYuTYuAdjYd - 150*traceYdAdjYuTYuAdjTYd - 150*mHd2*
      traceYdAdjYuYuAdjYd - 150*mHu2*traceYdAdjYuYuAdjYd - 150*
      traceYuAdjTYdTYdAdjYu - 900*traceYuAdjTYuTYuAdjYu - 150*
      traceYuAdjYdTYdAdjTYu - 900*traceYuAdjYuTYuAdjTYu - 900*mHu2*
      traceYuAdjYuYuAdjYu - 150*traceconjTYdTpTYd*AbsSqr(Lambdax) - 50*
      traceconjTYeTpTYe*AbsSqr(Lambdax) - 150*tracemd2YdAdjYd*AbsSqr(Lambdax) -
      50*traceme2YeAdjYe*AbsSqr(Lambdax) - 50*traceml2AdjYeYe*AbsSqr(Lambdax) -
      150*tracemq2AdjYdYd*AbsSqr(Lambdax) - 300*mHd2*traceYdAdjYd*AbsSqr(
      Lambdax) - 150*mHu2*traceYdAdjYd*AbsSqr(Lambdax) - 150*ms2*traceYdAdjYd*
      AbsSqr(Lambdax) - 100*mHd2*traceYeAdjYe*AbsSqr(Lambdax) - 50*mHu2*
      traceYeAdjYe*AbsSqr(Lambdax) - 50*ms2*traceYeAdjYe*AbsSqr(Lambdax) - 100*
      mHd2*AbsSqr(Kappa)*AbsSqr(Lambdax) - 100*mHu2*AbsSqr(Kappa)*AbsSqr(
      Lambdax) - 400*ms2*AbsSqr(Kappa)*AbsSqr(Lambdax) - 100*AbsSqr(Lambdax)*
      AbsSqr(TKappa) - 150*traceYdAdjYd*AbsSqr(TLambdax) - 50*traceYeAdjYe*
      AbsSqr(TLambdax) - 100*AbsSqr(Kappa)*AbsSqr(TLambdax) - 600*AbsSqr(
      Lambdax)*AbsSqr(TLambdax) - 150*traceAdjYdTYd*Conj(TLambdax)*Lambdax - 50
      *traceAdjYeTYe*Conj(TLambdax)*Lambdax + 621*AbsSqr(MassB)*Quad(g1) + 150*
      Tr22*Quad(g2) + 825*AbsSqr(MassWB)*Quad(g2) + 30*Tr2U111*Sqr(g1) + 40*
      traceconjTYuTpTYu*Sqr(g1) - 40*MassB*traceconjTYuTpYu*Sqr(g1) + 40*
      tracemq2AdjYuYu*Sqr(g1) + 40*tracemu2YuAdjYu*Sqr(g1) + 40*mHu2*
      traceYuAdjYu*Sqr(g1) + 80*traceYuAdjYu*AbsSqr(MassB)*Sqr(g1) - 40*
      traceAdjYuTYu*Conj(MassB)*Sqr(g1) + 90*AbsSqr(MassB)*Sqr(g1)*Sqr(g2) + 90
      *AbsSqr(MassWB)*Sqr(g1)*Sqr(g2) + 45*MassWB*Conj(MassB)*Sqr(g1)*Sqr(g2) +
      45*MassB*Conj(MassWB)*Sqr(g1)*Sqr(g2) + 800*traceconjTYuTpTYu*Sqr(g3) -
      800*MassG*traceconjTYuTpYu*Sqr(g3) + 800*tracemq2AdjYuYu*Sqr(g3) + 800*
      tracemu2YuAdjYu*Sqr(g3) + 800*mHu2*traceYuAdjYu*Sqr(g3) + 1600*
      traceYuAdjYu*AbsSqr(MassG)*Sqr(g3) - 800*traceAdjYuTYu*Conj(MassG)*Sqr(g3
      ) - 300*mHd2*Sqr(Conj(Lambdax))*Sqr(Lambdax) - 300*mHu2*Sqr(Conj(Lambdax)
      )*Sqr(Lambdax) - 300*ms2*Sqr(Conj(Lambdax))*Sqr(Lambdax) - 100*Conj(Kappa
      )*Conj(TLambdax)*Lambdax*TKappa - 150*traceconjTYdTpYd*Conj(Lambdax)*
      TLambdax - 50*traceconjTYeTpYe*Conj(Lambdax)*TLambdax - 100*Conj(Lambdax)
      *Conj(TKappa)*Kappa*TLambdax));


   return beta_mHu2;
}

/**
 * Calculates the 3-loop beta function of mHu2.
 *
 * @return 3-loop beta function
 */
double CNMSSM_soft_parameters::calc_beta_mHu2_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mHu2;

   beta_mHu2 = 0;


   return beta_mHu2;
}

/**
 * Calculates the 4-loop beta function of mHu2.
 *
 * @return 4-loop beta function
 */
double CNMSSM_soft_parameters::calc_beta_mHu2_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mHu2;

   beta_mHu2 = 0;


   return beta_mHu2;
}

/**
 * Calculates the 5-loop beta function of mHu2.
 *
 * @return 5-loop beta function
 */
double CNMSSM_soft_parameters::calc_beta_mHu2_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mHu2;

   beta_mHu2 = 0;


   return beta_mHu2;
}

} // namespace flexiblesusy
