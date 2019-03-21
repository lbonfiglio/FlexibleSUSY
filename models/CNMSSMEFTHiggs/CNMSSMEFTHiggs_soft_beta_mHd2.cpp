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

// File generated at Thu 21 Mar 2019 14:09:23

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
 * Calculates the 1-loop beta function of mHd2.
 *
 * @return 1-loop beta function
 */
double CNMSSMEFTHiggs_soft_parameters::calc_beta_mHd2_1_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceconjTYdTpTYd = TRACE_STRUCT.traceconjTYdTpTYd;
   const double traceconjTYeTpTYe = TRACE_STRUCT.traceconjTYeTpTYe;
   const double tracemd2YdAdjYd = TRACE_STRUCT.tracemd2YdAdjYd;
   const double traceme2YeAdjYe = TRACE_STRUCT.traceme2YeAdjYe;
   const double traceml2AdjYeYe = TRACE_STRUCT.traceml2AdjYeYe;
   const double tracemq2AdjYdYd = TRACE_STRUCT.tracemq2AdjYdYd;
   const double Tr11 = TRACE_STRUCT.Tr11;


   double beta_mHd2;

   beta_mHd2 = Re(0.2*oneOver16PiSqr*(-3.872983346207417*g1*Tr11 + 30*
      traceconjTYdTpTYd + 10*traceconjTYeTpTYe + 30*tracemd2YdAdjYd + 10*
      traceme2YeAdjYe + 10*traceml2AdjYeYe + 30*tracemq2AdjYdYd + 30*mHd2*
      traceYdAdjYd + 10*mHd2*traceYeAdjYe + 10*mHd2*AbsSqr(Lambdax) + 10*mHu2*
      AbsSqr(Lambdax) + 10*ms2*AbsSqr(Lambdax) + 10*AbsSqr(TLambdax) - 6*AbsSqr
      (MassB)*Sqr(g1) - 30*AbsSqr(MassWB)*Sqr(g2)));


   return beta_mHd2;
}

/**
 * Calculates the 2-loop beta function of mHd2.
 *
 * @return 2-loop beta function
 */
double CNMSSMEFTHiggs_soft_parameters::calc_beta_mHd2_2_loop(const Soft_traces& soft_traces) const
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
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYdTYdAdjTYd = TRACE_STRUCT.traceYdAdjYdTYdAdjTYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYdAdjYuTYuAdjTYd = TRACE_STRUCT.traceYdAdjYuTYuAdjTYd;
   const double traceYdAdjTYdTYdAdjYd = TRACE_STRUCT.traceYdAdjTYdTYdAdjYd;
   const double traceYdAdjTYuTYuAdjYd = TRACE_STRUCT.traceYdAdjTYuTYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYeAdjYeTYeAdjTYe = TRACE_STRUCT.traceYeAdjYeTYeAdjTYe;
   const double traceYeAdjTYeTYeAdjYe = TRACE_STRUCT.traceYeAdjTYeTYeAdjYe;
   const double traceYuAdjYdTYdAdjTYu = TRACE_STRUCT.traceYuAdjYdTYdAdjTYu;
   const double traceYuAdjTYdTYdAdjYu = TRACE_STRUCT.traceYuAdjTYdTYdAdjYu;
   const double tracemd2YdAdjYdYdAdjYd = TRACE_STRUCT.tracemd2YdAdjYdYdAdjYd;
   const double tracemd2YdAdjYuYuAdjYd = TRACE_STRUCT.tracemd2YdAdjYuYuAdjYd;
   const double traceme2YeAdjYeYeAdjYe = TRACE_STRUCT.traceme2YeAdjYeYeAdjYe;
   const double traceml2AdjYeYeAdjYeYe = TRACE_STRUCT.traceml2AdjYeYeAdjYeYe;
   const double tracemq2AdjYdYdAdjYdYd = TRACE_STRUCT.tracemq2AdjYdYdAdjYdYd;
   const double tracemq2AdjYdYdAdjYuYu = TRACE_STRUCT.tracemq2AdjYdYdAdjYuYu;
   const double tracemq2AdjYuYuAdjYdYd = TRACE_STRUCT.tracemq2AdjYuYuAdjYdYd;
   const double tracemu2YuAdjYdYdAdjYu = TRACE_STRUCT.tracemu2YuAdjYdYdAdjYu;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr22 = TRACE_STRUCT.Tr22;


   double beta_mHd2;

   beta_mHd2 = Re(0.04*twoLoop*(-77.45966692414834*g1*Tr31 - 900*
      tracemd2YdAdjYdYdAdjYd - 150*tracemd2YdAdjYuYuAdjYd - 300*
      traceme2YeAdjYeYeAdjYe - 300*traceml2AdjYeYeAdjYeYe - 900*
      tracemq2AdjYdYdAdjYdYd - 150*tracemq2AdjYdYdAdjYuYu - 150*
      tracemq2AdjYuYuAdjYdYd - 150*tracemu2YuAdjYdYdAdjYu - 900*
      traceYdAdjTYdTYdAdjYd - 150*traceYdAdjTYuTYuAdjYd - 900*
      traceYdAdjYdTYdAdjTYd - 900*mHd2*traceYdAdjYdYdAdjYd - 150*
      traceYdAdjYuTYuAdjTYd - 150*mHd2*traceYdAdjYuYuAdjYd - 150*mHu2*
      traceYdAdjYuYuAdjYd - 300*traceYeAdjTYeTYeAdjYe - 300*
      traceYeAdjYeTYeAdjTYe - 300*mHd2*traceYeAdjYeYeAdjYe - 150*
      traceYuAdjTYdTYdAdjYu - 150*traceYuAdjYdTYdAdjTYu - 150*traceconjTYuTpTYu
      *AbsSqr(Lambdax) - 150*tracemq2AdjYuYu*AbsSqr(Lambdax) - 150*
      tracemu2YuAdjYu*AbsSqr(Lambdax) - 150*mHd2*traceYuAdjYu*AbsSqr(Lambdax) -
      300*mHu2*traceYuAdjYu*AbsSqr(Lambdax) - 150*ms2*traceYuAdjYu*AbsSqr(
      Lambdax) - 100*mHd2*AbsSqr(Kappa)*AbsSqr(Lambdax) - 100*mHu2*AbsSqr(Kappa
      )*AbsSqr(Lambdax) - 400*ms2*AbsSqr(Kappa)*AbsSqr(Lambdax) - 100*AbsSqr(
      Lambdax)*AbsSqr(TKappa) - 150*traceYuAdjYu*AbsSqr(TLambdax) - 100*AbsSqr(
      Kappa)*AbsSqr(TLambdax) - 600*AbsSqr(Lambdax)*AbsSqr(TLambdax) - 150*
      traceAdjYuTYu*Conj(TLambdax)*Lambdax + 621*AbsSqr(MassB)*Quad(g1) + 150*
      Tr22*Quad(g2) + 825*AbsSqr(MassWB)*Quad(g2) + 30*Tr2U111*Sqr(g1) - 20*
      traceconjTYdTpTYd*Sqr(g1) + 20*MassB*traceconjTYdTpYd*Sqr(g1) + 60*
      traceconjTYeTpTYe*Sqr(g1) - 60*MassB*traceconjTYeTpYe*Sqr(g1) - 20*
      tracemd2YdAdjYd*Sqr(g1) + 60*traceme2YeAdjYe*Sqr(g1) + 60*traceml2AdjYeYe
      *Sqr(g1) - 20*tracemq2AdjYdYd*Sqr(g1) - 20*mHd2*traceYdAdjYd*Sqr(g1) + 60
      *mHd2*traceYeAdjYe*Sqr(g1) - 40*traceYdAdjYd*AbsSqr(MassB)*Sqr(g1) + 120*
      traceYeAdjYe*AbsSqr(MassB)*Sqr(g1) + 20*traceAdjYdTYd*Conj(MassB)*Sqr(g1)
      - 60*traceAdjYeTYe*Conj(MassB)*Sqr(g1) + 90*AbsSqr(MassB)*Sqr(g1)*Sqr(g2)
      + 90*AbsSqr(MassWB)*Sqr(g1)*Sqr(g2) + 45*MassWB*Conj(MassB)*Sqr(g1)*Sqr(
      g2) + 45*MassB*Conj(MassWB)*Sqr(g1)*Sqr(g2) + 800*traceconjTYdTpTYd*Sqr(
      g3) - 800*MassG*traceconjTYdTpYd*Sqr(g3) + 800*tracemd2YdAdjYd*Sqr(g3) +
      800*tracemq2AdjYdYd*Sqr(g3) + 800*mHd2*traceYdAdjYd*Sqr(g3) + 1600*
      traceYdAdjYd*AbsSqr(MassG)*Sqr(g3) - 800*traceAdjYdTYd*Conj(MassG)*Sqr(g3
      ) - 300*mHd2*Sqr(Conj(Lambdax))*Sqr(Lambdax) - 300*mHu2*Sqr(Conj(Lambdax)
      )*Sqr(Lambdax) - 300*ms2*Sqr(Conj(Lambdax))*Sqr(Lambdax) - 100*Conj(Kappa
      )*Conj(TLambdax)*Lambdax*TKappa - 150*traceconjTYuTpYu*Conj(Lambdax)*
      TLambdax - 100*Conj(Lambdax)*Conj(TKappa)*Kappa*TLambdax));


   return beta_mHd2;
}

/**
 * Calculates the 3-loop beta function of mHd2.
 *
 * @return 3-loop beta function
 */
double CNMSSMEFTHiggs_soft_parameters::calc_beta_mHd2_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mHd2;

   beta_mHd2 = 0;


   return beta_mHd2;
}

/**
 * Calculates the 4-loop beta function of mHd2.
 *
 * @return 4-loop beta function
 */
double CNMSSMEFTHiggs_soft_parameters::calc_beta_mHd2_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mHd2;

   beta_mHd2 = 0;


   return beta_mHd2;
}

/**
 * Calculates the 5-loop beta function of mHd2.
 *
 * @return 5-loop beta function
 */
double CNMSSMEFTHiggs_soft_parameters::calc_beta_mHd2_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mHd2;

   beta_mHd2 = 0;


   return beta_mHd2;
}

} // namespace flexiblesusy
