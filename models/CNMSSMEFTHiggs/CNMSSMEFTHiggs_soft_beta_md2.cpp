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

// File generated at Thu 21 Mar 2019 14:09:25

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
 * Calculates the 1-loop beta function of md2.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> CNMSSMEFTHiggs_soft_parameters::calc_beta_md2_1_loop(const Soft_traces& soft_traces) const
{
   const double Tr11 = TRACE_STRUCT.Tr11;


   Eigen::Matrix<double,3,3> beta_md2;

   beta_md2 = (oneOver16PiSqr*(4*mHd2*(Yd*Yd.adjoint()) + 4*(TYd*(TYd).adjoint(
      )) + 2*(md2*Yd*Yd.adjoint()) + 4*(Yd*mq2*Yd.adjoint()) + 2*(Yd*Yd.adjoint
      ()*md2) + 0.13333333333333333*(3.872983346207417*g1*Tr11 - 4*AbsSqr(MassB
      )*Sqr(g1) - 80*AbsSqr(MassG)*Sqr(g3))*UNITMATRIX(3))).real();


   return beta_md2;
}

/**
 * Calculates the 2-loop beta function of md2.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> CNMSSMEFTHiggs_soft_parameters::calc_beta_md2_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceconjTYdTpTYd = TRACE_STRUCT.traceconjTYdTpTYd;
   const double traceconjTYeTpTYe = TRACE_STRUCT.traceconjTYeTpTYe;
   const double tracemd2YdAdjYd = TRACE_STRUCT.tracemd2YdAdjYd;
   const double traceme2YeAdjYe = TRACE_STRUCT.traceme2YeAdjYe;
   const double traceml2AdjYeYe = TRACE_STRUCT.traceml2AdjYeYe;
   const double tracemq2AdjYdYd = TRACE_STRUCT.tracemq2AdjYdYd;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceconjTYdTpYd = TRACE_STRUCT.traceconjTYdTpYd;
   const double traceconjTYeTpYe = TRACE_STRUCT.traceconjTYeTpYe;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr23 = TRACE_STRUCT.Tr23;


   Eigen::Matrix<double,3,3> beta_md2;

   beta_md2 = (twoLoop*(0.8*(-15*traceconjTYdTpTYd - 5*traceconjTYeTpTYe - 15*
      tracemd2YdAdjYd - 5*traceme2YeAdjYe - 5*traceml2AdjYeYe - 15*
      tracemq2AdjYdYd - 30*mHd2*traceYdAdjYd - 10*mHd2*traceYeAdjYe - 10*mHd2*
      AbsSqr(Lambdax) - 5*mHu2*AbsSqr(Lambdax) - 5*ms2*AbsSqr(Lambdax) - 5*
      AbsSqr(TLambdax) + mHd2*Sqr(g1) + 2*AbsSqr(MassB)*Sqr(g1) + 15*mHd2*Sqr(
      g2) + 30*AbsSqr(MassWB)*Sqr(g2))*(Yd*Yd.adjoint()) - 0.8*(15*
      traceAdjYdTYd + 5*traceAdjYeTYe + MassB*Sqr(g1) + 15*MassWB*Sqr(g2) + 5*
      Conj(Lambdax)*TLambdax)*(Yd*(TYd).adjoint()) - 0.8*(15*traceconjTYdTpYd +
      5*traceconjTYeTpYe + 5*Conj(TLambdax)*Lambdax + Conj(MassB)*Sqr(g1) + 15*
      Conj(MassWB)*Sqr(g2))*(TYd*Yd.adjoint()) + 0.8*(-15*traceYdAdjYd - 5*
      traceYeAdjYe - 5*AbsSqr(Lambdax) + Sqr(g1) + 15*Sqr(g2))*(TYd*(TYd).
      adjoint()) + 0.4*(-15*traceYdAdjYd - 5*traceYeAdjYe - 5*AbsSqr(Lambdax) +
      Sqr(g1) + 15*Sqr(g2))*(md2*Yd*Yd.adjoint()) + 0.8*(-15*traceYdAdjYd - 5*
      traceYeAdjYe - 5*AbsSqr(Lambdax) + Sqr(g1) + 15*Sqr(g2))*(Yd*mq2*Yd.
      adjoint()) + 0.4*(-15*traceYdAdjYd - 5*traceYeAdjYe - 5*AbsSqr(Lambdax) +
      Sqr(g1) + 15*Sqr(g2))*(Yd*Yd.adjoint()*md2) - 8*mHd2*(Yd*Yd.adjoint()*Yd*
      Yd.adjoint()) - 4*(Yd*Yd.adjoint()*TYd*(TYd).adjoint()) - 4*(mHd2 + mHu2)
      *(Yd*Yu.adjoint()*Yu*Yd.adjoint()) - 4*(Yd*Yu.adjoint()*TYu*(TYd).adjoint
      ()) - 4*(Yd*(TYd).adjoint()*TYd*Yd.adjoint()) - 4*(Yd*(TYu).adjoint()*TYu
      *Yd.adjoint()) - 4*(TYd*Yd.adjoint()*Yd*(TYd).adjoint()) - 4*(TYd*Yu.
      adjoint()*Yu*(TYd).adjoint()) - 4*(TYd*(TYd).adjoint()*Yd*Yd.adjoint()) -
      4*(TYd*(TYu).adjoint()*Yu*Yd.adjoint()) - 2*(md2*Yd*Yd.adjoint()*Yd*Yd.
      adjoint()) - 2*(md2*Yd*Yu.adjoint()*Yu*Yd.adjoint()) - 4*(Yd*mq2*Yd.
      adjoint()*Yd*Yd.adjoint()) - 4*(Yd*mq2*Yu.adjoint()*Yu*Yd.adjoint()) - 4*
      (Yd*Yd.adjoint()*md2*Yd*Yd.adjoint()) - 4*(Yd*Yd.adjoint()*Yd*mq2*Yd.
      adjoint()) - 2*(Yd*Yd.adjoint()*Yd*Yd.adjoint()*md2) - 4*(Yd*Yu.adjoint()
      *mu2*Yu*Yd.adjoint()) - 4*(Yd*Yu.adjoint()*Yu*mq2*Yd.adjoint()) - 2*(Yd*
      Yu.adjoint()*Yu*Yd.adjoint()*md2) + 0.035555555555555556*(
      58.09475019311125*g1*Tr31 + 303*AbsSqr(MassB)*Quad(g1) + 300*Tr23*Quad(g3
      ) - 1200*AbsSqr(MassG)*Quad(g3) + 15*Tr2U111*Sqr(g1) + 80*AbsSqr(MassB)*
      Sqr(g1)*Sqr(g3) + 80*AbsSqr(MassG)*Sqr(g1)*Sqr(g3) + 40*MassG*Conj(MassB)
      *Sqr(g1)*Sqr(g3) + 40*MassB*Conj(MassG)*Sqr(g1)*Sqr(g3))*UNITMATRIX(3))).
      real();


   return beta_md2;
}

/**
 * Calculates the 3-loop beta function of md2.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> CNMSSMEFTHiggs_soft_parameters::calc_beta_md2_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_md2;

   beta_md2 = ZEROMATRIX(3,3);


   return beta_md2;
}

/**
 * Calculates the 4-loop beta function of md2.
 *
 * @return 4-loop beta function
 */
Eigen::Matrix<double,3,3> CNMSSMEFTHiggs_soft_parameters::calc_beta_md2_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_md2;

   beta_md2 = ZEROMATRIX(3,3);


   return beta_md2;
}

/**
 * Calculates the 5-loop beta function of md2.
 *
 * @return 5-loop beta function
 */
Eigen::Matrix<double,3,3> CNMSSMEFTHiggs_soft_parameters::calc_beta_md2_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_md2;

   beta_md2 = ZEROMATRIX(3,3);


   return beta_md2;
}

} // namespace flexiblesusy
