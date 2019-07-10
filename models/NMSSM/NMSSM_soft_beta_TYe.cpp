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

// File generated at Wed 10 Jul 2019 12:33:42

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
 * Calculates the 1-loop beta function of TYe.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> NMSSM_soft_parameters::calc_beta_TYe_1_loop(const Soft_traces& soft_traces) const
{
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   Eigen::Matrix<double,3,3> beta_TYe;

   beta_TYe = (oneOver16PiSqr*(0.2*(30*traceAdjYdTYd*Ye + 10*traceAdjYeTYe*Ye +
      18*MassB*Ye*Sqr(g1) + 30*MassWB*Ye*Sqr(g2) + 15*traceYdAdjYd*TYe + 5*
      traceYeAdjYe*TYe + 5*AbsSqr(Lambdax)*TYe - 9*Sqr(g1)*TYe - 15*Sqr(g2)*TYe
       + 10*Ye*Conj(Lambdax)*TLambdax) + 4*(Ye*Ye.adjoint()*TYe) + 5*(TYe*Ye.
      adjoint()*Ye))).real();


   return beta_TYe;
}

/**
 * Calculates the 2-loop beta function of TYe.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> NMSSM_soft_parameters::calc_beta_TYe_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceYdAdjYdTYdAdjYd = TRACE_STRUCT.traceYdAdjYdTYdAdjYd;
   const double traceYdAdjYuTYuAdjYd = TRACE_STRUCT.traceYdAdjYuTYuAdjYd;
   const double traceYeAdjYeTYeAdjYe = TRACE_STRUCT.traceYeAdjYeTYeAdjYe;
   const double traceYuAdjYdTYdAdjYu = TRACE_STRUCT.traceYuAdjYdTYdAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;


   Eigen::Matrix<double,3,3> beta_TYe;

   beta_TYe = (twoLoop*(0.1*(-360*traceYdAdjYdTYdAdjYd*Ye - 60*
      traceYdAdjYuTYuAdjYd*Ye - 120*traceYeAdjYeTYeAdjYe*Ye - 60*
      traceYuAdjYdTYdAdjYu*Ye - 60*traceAdjYuTYu*Ye*AbsSqr(Lambdax) - 540*MassB
      *Ye*Quad(g1) - 300*MassWB*Ye*Quad(g2) - 8*traceAdjYdTYd*Ye*Sqr(g1) + 24*
      traceAdjYeTYe*Ye*Sqr(g1) + 8*MassB*traceYdAdjYd*Ye*Sqr(g1) - 24*MassB*
      traceYeAdjYe*Ye*Sqr(g1) - 36*MassB*Ye*Sqr(g1)*Sqr(g2) - 36*MassWB*Ye*Sqr(
      g1)*Sqr(g2) + 320*traceAdjYdTYd*Ye*Sqr(g3) - 320*MassG*traceYdAdjYd*Ye*
      Sqr(g3) - 90*traceYdAdjYdYdAdjYd*TYe - 30*traceYdAdjYuYuAdjYd*TYe - 30*
      traceYeAdjYeYeAdjYe*TYe - 30*traceYuAdjYu*AbsSqr(Lambdax)*TYe - 20*AbsSqr
      (Kappa)*AbsSqr(Lambdax)*TYe + 135*Quad(g1)*TYe + 75*Quad(g2)*TYe - 4*
      traceYdAdjYd*Sqr(g1)*TYe + 12*traceYeAdjYe*Sqr(g1)*TYe + 18*Sqr(g1)*Sqr(
      g2)*TYe + 160*traceYdAdjYd*Sqr(g3)*TYe - 30*Sqr(Conj(Lambdax))*Sqr(
      Lambdax)*TYe - 40*Ye*AbsSqr(Lambdax)*Conj(Kappa)*TKappa - 60*traceYuAdjYu
      *Ye*Conj(Lambdax)*TLambdax - 40*Ye*AbsSqr(Kappa)*Conj(Lambdax)*TLambdax -
      120*Ye*Lambdax*Sqr(Conj(Lambdax))*TLambdax) - 6*(3*traceAdjYdTYd +
      traceAdjYeTYe + 2*MassWB*Sqr(g2) + Conj(Lambdax)*TLambdax)*(Ye*Ye.adjoint
      ()*Ye) + 0.4*(-30*traceYdAdjYd - 10*traceYeAdjYe - 10*AbsSqr(Lambdax) + 3
      *Sqr(g1) + 15*Sqr(g2))*(Ye*Ye.adjoint()*TYe) + 0.2*(-75*traceYdAdjYd - 25
      *traceYeAdjYe - 25*AbsSqr(Lambdax) - 6*Sqr(g1) + 60*Sqr(g2))*(TYe*Ye.
      adjoint()*Ye) - 6*(Ye*Ye.adjoint()*Ye*Ye.adjoint()*TYe) - 8*(Ye*Ye.
      adjoint()*TYe*Ye.adjoint()*Ye) - 6*(TYe*Ye.adjoint()*Ye*Ye.adjoint()*Ye))
      ).real();


   return beta_TYe;
}

/**
 * Calculates the 3-loop beta function of TYe.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> NMSSM_soft_parameters::calc_beta_TYe_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_TYe;

   beta_TYe = ZEROMATRIX(3,3);


   return beta_TYe;
}

/**
 * Calculates the 4-loop beta function of TYe.
 *
 * @return 4-loop beta function
 */
Eigen::Matrix<double,3,3> NMSSM_soft_parameters::calc_beta_TYe_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_TYe;

   beta_TYe = ZEROMATRIX(3,3);


   return beta_TYe;
}

/**
 * Calculates the 5-loop beta function of TYe.
 *
 * @return 5-loop beta function
 */
Eigen::Matrix<double,3,3> NMSSM_soft_parameters::calc_beta_TYe_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_TYe;

   beta_TYe = ZEROMATRIX(3,3);


   return beta_TYe;
}

} // namespace flexiblesusy
