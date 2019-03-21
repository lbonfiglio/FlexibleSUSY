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

// File generated at Thu 21 Mar 2019 14:09:28

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
 * Calculates the 1-loop beta function of MassG.
 *
 * @return 1-loop beta function
 */
double CNMSSMEFTHiggs_soft_parameters::calc_beta_MassG_1_loop(const Soft_traces& soft_traces) const
{


   double beta_MassG;

   beta_MassG = Re(-6*MassG*oneOver16PiSqr*Sqr(g3));


   return beta_MassG;
}

/**
 * Calculates the 2-loop beta function of MassG.
 *
 * @return 2-loop beta function
 */
double CNMSSMEFTHiggs_soft_parameters::calc_beta_MassG_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;


   double beta_MassG;

   beta_MassG = Re(0.4*twoLoop*Sqr(g3)*(20*traceAdjYdTYd + 20*traceAdjYuTYu -
      20*MassG*traceYdAdjYd - 20*MassG*traceYuAdjYu + 11*MassB*Sqr(g1) + 11*
      MassG*Sqr(g1) + 45*MassG*Sqr(g2) + 45*MassWB*Sqr(g2) + 140*MassG*Sqr(g3))
      );


   return beta_MassG;
}

/**
 * Calculates the 3-loop beta function of MassG.
 *
 * @return 3-loop beta function
 */
double CNMSSMEFTHiggs_soft_parameters::calc_beta_MassG_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_MassG;

   beta_MassG = 0;


   return beta_MassG;
}

/**
 * Calculates the 4-loop beta function of MassG.
 *
 * @return 4-loop beta function
 */
double CNMSSMEFTHiggs_soft_parameters::calc_beta_MassG_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_MassG;

   beta_MassG = 0;


   return beta_MassG;
}

/**
 * Calculates the 5-loop beta function of MassG.
 *
 * @return 5-loop beta function
 */
double CNMSSMEFTHiggs_soft_parameters::calc_beta_MassG_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_MassG;

   beta_MassG = 0;


   return beta_MassG;
}

} // namespace flexiblesusy
