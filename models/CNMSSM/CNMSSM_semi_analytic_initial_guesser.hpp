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

// File generated at Tue 23 Jul 2019 15:26:06

#ifndef CNMSSM_SEMI_ANALYTIC_INITIAL_GUESSER_H
#define CNMSSM_SEMI_ANALYTIC_INITIAL_GUESSER_H

#include "CNMSSM_initial_guesser.hpp"
#include "CNMSSM_semi_analytic_low_scale_constraint.hpp"
#include "CNMSSM_semi_analytic_susy_scale_constraint.hpp"
#include "CNMSSM_semi_analytic_high_scale_constraint.hpp"
#include "initial_guesser.hpp"

#include <Eigen/Core>

namespace flexiblesusy {

template <class T>
class CNMSSM;

class Semi_analytic;

/**
 * @class CNMSSM_initial_guesser<Semi_analytic>
 * @brief initial guesser for the CNMSSM
 */

template<>
class CNMSSM_initial_guesser<Semi_analytic> : public Initial_guesser {
public:
   CNMSSM_initial_guesser(CNMSSM<Semi_analytic>*,
                               const softsusy::QedQcd&,
                               CNMSSM_low_scale_constraint<Semi_analytic>&,
                               CNMSSM_susy_scale_constraint<Semi_analytic>&,
                               CNMSSM_high_scale_constraint<Semi_analytic>&);
   virtual ~CNMSSM_initial_guesser() = default;
   virtual void guess() override; ///< initial guess

   void set_running_precision(double p) { running_precision = p; }

private:
   CNMSSM<Semi_analytic>* model{nullptr}; ///< pointer to model class
   softsusy::QedQcd qedqcd{};           ///< Standard Model low-energy data
   double mu_guess{0.}; ///< guessed DR-bar mass of up-quark
   double mc_guess{0.}; ///< guessed DR-bar mass of charm-quark
   double mt_guess{0.}; ///< guessed DR-bar mass of top-quark
   double md_guess{0.}; ///< guessed DR-bar mass of down-quark
   double ms_guess{0.}; ///< guessed DR-bar mass of strange-quark
   double mb_guess{0.}; ///< guessed DR-bar mass of bottom-quark
   double me_guess{0.}; ///< guessed DR-bar mass of electron
   double mm_guess{0.}; ///< guessed DR-bar mass of muon
   double mtau_guess{0.}; ///< guessed DR-bar mass of tau
   double running_precision{1.e-3}; ///< Runge-Kutta RG running precision
   CNMSSM_low_scale_constraint<Semi_analytic>& low_constraint;
   CNMSSM_susy_scale_constraint<Semi_analytic>& susy_constraint;
   CNMSSM_high_scale_constraint<Semi_analytic>& high_constraint;
   Eigen::Matrix<std::complex<double>,3,3> upQuarksDRbar{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> downQuarksDRbar{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> downLeptonsDRbar{Eigen::Matrix<std::complex<double>,3,3>::Zero()};

   void initial_guess_low_scale_parameters();
   void initial_guess_high_scale_parameters();
   void solve_susy_parameters();
   void guess_susy_parameters();
   void guess_soft_parameters();
   void calculate_DRbar_yukawa_couplings();
   void calculate_Yu_DRbar();
   void calculate_Yd_DRbar();
   void calculate_Ye_DRbar();
   void calculate_running_SM_masses();
};

} // namespace flexiblesusy

#endif
