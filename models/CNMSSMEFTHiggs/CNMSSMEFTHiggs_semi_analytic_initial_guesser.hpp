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

// File generated at Thu 21 Mar 2019 14:21:43

#ifndef CNMSSMEFTHiggs_STANDARD_MODEL_SEMI_ANALYTIC_INITIAL_GUESSER_H
#define CNMSSMEFTHiggs_STANDARD_MODEL_SEMI_ANALYTIC_INITIAL_GUESSER_H

#include "CNMSSMEFTHiggs_initial_guesser.hpp"
#include "CNMSSMEFTHiggs_semi_analytic_susy_scale_constraint.hpp"
#include "CNMSSMEFTHiggs_semi_analytic_high_scale_constraint.hpp"
#include "CNMSSMEFTHiggs_semi_analytic_matching_constraint.hpp"
#include "CNMSSMEFTHiggs_standard_model_semi_analytic_matching.hpp"
#include "standard_model_two_scale_low_scale_constraint.hpp"
#include "initial_guesser.hpp"
#include "lowe.h"

#include <sstream>

namespace flexiblesusy {

class Semi_analytic;
class Two_scale;

template <class T>
class CNMSSMEFTHiggs;

template <class T>
class StandardModel;

template <class T>
class CNMSSMEFTHiggs_standard_model_initial_guesser;

/**
 * @class CNMSSMEFTHiggs_standard_model_initial_guesser<Semi_analytic>
 * @brief initial guesser for the CNMSSMEFTHiggs tower
 */

template<>
class CNMSSMEFTHiggs_standard_model_initial_guesser<Semi_analytic> : public Initial_guesser {
public:
   CNMSSMEFTHiggs_standard_model_initial_guesser(CNMSSMEFTHiggs<Semi_analytic>*,
                               standard_model::StandardModel<Two_scale>*,
                               const softsusy::QedQcd&,
                               const standard_model::Standard_model_low_scale_constraint<Two_scale>&,
                               CNMSSMEFTHiggs_susy_scale_constraint<Semi_analytic>&,
                               CNMSSMEFTHiggs_high_scale_constraint<Semi_analytic>&,
                               CNMSSMEFTHiggs_matching_constraint<Semi_analytic>&,
                               CNMSSMEFTHiggs_standard_model_matching_up<Semi_analytic>&);
   virtual ~CNMSSMEFTHiggs_standard_model_initial_guesser() = default;
   virtual void guess() override; ///< initial guess

   void set_running_precision(double p) { running_precision = p; }

private:
   CNMSSMEFTHiggs<Semi_analytic>* model{nullptr}; ///< pointer to model class
   standard_model::StandardModel<Two_scale>* eft{nullptr}; ///< pointer to effective low energy model
   softsusy::QedQcd qedqcd{}; ///< Standard Model low-energy data
   double mu_guess{0.}; ///< guessed DR-bar mass of up-quark
   double mc_guess{0.}; ///< guessed DR-bar mass of charm-quark
   double mt_guess{0.}; ///< guessed DR-bar mass of top-quark
   double md_guess{0.}; ///< guessed DR-bar mass of down-quark
   double ms_guess{0.}; ///< guessed DR-bar mass of strange-quark
   double mb_guess{0.}; ///< guessed DR-bar mass of bottom-quark
   double me_guess{0.}; ///< guessed DR-bar mass of electron
   double mm_guess{0.}; ///< guessed DR-bar mass of muon
   double mtau_guess{0.}; ///< guessed DR-bar mass of tau
   double running_precision{1.0e-3}; ///< Runge-Kutta RG running precision
   standard_model::Standard_model_low_scale_constraint<Two_scale> low_constraint{};
   CNMSSMEFTHiggs_susy_scale_constraint<Semi_analytic>& susy_constraint;
   CNMSSMEFTHiggs_high_scale_constraint<Semi_analytic>& high_constraint;
   CNMSSMEFTHiggs_matching_constraint<Semi_analytic>& matching_constraint;
   CNMSSMEFTHiggs_standard_model_matching_up<Semi_analytic>& matching_up;

   void guess_eft_parameters();
   void guess_model_parameters();
   void calculate_DRbar_yukawa_couplings();
   void calculate_Yu_DRbar();
   void calculate_Yd_DRbar();
   void calculate_Ye_DRbar();
   void initial_guess_susy_scale_parameters();
   void initial_guess_high_scale_parameters();
   void solve_susy_parameters();
   void guess_soft_parameters();
};

} // namespace flexiblesusy

#endif
