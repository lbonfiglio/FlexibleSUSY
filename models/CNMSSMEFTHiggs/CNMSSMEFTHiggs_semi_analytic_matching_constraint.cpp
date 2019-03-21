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

#include "CNMSSMEFTHiggs_semi_analytic_matching_constraint.hpp"
#include "CNMSSMEFTHiggs_semi_analytic_model.hpp"
#include "CNMSSMEFTHiggs_standard_model_semi_analytic_matching.hpp"
#include "wrappers.hpp"
#include "logger.hpp"
#include "ew_input.hpp"

#include <cmath>

namespace flexiblesusy {

#define DERIVEDPARAMETER(p) model->p()
#define INPUTPARAMETER(p) model->get_input().p
#define MODELPARAMETER(p) matching_condition->get_##p()
#define PHASE(p) matching_condition->get_##p()
#define BETAPARAMETER(p) beta_functions.get_##p()
#define BETAPARAMETER1(l,p) beta_functions_##l##L.get_##p()
#define BETA(p) beta_##p
#define BETA1(l,p) beta_##l##L_##p
#define LowEnergyConstant(p) Electroweak_constants::p
#define MZPole qedqcd.displayPoleMZ()
#define STANDARDDEVIATION(p) Electroweak_constants::Error_##p
#define Pole(p) matching_model.get_physical().p
#define SCALE model->get_scale()
#define THRESHOLD static_cast<int>(model->get_thresholds())
#define MODEL model
#define MODELCLASSNAME CNMSSMEFTHiggs<Semi_analytic>

CNMSSMEFTHiggs_matching_constraint<Semi_analytic>::CNMSSMEFTHiggs_matching_constraint(
   CNMSSMEFTHiggs<Semi_analytic>* model_)
   : model(model_)
{
}

void CNMSSMEFTHiggs_matching_constraint<Semi_analytic>::apply()
{
   check_model_ptr();
   check_matching_condition_ptr();

   // apply matched parameter values
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto g3 = MODELPARAMETER(g3);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto Yd = MODELPARAMETER(Yd);
   const auto Ye = MODELPARAMETER(Ye);
   const auto Yu = MODELPARAMETER(Yu);

   MODEL->set_g1(Re(g1));
   MODEL->set_g2(Re(g2));
   MODEL->set_g3(Re(g3));
   MODEL->set_vd(Re(vd));
   MODEL->set_vu(Re(vu));
   MODEL->set_Yd((Yd).real());
   MODEL->set_Ye((Ye).real());
   MODEL->set_Yu((Yu).real());

}

void CNMSSMEFTHiggs_matching_constraint<Semi_analytic>::set_model(Model* model_)
{
   model = cast_model<CNMSSMEFTHiggs<Semi_analytic>*>(model_);
}

void CNMSSMEFTHiggs_matching_constraint<Semi_analytic>::check_model_ptr() const
{
   if (!model)
      throw SetupError("CNMSSMEFTHiggs_matching_constraint<Semi_analytic>: "
                       "model pointer is zero!");
}

void CNMSSMEFTHiggs_matching_constraint<Semi_analytic>::check_matching_condition_ptr() const
{
   if (!matching_condition)
      throw SetupError("CNMSSMEFTHiggs_matching_constraint<Semi_analytic>: "
                       "matching condition pointer is zero!");
}

} // namespace flexiblesusy
