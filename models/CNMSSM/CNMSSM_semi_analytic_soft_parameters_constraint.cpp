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

// File generated at Tue 23 Jul 2019 15:26:05

#include "CNMSSM_semi_analytic_soft_parameters_constraint.hpp"
#include "CNMSSM_semi_analytic_model.hpp"
#include "wrappers.hpp"
#include "logger.hpp"
#include "ew_input.hpp"
#include "gsl_utils.hpp"
#include "minimizer.hpp"
#include "root_finder.hpp"
#include "threshold_loop_functions.hpp"

#include <cmath>

namespace flexiblesusy {

#define DERIVEDPARAMETER(p) model->p()
#define INPUTPARAMETER(p) model->get_input().p
#define MODELPARAMETER(p) model->get_##p()
#define PHASE(p) model->get_##p()
#define BETAPARAMETER(p) beta_functions.get_##p()
#define BETAPARAMETER1(l,p) beta_functions_##l##L.get_##p()
#define BETA(p) beta_##p
#define BETA1(l,p) beta_##l##L_##p
#define LowEnergyConstant(p) Electroweak_constants::p
#define MZPole qedqcd.displayPoleMZ()
#define STANDARDDEVIATION(p) Electroweak_constants::Error_##p
#define Pole(p) model->get_physical().p
#define SCALE model->get_scale()
#define THRESHOLD static_cast<int>(model->get_thresholds())
#define MODEL model
#define MODELCLASSNAME CNMSSM<Semi_analytic>

CNMSSM_soft_parameters_constraint<Semi_analytic>::CNMSSM_soft_parameters_constraint(
   CNMSSM<Semi_analytic>* model_, const softsusy::QedQcd& qedqcd_)
   : model(model_)
   , qedqcd(qedqcd_)
{
   initialize();
}

void CNMSSM_soft_parameters_constraint<Semi_analytic>::apply()
{
   check_model_ptr();

   


   const double boundary_scale = get_boundary_scale();
   model->calculate_semi_analytic_solutions(boundary_scale);
   model->calculate_DRbar_masses();
   update_scale();

   // apply user-defined constraints
   
   MODEL->solve_ewsb();


}

   const CNMSSM_input_parameters& CNMSSM_soft_parameters_constraint<Semi_analytic>::get_input_parameters() const
{
   check_model_ptr();

   return model->get_input();
}

void CNMSSM_soft_parameters_constraint<Semi_analytic>::set_model(Model* model_)
{
   model = cast_model<CNMSSM<Semi_analytic>*>(model_);
}

void CNMSSM_soft_parameters_constraint<Semi_analytic>::clear()
{
   scale = 0.;
   initial_scale_guess = 0.;
   model = nullptr;
   qedqcd = softsusy::QedQcd();
}

void CNMSSM_soft_parameters_constraint<Semi_analytic>::initialize()
{
   check_model_ptr();

   const auto Azero = INPUTPARAMETER(Azero);
   const auto m12 = INPUTPARAMETER(m12);

   initial_scale_guess = Sqrt(-3*Azero*m12 + Sqr(Azero) + 14*Sqr(m12));

   scale = initial_scale_guess;
}

void CNMSSM_soft_parameters_constraint<Semi_analytic>::update_scale()
{
   check_model_ptr();

   const auto ZU = MODELPARAMETER(ZU);
   const auto MSu = MODELPARAMETER(MSu);

   scale = Sqrt(Power(MSu(0),Sqr(Abs(ZU(0,2))) + Sqr(Abs(ZU(0,5))))*Power(MSu(1),
      Sqr(Abs(ZU(1,2))) + Sqr(Abs(ZU(1,5))))*Power(MSu(2),Sqr(Abs(ZU(2,2))) + Sqr(
      Abs(ZU(2,5))))*Power(MSu(3),Sqr(Abs(ZU(3,2))) + Sqr(Abs(ZU(3,5))))*Power(MSu
      (4),Sqr(Abs(ZU(4,2))) + Sqr(Abs(ZU(4,5))))*Power(MSu(5),Sqr(Abs(ZU(5,2))) +
      Sqr(Abs(ZU(5,5)))));


}

double CNMSSM_soft_parameters_constraint<Semi_analytic>::get_boundary_scale() const
{
   if (!boundary_scale_getter)
      return scale;

   return boundary_scale_getter();
}

void CNMSSM_soft_parameters_constraint<Semi_analytic>::set_boundary_scale(const Scale_getter& scale_getter_)
{
   boundary_scale_getter = scale_getter_;
}

void CNMSSM_soft_parameters_constraint<Semi_analytic>::set_boundary_scale(Scale_getter&& scale_getter_)
{
   boundary_scale_getter = std::move(scale_getter_);
}

void CNMSSM_soft_parameters_constraint<Semi_analytic>::check_model_ptr() const
{
   if (!model)
      throw SetupError("CNMSSM_soft_parameters_constraint<Semi_analytic>: "
                       "model pointer is zero!");
}

} // namespace flexiblesusy
