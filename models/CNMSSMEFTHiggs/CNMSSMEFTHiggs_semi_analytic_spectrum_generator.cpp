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

// File generated at Thu 21 Mar 2019 14:21:46

#include "CNMSSMEFTHiggs_semi_analytic_spectrum_generator.hpp"
#include "CNMSSMEFTHiggs_input_parameters.hpp"
#include "CNMSSMEFTHiggs_semi_analytic_convergence_tester.hpp"
#include "CNMSSMEFTHiggs_semi_analytic_ewsb_solver.hpp"
#include "CNMSSMEFTHiggs_semi_analytic_high_scale_constraint.hpp"
#include "CNMSSMEFTHiggs_semi_analytic_initial_guesser.hpp"
#include "CNMSSMEFTHiggs_semi_analytic_matching_constraint.hpp"
#include "CNMSSMEFTHiggs_semi_analytic_soft_parameters_constraint.hpp"
#include "CNMSSMEFTHiggs_semi_analytic_susy_convergence_tester.hpp"
#include "CNMSSMEFTHiggs_semi_analytic_susy_scale_constraint.hpp"
#include "CNMSSMEFTHiggs_standard_model_matching.hpp"
#include "CNMSSMEFTHiggs_standard_model_semi_analytic_matching.hpp"
#include "standard_model_two_scale_convergence_tester.hpp"
#include "standard_model_two_scale_low_scale_constraint.hpp"

#include "composite_convergence_tester.hpp"
#include "error.hpp"
#include "lowe.h"
#include "numerics2.hpp"
#include "semi_analytic_solver.hpp"
#include "two_scale_running_precision.hpp"

#include <algorithm>
#include <limits>

namespace flexiblesusy {

double CNMSSMEFTHiggs_spectrum_generator<Semi_analytic>::get_pole_mass_scale(double susy_scale) const
{
   return settings.get(Spectrum_generator_settings::pole_mass_scale) != 0. ?
      settings.get(Spectrum_generator_settings::pole_mass_scale) :
      susy_scale;
}

double CNMSSMEFTHiggs_spectrum_generator<Semi_analytic>::get_eft_pole_mass_scale(double susy_scale, double Mt) const
{
   double Q_higgs = settings.get(
      Spectrum_generator_settings::eft_pole_mass_scale);

   if (Q_higgs == 0.)
      Q_higgs = std::min(susy_scale, Mt);

   return Q_higgs;
}

/**
 * @brief Run's the RG solver with the given input parameters
 *
 * This function sets up the RG solver using a high-scale, susy-scale
 * and low-scale constraint.  Afterwards the solver is run until
 * convergence is reached or an error occours.  Finally the particle
 * spectrum (pole masses) is calculated.
 *
 * @param qedqcd Standard Model input parameters
 * @param input model input parameters
 */
void CNMSSMEFTHiggs_spectrum_generator<Semi_analytic>::run_except(const softsusy::QedQcd& qedqcd,
                                const CNMSSMEFTHiggs_input_parameters& input)
{
   VERBOSE_MSG("Solving BVP using semi-analytic solver");

   problems.set_bvp_solver_problems({ BVP_solver_problems("SemiAnalyticSolver") });

   auto& model = this->model;
   model.clear();
   model.set_input_parameters(input);
   model.do_calculate_sm_pole_masses(
      settings.get(Spectrum_generator_settings::calculate_sm_masses));
   model.do_calculate_bsm_pole_masses(
      settings.get(Spectrum_generator_settings::calculate_bsm_masses));
   model.do_force_output(
      settings.get(Spectrum_generator_settings::force_output));
   model.set_loops(
      settings.get(Spectrum_generator_settings::beta_loop_order));
   model.set_thresholds(
      settings.get(
         Spectrum_generator_settings::threshold_corrections_loop_order));
   model.set_zero_threshold(
      settings.get(Spectrum_generator_settings::beta_zero_threshold));

   eft.clear();
   eft.do_force_output(
      settings.get(Spectrum_generator_settings::force_output));
   eft.set_loops(
      settings.get(Spectrum_generator_settings::beta_loop_order));
   eft.set_thresholds(
      settings.get(
         Spectrum_generator_settings::threshold_corrections_loop_order));
   eft.set_zero_threshold(
      settings.get(Spectrum_generator_settings::beta_zero_threshold));
   eft.set_pole_mass_loop_order(this->model.get_pole_mass_loop_order());
   eft.set_ewsb_loop_order(this->model.get_ewsb_loop_order());
   eft.set_ewsb_iteration_precision(this->model.get_ewsb_iteration_precision());
   eft.set_loop_corrections(this->model.get_loop_corrections());
   eft.set_threshold_corrections(this->model.get_threshold_corrections());

   CNMSSMEFTHiggs_semi_analytic_solutions& solutions(
      model.get_semi_analytic_solutions());

   CNMSSMEFTHiggs_ewsb_solver<Semi_analytic> ewsb_solver;
   ewsb_solver.set_semi_analytic_solutions(&solutions);
   model.set_ewsb_solver(
      std::make_shared<CNMSSMEFTHiggs_ewsb_solver<Semi_analytic> >(ewsb_solver));

   CNMSSMEFTHiggs_high_scale_constraint<Semi_analytic> high_scale_constraint(&model);
   CNMSSMEFTHiggs_susy_scale_constraint<Semi_analytic> susy_scale_constraint(&model, qedqcd);
   CNMSSMEFTHiggs_soft_parameters_constraint<Semi_analytic> soft_constraint(&model, qedqcd);
   CNMSSMEFTHiggs_matching_constraint<Semi_analytic> matching_constraint(&model);
   standard_model::Standard_model_low_scale_constraint<Two_scale> low_scale_constraint(&eft, qedqcd);

   // note: to avoid large logarithms the downwards matching loop order
   // is used for both matching conditions
   const int matching_loop_order_up = settings.get(
      Spectrum_generator_settings::eft_matching_loop_order_down);
   const int matching_loop_order_down = settings.get(
      Spectrum_generator_settings::eft_matching_loop_order_down);
   const int index = settings.get(
      Spectrum_generator_settings::eft_higgs_index);
   const auto scale_getter = [this,&susy_scale_constraint] () {
      return susy_scale_constraint.get_scale(); };
   const auto boundary_scale_getter = [this,&high_scale_constraint] () {
      return high_scale_constraint.get_scale(); };

   CNMSSMEFTHiggs_standard_model_matching_up<Semi_analytic> matching_up(
      &eft, &model, scale_getter, matching_loop_order_up, index);
   CNMSSMEFTHiggs_standard_model_matching_down<Semi_analytic> matching_down(
      &eft, &model, scale_getter, matching_loop_order_down, index);

   matching_up.set_scale(
      settings.get(Spectrum_generator_settings::eft_matching_scale));
   matching_down.set_scale(
      settings.get(Spectrum_generator_settings::eft_matching_scale));

   soft_constraint.set_boundary_scale(boundary_scale_getter);
   susy_scale_constraint.set_scale(
      [this,&soft_constraint] () {
         return soft_constraint.get_scale(); });

   matching_constraint.set_matching_condition(&matching_up);
   matching_up.set_boundary_scale(boundary_scale_getter);
   matching_up.set_matching_constraint(&matching_constraint);
   matching_down.set_boundary_scale(boundary_scale_getter);

   high_scale_constraint.initialize();
   susy_scale_constraint.initialize();
   soft_constraint.initialize();
   low_scale_constraint.initialize();

   // convergence tester for CNMSSMEFTHiggs
   CNMSSMEFTHiggs_susy_convergence_tester<Semi_analytic> inner_ct(
      &model, settings.get(Spectrum_generator_settings::precision));
   CNMSSMEFTHiggs_convergence_tester<Semi_analytic> outer_ct(
      &model, settings.get(Spectrum_generator_settings::precision),
      [this,&susy_scale_constraint](){
         return get_pole_mass_scale(susy_scale_constraint.get_scale()); });

   // convergence tester for Standard Model
   const double Mt = qedqcd.displayPoleMt();
   standard_model::Standard_model_convergence_tester<Two_scale> eft_ct(
      &eft, settings.get(Spectrum_generator_settings::precision),
      [this,&susy_scale_constraint,Mt](){
         return get_eft_pole_mass_scale(susy_scale_constraint.get_scale(), Mt);
      });

   if (settings.get(Spectrum_generator_settings::max_iterations) > 0) {
      inner_ct.set_max_iterations(
         settings.get(Spectrum_generator_settings::max_iterations));
      outer_ct.set_max_iterations(
         settings.get(Spectrum_generator_settings::max_iterations));
      eft_ct.set_max_iterations(
         settings.get(Spectrum_generator_settings::max_iterations));
   }

   Composite_convergence_tester cct;
   cct.add_convergence_tester(&outer_ct);
   cct.add_convergence_tester(&eft_ct);

   CNMSSMEFTHiggs_standard_model_initial_guesser<Semi_analytic> initial_guesser(&model, &eft, qedqcd,
                                                                             low_scale_constraint,
                                                                             susy_scale_constraint,
                                                                             high_scale_constraint,
                                                                             matching_constraint,
                                                                             matching_up);

   Two_scale_increasing_precision precision(
      10.0, settings.get(Spectrum_generator_settings::precision));

   RGFlow<Semi_analytic> solver;
   solver.set_inner_convergence_tester(&inner_ct);
   solver.set_outer_convergence_tester(&cct);
   solver.set_running_precision(&precision);
   solver.set_initial_guesser(&initial_guesser);
   solver.add_inner(&matching_constraint, &model);
   solver.add_inner(&high_scale_constraint, &model);
   solver.add_inner(&susy_scale_constraint, &model);
   solver.add_outer(&soft_constraint, &model);
   solver.add_outer(&matching_down, &model, &eft);
   solver.add_outer(&low_scale_constraint, &eft);
   solver.add_outer(&matching_up, &eft, &model);

   high_scale = susy_scale = low_scale = 0.;
   reached_precision = std::numeric_limits<double>::infinity();

   solver.solve();

   // impose low-scale constraint one last time
   eft.run_to(low_scale_constraint.get_scale());
   low_scale_constraint.apply();

   high_scale = high_scale_constraint.get_scale();
   susy_scale = susy_scale_constraint.get_scale();
   low_scale  = low_scale_constraint.get_scale();
   reached_precision = std::max(outer_ct.get_current_accuracy(),
                                eft_ct.get_current_accuracy());

   calculate_spectrum(Mt, low_scale_constraint.get_sm_parameters().displayPoleMW());

   // run to output scale (if scale > 0)
   if (!is_zero(parameter_output_scale)) {
      model.run_to(parameter_output_scale);
      eft.run_to(parameter_output_scale);
   }
}

/**
 * Create a text file which contains the values of all model
 * parameters at all scales between the low-scale and the high-scale.
 *
 * @param filename name of output file
 */
void CNMSSMEFTHiggs_spectrum_generator<Semi_analytic>::write_running_couplings(
   const std::string& filename) const
{
   CNMSSMEFTHiggs_spectrum_generator_interface<
      Semi_analytic, standard_model::StandardModel<Two_scale> >::write_running_couplings(
      filename, get_low_scale(), get_high_scale());
}

void CNMSSMEFTHiggs_spectrum_generator<Semi_analytic>::calculate_spectrum(double Mt, double MW)
{
   model.run_to(get_pole_mass_scale(susy_scale));
   model.calculate_semi_analytic_solutions(get_high_scale());
   model.solve_ewsb();
   model.calculate_spectrum();

   eft.run_to(get_eft_pole_mass_scale(susy_scale, Mt));

   // computation of pole mass spectrum in the SM
   const int eft_pole_loops = eft.get_pole_mass_loop_order();
   const int eft_ewsb_loops = eft.get_ewsb_loop_order();

   eft.calculate_DRbar_masses();
   eft.solve_ewsb();
   eft.calculate_spectrum();

   eft.set_pole_mass_loop_order(eft_pole_loops);
   eft.set_ewsb_loop_order(eft_ewsb_loops);

   const int index = settings.get(Spectrum_generator_settings::eft_higgs_index);

   model.get_physical().Mhh(index) = eft.get_physical().Mhh;
   model.get_physical().MVZ = eft.get_physical().MVZ;
   model.get_physical().MVWm = MW;
   this->model.get_physical().MFu(0) = eft.get_physical().MFu(0);
   this->model.get_physical().MFu(1) = eft.get_physical().MFu(1);
   this->model.get_physical().MFu(2) = eft.get_physical().MFu(2);
   this->model.get_physical().MFd(0) = eft.get_physical().MFd(0);
   this->model.get_physical().MFd(1) = eft.get_physical().MFd(1);
   this->model.get_physical().MFd(2) = eft.get_physical().MFd(2);
   this->model.get_physical().MFe(0) = eft.get_physical().MFe(0);
   this->model.get_physical().MFe(1) = eft.get_physical().MFe(1);
   this->model.get_physical().MFe(2) = eft.get_physical().MFe(2);

   if (eft.get_problems().is_running_tachyon(standard_model_info::hh))
      model.get_problems().flag_running_tachyon(CNMSSMEFTHiggs_info::hh);
   if (eft.get_problems().is_pole_tachyon(standard_model_info::hh))
      model.get_problems().flag_pole_tachyon(CNMSSMEFTHiggs_info::hh);
   if (eft.get_problems().is_running_tachyon(standard_model_info::VZ))
      model.get_problems().flag_running_tachyon(CNMSSMEFTHiggs_info::VZ);
   if (eft.get_problems().is_pole_tachyon(standard_model_info::VZ))
      model.get_problems().flag_pole_tachyon(CNMSSMEFTHiggs_info::VZ);
   if (eft.get_problems().is_running_tachyon(standard_model_info::VWp))
      model.get_problems().flag_running_tachyon(CNMSSMEFTHiggs_info::VWm);
   if (eft.get_problems().is_pole_tachyon(standard_model_info::VWp))
      model.get_problems().flag_pole_tachyon(CNMSSMEFTHiggs_info::VWm);
}

} // namespace flexiblesusy
