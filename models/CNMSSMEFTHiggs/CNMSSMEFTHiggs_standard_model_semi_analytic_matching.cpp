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

#include "CNMSSMEFTHiggs_standard_model_semi_analytic_matching.hpp"
#include "CNMSSMEFTHiggs_standard_model_matching.hpp"
#include "CNMSSMEFTHiggs_semi_analytic_matching_constraint.hpp"
#include "CNMSSMEFTHiggs_semi_analytic_model.hpp"
#include "standard_model_two_scale_model.hpp"
#include "error.hpp"

namespace flexiblesusy {

#define CLASSNAME CNMSSMEFTHiggs_standard_model_matching_up<Semi_analytic>

CLASSNAME::CNMSSMEFTHiggs_standard_model_matching_up(
   standard_model::StandardModel<Two_scale>* low_,
   CNMSSMEFTHiggs<Semi_analytic>* high_,
   const Scale_getter& scale_getter_,
   int loop_order_,
   int higgs_idx_)
   : model(high_)
   , eft(low_)
   , scale_getter(scale_getter_)
   , scale(0.)
   , higgs_idx(higgs_idx_)
{
   set_loop_order(loop_order_);
}

double CLASSNAME::get_scale() const
{
   if (scale != 0.)
      return scale;

   if (!scale_getter)
      throw SetupError("Scale getter is not set!");

   return scale_getter();
}

void CLASSNAME::set_models(Model* low, Model* high)
{
   eft = cast_model<standard_model::StandardModel<Two_scale>*>(low);
   model = cast_model<CNMSSMEFTHiggs<Semi_analytic>*>(high);
}

void CLASSNAME::match()
{
   if (!model || !eft)
      throw SetupError("Model pointer in matching class is NULL!");

   eft->run_to(get_scale());
   model->run_to(get_scale());

   if (matching_constraint)
      matching_constraint->set_scale(get_scale());

   model->calculate_semi_analytic_solutions(get_boundary_scale());

   if (model->get_thresholds() && loop_order)
      CNMSSMEFTHiggs_standard_model_matching::match_low_to_high_scale_model(*model, *eft, loop_order, higgs_idx);
   else
      CNMSSMEFTHiggs_standard_model_matching::match_low_to_high_scale_model_tree_level(*model, *eft);

   save_matched_parameters();
}

void CLASSNAME::match_tree_level()
{
   if (!model || !eft)
      throw SetupError("Model pointer in matching class is NULL!");

   eft->run_to(get_scale());
   model->run_to(get_scale());

   if (matching_constraint)
      matching_constraint->set_scale(get_scale());

   model->calculate_semi_analytic_solutions(get_boundary_scale());

   CNMSSMEFTHiggs_standard_model_matching::match_low_to_high_scale_model_tree_level(*model, *eft);

   save_matched_parameters();
}

void CLASSNAME::save_matched_parameters()
{
   g1 = model->get_g1();
   g2 = model->get_g2();
   g3 = model->get_g3();
   vd = model->get_vd();
   vu = model->get_vu();
   Yd = model->get_Yd();
   Ye = model->get_Ye();
   Yu = model->get_Yu();

}

int CLASSNAME::get_higgs_index() const
{
   return higgs_idx;
}

int CLASSNAME::get_loop_order() const
{
   return loop_order;
}

double CLASSNAME::get_boundary_scale() const
{
   if (!boundary_scale_getter)
      return get_scale();
   return boundary_scale_getter();
}

void CLASSNAME::set_higgs_index(int idx)
{
   higgs_idx = idx;
}

void CLASSNAME::set_loop_order(int loop_order_)
{
   loop_order = loop_order_;
}

void CLASSNAME::set_scale(const Scale_getter& scale_getter_)
{
   scale_getter = scale_getter_;
}

void CLASSNAME::set_scale(Scale_getter&& scale_getter_)
{
   scale_getter = std::move(scale_getter_);
}

void CLASSNAME::set_scale(double scale_)
{
   scale = scale_;
}

void CLASSNAME::set_boundary_scale(const Scale_getter& scale_getter_)
{
   boundary_scale_getter = scale_getter_;
}

void CLASSNAME::set_boundary_scale(Scale_getter&& scale_getter_)
{
   boundary_scale_getter = std::move(scale_getter_);
}

#ifdef CLASSNAME
#undef CLASSNAME
#endif

#define CLASSNAME CNMSSMEFTHiggs_standard_model_matching_down<Semi_analytic>

CLASSNAME::CNMSSMEFTHiggs_standard_model_matching_down(
   standard_model::StandardModel<Two_scale>* low_,
   CNMSSMEFTHiggs<Semi_analytic>* high_,
   const Scale_getter& scale_getter_,
   int loop_order_,
   int higgs_idx_)
   : model(high_)
   , eft(low_)
   , scale_getter(scale_getter_)
   , scale(0.)
   , higgs_idx(higgs_idx_)
{
   set_loop_order(loop_order_);
}

double CLASSNAME::get_scale() const
{
   if (scale != 0.)
      return scale;

   if (!scale_getter)
      throw SetupError("Scale getter is not set!");

   return scale_getter();
}

void CLASSNAME::set_models(Model* high, Model* low)
{
   eft = cast_model<standard_model::StandardModel<Two_scale>*>(low);
   model = cast_model<CNMSSMEFTHiggs<Semi_analytic>*>(high);
}

void CLASSNAME::match()
{
   if (!model || !eft)
      throw SetupError("Model pointer in matching class is NULL!");

   eft->run_to(get_scale());
   model->run_to(get_scale());

   model->calculate_semi_analytic_solutions(get_boundary_scale());

   if (model->get_thresholds() && loop_order)
      CNMSSMEFTHiggs_standard_model_matching::match_high_to_low_scale_model(*eft, *model, loop_order, higgs_idx);
   else
      CNMSSMEFTHiggs_standard_model_matching::match_high_to_low_scale_model_tree_level(*eft, *model, higgs_idx);
}

void CLASSNAME::match_tree_level()
{
   if (!model || !eft)
      throw SetupError("Model pointer in matching class is NULL!");

   eft->run_to(get_scale());
   model->run_to(get_scale());

   model->calculate_semi_analytic_solutions(get_boundary_scale());

   CNMSSMEFTHiggs_standard_model_matching::match_high_to_low_scale_model_tree_level(*eft, *model, higgs_idx);
}

int CLASSNAME::get_higgs_index() const
{
   return higgs_idx;
}

int CLASSNAME::get_loop_order() const
{
   return loop_order;
}

double CLASSNAME::get_boundary_scale() const
{
   if (!boundary_scale_getter)
      return get_scale();
   return boundary_scale_getter();
}

void CLASSNAME::set_higgs_index(int idx)
{
   higgs_idx = idx;
}

void CLASSNAME::set_loop_order(int loop_order_)
{
   if (loop_order_ > 1) {
      WARNING("Matching loop order " << loop_order_
              << " for downwards matching currently not"
              " supported!  I'm using 1-loop matching.");
      loop_order_ = 1;
   }

   loop_order = loop_order_;
}

void CLASSNAME::set_scale(const Scale_getter& scale_getter_)
{
   scale_getter = scale_getter_;
}

void CLASSNAME::set_scale(Scale_getter&& scale_getter_)
{
   scale_getter = std::move(scale_getter_);
}

void CLASSNAME::set_scale(double scale_)
{
   scale = scale_;
}

void CLASSNAME::set_boundary_scale(const Scale_getter& scale_getter_)
{
   boundary_scale_getter = scale_getter_;
}

void CLASSNAME::set_boundary_scale(Scale_getter&& scale_getter_)
{
   boundary_scale_getter = std::move(scale_getter_);
}

} // namespace flexiblesusy