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

#ifndef CNMSSMEFTHiggs_SEMI_ANALYTIC_MATCHING_CONSTRAINT_H
#define CNMSSMEFTHiggs_SEMI_ANALYTIC_MATCHING_CONSTRAINT_H

#include "CNMSSMEFTHiggs_matching_constraint.hpp"
#include "single_scale_constraint.hpp"

namespace flexiblesusy {

template <class T>
class CNMSSMEFTHiggs;

template <class T>
class CNMSSMEFTHiggs_standard_model_matching_up;

class Semi_analytic;

template<>
class CNMSSMEFTHiggs_matching_constraint<Semi_analytic> : public Single_scale_constraint {
public:
   CNMSSMEFTHiggs_matching_constraint() = default;
   CNMSSMEFTHiggs_matching_constraint(CNMSSMEFTHiggs<Semi_analytic>*);
   virtual ~CNMSSMEFTHiggs_matching_constraint() = default;
   virtual void apply() override;
   virtual double get_scale() const override { return scale; }
   virtual std::string name() const override { return "CNMSSMEFTHiggs matching constraint"; }
   virtual void set_model(Model*) override;

   void set_matching_condition(CNMSSMEFTHiggs_standard_model_matching_up<Semi_analytic>* mc) { matching_condition = mc; }
   void set_scale(double s) { scale = s; }
private:
   double scale{0.};
   double initial_scale_guess{0.};
   CNMSSMEFTHiggs<Semi_analytic>* model{nullptr};
   CNMSSMEFTHiggs_standard_model_matching_up<Semi_analytic>* matching_condition{nullptr};

   void check_model_ptr() const;
   void check_matching_condition_ptr() const;
};

} // namespace flexiblesusy

#endif
