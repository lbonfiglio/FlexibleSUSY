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

// File generated at Thu 21 Mar 2019 14:21:41

#ifndef CNMSSMEFTHiggs_SEMI_ANALYTIC_HIGH_SCALE_CONSTRAINT_H
#define CNMSSMEFTHiggs_SEMI_ANALYTIC_HIGH_SCALE_CONSTRAINT_H

#include "CNMSSMEFTHiggs_high_scale_constraint.hpp"
#include "CNMSSMEFTHiggs_input_parameters.hpp"
#include "single_scale_constraint.hpp"

namespace flexiblesusy {

template <class T>
class CNMSSMEFTHiggs;

class Semi_analytic;

template<>
class CNMSSMEFTHiggs_high_scale_constraint<Semi_analytic> : public Single_scale_constraint {
public:

   CNMSSMEFTHiggs_high_scale_constraint() = default;
   CNMSSMEFTHiggs_high_scale_constraint(CNMSSMEFTHiggs<Semi_analytic>*);
   virtual ~CNMSSMEFTHiggs_high_scale_constraint() = default;
   virtual void apply() override;
   virtual double get_scale() const override;
   virtual std::string name() const override { return "CNMSSMEFTHiggs high-scale constraint"; }
   virtual void set_model(Model*) override;

   void clear();
   double get_initial_scale_guess() const { return initial_scale_guess; }
   const CNMSSMEFTHiggs_input_parameters& get_input_parameters() const;
   CNMSSMEFTHiggs<Semi_analytic>* get_model() const { return model; }
   void initialize();
   void set_scale(double s) { scale = s; }

protected:
   void update_scale();
   bool check_non_perturbative();

private:
   double scale{0.};
   double initial_scale_guess{0.};
   CNMSSMEFTHiggs<Semi_analytic>* model{nullptr};

   void check_model_ptr() const;
};

} // namespace flexiblesusy

#endif
