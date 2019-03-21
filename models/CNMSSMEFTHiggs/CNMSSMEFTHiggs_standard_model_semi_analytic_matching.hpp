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

#ifndef CNMSSMEFTHiggs_SEMI_ANALYTIC_MATCHING_H
#define CNMSSMEFTHiggs_SEMI_ANALYTIC_MATCHING_H

#include "single_scale_matching.hpp"

#include <Eigen/Core>

#include <functional>

namespace flexiblesusy {

class Model;
class Semi_analytic;
class Two_scale;

template <class T> class CNMSSMEFTHiggs;
template <class T> class CNMSSMEFTHiggs_standard_model_matching_up;
template <class T> class CNMSSMEFTHiggs_standard_model_matching_down;
template <class T> class CNMSSMEFTHiggs_matching_constraint;

namespace standard_model {
template <class T> class StandardModel;
}

template<>
class CNMSSMEFTHiggs_standard_model_matching_up<Semi_analytic> : public Single_scale_matching
{
public:
   using Scale_getter = std::function<double()>;

   CNMSSMEFTHiggs_standard_model_matching_up() = default;

   CNMSSMEFTHiggs_standard_model_matching_up(standard_model::StandardModel<Two_scale>*,
                                          CNMSSMEFTHiggs<Semi_analytic>*,
                                          const Scale_getter&,
                                          int, int);

   virtual ~CNMSSMEFTHiggs_standard_model_matching_up() = default;

   virtual double get_scale() const override;
   virtual void set_models(Model*, Model*) override;
   virtual void match() override;

   int get_higgs_index() const;
   int get_loop_order() const;
   double get_boundary_scale() const;
   void set_higgs_index(int);
   void set_loop_order(int);
   void set_scale(double);
   void set_scale(const Scale_getter&);
   void set_scale(Scale_getter&&);
   void set_boundary_scale(const Scale_getter&);
   void set_boundary_scale(Scale_getter&&);
   void match_tree_level();

   double get_g1() const { return g1; }
   double get_g2() const { return g2; }
   double get_g3() const { return g3; }
   double get_vd() const { return vd; }
   double get_vu() const { return vu; }
   const Eigen::Matrix<double,3,3>& get_Yd() const { return Yd; }
   const Eigen::Matrix<double,3,3>& get_Ye() const { return Ye; }
   const Eigen::Matrix<double,3,3>& get_Yu() const { return Yu; }

   void set_matching_constraint(CNMSSMEFTHiggs_matching_constraint<Semi_analytic>* mc) { matching_constraint = mc; }

private:
   CNMSSMEFTHiggs<Semi_analytic>* model{nullptr};
   standard_model::StandardModel<Two_scale>* eft{nullptr};
   Scale_getter scale_getter{};
   Scale_getter boundary_scale_getter{};
   double scale{0.};
   int loop_order{0};
   int higgs_idx{0}; ///< Higgs index
   double g1{};
   double g2{};
   double g3{};
   double vd{};
   double vu{};
   Eigen::Matrix<double,3,3> Yd{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> Ye{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> Yu{Eigen::Matrix<double,3,3>::Zero()};

   CNMSSMEFTHiggs_matching_constraint<Semi_analytic>* matching_constraint{nullptr};

   void save_matched_parameters();
};

template<>
class CNMSSMEFTHiggs_standard_model_matching_down<Semi_analytic> : public Single_scale_matching
{
public:
   using Scale_getter = std::function<double()>;

   CNMSSMEFTHiggs_standard_model_matching_down() = default;

   CNMSSMEFTHiggs_standard_model_matching_down(standard_model::StandardModel<Two_scale>*,
                                            CNMSSMEFTHiggs<Semi_analytic>*,
                                            const Scale_getter&,
                                            int, int);

   virtual ~CNMSSMEFTHiggs_standard_model_matching_down() = default;

   virtual double get_scale() const override;
   virtual void set_models(Model*, Model*) override;
   virtual void match() override;

   int get_higgs_index() const;
   int get_loop_order() const;
   double get_boundary_scale() const;
   void set_higgs_index(int);
   void set_loop_order(int);
   void set_scale(double);
   void set_scale(const Scale_getter&);
   void set_scale(Scale_getter&&);
   void set_boundary_scale(const Scale_getter&);
   void set_boundary_scale(Scale_getter&&);
   void match_tree_level();

private:
   CNMSSMEFTHiggs<Semi_analytic>* model{nullptr};
   standard_model::StandardModel<Two_scale>* eft{nullptr};
   Scale_getter scale_getter{};
   Scale_getter boundary_scale_getter{};
   double scale{0.};
   int loop_order{0};
   int higgs_idx{0}; ///< Higgs index
};

} // namespace flexiblesusy

#endif
