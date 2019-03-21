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

// File generated at Thu 21 Mar 2019 19:41:49

/**
 * @file CNMSSM_semi_analytic_model.hpp
 * @brief contains class for model with routines needed to solve boundary
 *        value problem using the semi_analytic solver by solving EWSB
 *        and determine the pole masses and mixings
 *
 * This file was generated at Thu 21 Mar 2019 19:41:49 with FlexibleSUSY
 * 2.3.0 (git commit: 532ab8d2088eac7a39c987f9276436f985afed75) and SARAH 4.12.3 .
 */

#ifndef CNMSSM_SEMI_ANALYTIC_MODEL_H
#define CNMSSM_SEMI_ANALYTIC_MODEL_H

#include "CNMSSM_model.hpp"
#include "CNMSSM_mass_eigenstates.hpp"
#include "CNMSSM_semi_analytic_solutions.hpp"

#include "model.hpp"

namespace flexiblesusy {

class Semi_analytic;
/**
 * @class CNMSSM<Semi_analytic>
 * @brief model class with routines for determining masses and mixings and EWSB
 */
template<>
class CNMSSM<Semi_analytic> : public Model, public CNMSSM_mass_eigenstates {
public:
   explicit CNMSSM(const CNMSSM_input_parameters& input_ = CNMSSM_input_parameters());
   CNMSSM(const CNMSSM&) = default;
   CNMSSM(CNMSSM&&) = default;
   virtual ~CNMSSM() = default;
   CNMSSM& operator=(const CNMSSM&) = default;
   CNMSSM& operator=(CNMSSM&&) = default;

   // interface functions
   virtual void calculate_spectrum() override;
   virtual void clear_problems() override;
   virtual std::string name() const override;
   virtual void run_to(double scale, double eps = -1.0) override;
   virtual void print(std::ostream& out = std::cerr) const override;
   virtual void set_precision(double) override;

   /**
    * @brief returns the current values of the semi-analytic coefficients
    *
    * @return the current set of semi-analytic solutions
    */
   const CNMSSM_semi_analytic_solutions& get_semi_analytic_solutions() const;

   /**
    * @brief returns the current values of the semi-analytic coefficients
    *
    * @return the current set of semi-analytic solutions
    */
   CNMSSM_semi_analytic_solutions& get_semi_analytic_solutions();

   /**
    * @brief calculates the semi-analytic solutions for the soft parameters
    *
    * @param[in] input_scale the scale where the boundary conditions hold
    */
   void calculate_semi_analytic_solutions(double);

   double get_MassBCoeff1() const { return solutions.get_MassBCoeff1(); }
   double get_MassBCoeff2() const { return solutions.get_MassBCoeff2(); }

   double get_MassGCoeff1() const { return solutions.get_MassGCoeff1(); }
   double get_MassGCoeff2() const { return solutions.get_MassGCoeff2(); }

   double get_MassWBCoeff1() const { return solutions.get_MassWBCoeff1(); }
   double get_MassWBCoeff2() const { return solutions.get_MassWBCoeff2(); }

   const Eigen::Matrix<double,3,3>& get_TYdCoeff1() const { return solutions.
      get_TYdCoeff1(); }
   const Eigen::Matrix<double,3,3>& get_TYdCoeff2() const { return solutions.
      get_TYdCoeff2(); }

   const Eigen::Matrix<double,3,3>& get_TYeCoeff1() const { return solutions.
      get_TYeCoeff1(); }
   const Eigen::Matrix<double,3,3>& get_TYeCoeff2() const { return solutions.
      get_TYeCoeff2(); }

   const Eigen::Matrix<double,3,3>& get_TYuCoeff1() const { return solutions.
      get_TYuCoeff1(); }
   const Eigen::Matrix<double,3,3>& get_TYuCoeff2() const { return solutions.
      get_TYuCoeff2(); }

   double get_TKappaCoeff1() const { return solutions.get_TKappaCoeff1(); }
   double get_TKappaCoeff2() const { return solutions.get_TKappaCoeff2(); }

   double get_TLambdaxCoeff1() const { return solutions.get_TLambdaxCoeff1(); }
   double get_TLambdaxCoeff2() const { return solutions.get_TLambdaxCoeff2(); }

   const Eigen::Matrix<double,3,3>& get_md2Coeff1() const { return solutions.
      get_md2Coeff1(); }
   const Eigen::Matrix<double,3,3>& get_md2Coeff2() const { return solutions.
      get_md2Coeff2(); }
   const Eigen::Matrix<double,3,3>& get_md2Coeff3() const { return solutions.
      get_md2Coeff3(); }
   const Eigen::Matrix<double,3,3>& get_md2Coeff4() const { return solutions.
      get_md2Coeff4(); }

   const Eigen::Matrix<double,3,3>& get_me2Coeff1() const { return solutions.
      get_me2Coeff1(); }
   const Eigen::Matrix<double,3,3>& get_me2Coeff2() const { return solutions.
      get_me2Coeff2(); }
   const Eigen::Matrix<double,3,3>& get_me2Coeff3() const { return solutions.
      get_me2Coeff3(); }
   const Eigen::Matrix<double,3,3>& get_me2Coeff4() const { return solutions.
      get_me2Coeff4(); }

   double get_mHd2Coeff1() const { return solutions.get_mHd2Coeff1(); }
   double get_mHd2Coeff2() const { return solutions.get_mHd2Coeff2(); }
   double get_mHd2Coeff3() const { return solutions.get_mHd2Coeff3(); }
   double get_mHd2Coeff4() const { return solutions.get_mHd2Coeff4(); }

   double get_mHu2Coeff1() const { return solutions.get_mHu2Coeff1(); }
   double get_mHu2Coeff2() const { return solutions.get_mHu2Coeff2(); }
   double get_mHu2Coeff3() const { return solutions.get_mHu2Coeff3(); }
   double get_mHu2Coeff4() const { return solutions.get_mHu2Coeff4(); }

   const Eigen::Matrix<double,3,3>& get_ml2Coeff1() const { return solutions.
      get_ml2Coeff1(); }
   const Eigen::Matrix<double,3,3>& get_ml2Coeff2() const { return solutions.
      get_ml2Coeff2(); }
   const Eigen::Matrix<double,3,3>& get_ml2Coeff3() const { return solutions.
      get_ml2Coeff3(); }
   const Eigen::Matrix<double,3,3>& get_ml2Coeff4() const { return solutions.
      get_ml2Coeff4(); }

   const Eigen::Matrix<double,3,3>& get_mq2Coeff1() const { return solutions.
      get_mq2Coeff1(); }
   const Eigen::Matrix<double,3,3>& get_mq2Coeff2() const { return solutions.
      get_mq2Coeff2(); }
   const Eigen::Matrix<double,3,3>& get_mq2Coeff3() const { return solutions.
      get_mq2Coeff3(); }
   const Eigen::Matrix<double,3,3>& get_mq2Coeff4() const { return solutions.
      get_mq2Coeff4(); }

   double get_ms2Coeff1() const { return solutions.get_ms2Coeff1(); }
   double get_ms2Coeff2() const { return solutions.get_ms2Coeff2(); }
   double get_ms2Coeff3() const { return solutions.get_ms2Coeff3(); }
   double get_ms2Coeff4() const { return solutions.get_ms2Coeff4(); }

   const Eigen::Matrix<double,3,3>& get_mu2Coeff1() const { return solutions.
      get_mu2Coeff1(); }
   const Eigen::Matrix<double,3,3>& get_mu2Coeff2() const { return solutions.
      get_mu2Coeff2(); }
   const Eigen::Matrix<double,3,3>& get_mu2Coeff3() const { return solutions.
      get_mu2Coeff3(); }
   const Eigen::Matrix<double,3,3>& get_mu2Coeff4() const { return solutions.
      get_mu2Coeff4(); }

private:
   /// semi-analytic solutions for the model soft parameters
   CNMSSM_semi_analytic_solutions solutions{};
};

std::ostream& operator<<(std::ostream&, const CNMSSM<Semi_analytic>&);

} // namespace flexiblesusy

#endif
