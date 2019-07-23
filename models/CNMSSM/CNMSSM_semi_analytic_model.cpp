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

// File generated at Tue 23 Jul 2019 15:26:08

/**
 * @file CNMSSM_semi_analytic_model.cpp
 * @brief implementation of the CNMSSM model class
 *
 * Contains the definition of the CNMSSM model class methods
 * which solve EWSB and calculate pole masses and mixings from DRbar
 * parameters.
 *
 * This file was generated at Tue 23 Jul 2019 15:26:08 with FlexibleSUSY
 * 2.3.0 (git commit: b132b3936f5e72e9f34d3f552dab7c81fb2c369e) and SARAH 4.12.3 .
 */

#include "CNMSSM_semi_analytic_model.hpp"

namespace flexiblesusy {

using namespace CNMSSM_info;

#define CLASSNAME CNMSSM<Semi_analytic>

#define COEFFICIENT(coefficient) solutions.get_##coefficient()
#define INPUTPARAMETER(parameter) model.get_input().parameter
#define MODELPARAMETER(parameter) model.get_##parameter()

CLASSNAME::CNMSSM(const CNMSSM_input_parameters& input_)
   : CNMSSM_mass_eigenstates(input_)
{
}

void CLASSNAME::calculate_spectrum()
{
   CNMSSM_mass_eigenstates::calculate_spectrum();
}

void CLASSNAME::clear_problems()
{
   CNMSSM_mass_eigenstates::clear_problems();
}

std::string CLASSNAME::name() const
{
   return CNMSSM_mass_eigenstates::name();
}

void CLASSNAME::run_to(double scale, double eps)
{
   CNMSSM_mass_eigenstates::run_to(scale, eps);
}

void CLASSNAME::print(std::ostream& out) const
{
   CNMSSM_mass_eigenstates::print(out);
   out << "----------------------------------------\n"
          "semi-analytic coefficients:\n"
          "----------------------------------------\n";
   out << "input_scale = " << solutions.get_input_scale() << '\n';
   out << "output_scale = " << solutions.get_output_scale() << '\n';
   out << "MassBCoeff1 = " << COEFFICIENT(MassBCoeff1) << '\n';
   out << "MassBCoeff2 = " << COEFFICIENT(MassBCoeff2) << '\n';
   out << "MassGCoeff1 = " << COEFFICIENT(MassGCoeff1) << '\n';
   out << "MassGCoeff2 = " << COEFFICIENT(MassGCoeff2) << '\n';
   out << "MassWBCoeff1 = " << COEFFICIENT(MassWBCoeff1) << '\n';
   out << "MassWBCoeff2 = " << COEFFICIENT(MassWBCoeff2) << '\n';
   out << "TYdCoeff1 = " << COEFFICIENT(TYdCoeff1) << '\n';
   out << "TYdCoeff2 = " << COEFFICIENT(TYdCoeff2) << '\n';
   out << "TYeCoeff1 = " << COEFFICIENT(TYeCoeff1) << '\n';
   out << "TYeCoeff2 = " << COEFFICIENT(TYeCoeff2) << '\n';
   out << "TYuCoeff1 = " << COEFFICIENT(TYuCoeff1) << '\n';
   out << "TYuCoeff2 = " << COEFFICIENT(TYuCoeff2) << '\n';
   out << "TKappaCoeff1 = " << COEFFICIENT(TKappaCoeff1) << '\n';
   out << "TKappaCoeff2 = " << COEFFICIENT(TKappaCoeff2) << '\n';
   out << "TLambdaxCoeff1 = " << COEFFICIENT(TLambdaxCoeff1) << '\n';
   out << "TLambdaxCoeff2 = " << COEFFICIENT(TLambdaxCoeff2) << '\n';
   out << "md2Coeff1 = " << COEFFICIENT(md2Coeff1) << '\n';
   out << "md2Coeff2 = " << COEFFICIENT(md2Coeff2) << '\n';
   out << "md2Coeff3 = " << COEFFICIENT(md2Coeff3) << '\n';
   out << "md2Coeff4 = " << COEFFICIENT(md2Coeff4) << '\n';
   out << "me2Coeff1 = " << COEFFICIENT(me2Coeff1) << '\n';
   out << "me2Coeff2 = " << COEFFICIENT(me2Coeff2) << '\n';
   out << "me2Coeff3 = " << COEFFICIENT(me2Coeff3) << '\n';
   out << "me2Coeff4 = " << COEFFICIENT(me2Coeff4) << '\n';
   out << "mHd2Coeff1 = " << COEFFICIENT(mHd2Coeff1) << '\n';
   out << "mHd2Coeff2 = " << COEFFICIENT(mHd2Coeff2) << '\n';
   out << "mHd2Coeff3 = " << COEFFICIENT(mHd2Coeff3) << '\n';
   out << "mHd2Coeff4 = " << COEFFICIENT(mHd2Coeff4) << '\n';
   out << "mHu2Coeff1 = " << COEFFICIENT(mHu2Coeff1) << '\n';
   out << "mHu2Coeff2 = " << COEFFICIENT(mHu2Coeff2) << '\n';
   out << "mHu2Coeff3 = " << COEFFICIENT(mHu2Coeff3) << '\n';
   out << "mHu2Coeff4 = " << COEFFICIENT(mHu2Coeff4) << '\n';
   out << "ml2Coeff1 = " << COEFFICIENT(ml2Coeff1) << '\n';
   out << "ml2Coeff2 = " << COEFFICIENT(ml2Coeff2) << '\n';
   out << "ml2Coeff3 = " << COEFFICIENT(ml2Coeff3) << '\n';
   out << "ml2Coeff4 = " << COEFFICIENT(ml2Coeff4) << '\n';
   out << "mq2Coeff1 = " << COEFFICIENT(mq2Coeff1) << '\n';
   out << "mq2Coeff2 = " << COEFFICIENT(mq2Coeff2) << '\n';
   out << "mq2Coeff3 = " << COEFFICIENT(mq2Coeff3) << '\n';
   out << "mq2Coeff4 = " << COEFFICIENT(mq2Coeff4) << '\n';
   out << "ms2Coeff1 = " << COEFFICIENT(ms2Coeff1) << '\n';
   out << "ms2Coeff2 = " << COEFFICIENT(ms2Coeff2) << '\n';
   out << "ms2Coeff3 = " << COEFFICIENT(ms2Coeff3) << '\n';
   out << "ms2Coeff4 = " << COEFFICIENT(ms2Coeff4) << '\n';
   out << "mu2Coeff1 = " << COEFFICIENT(mu2Coeff1) << '\n';
   out << "mu2Coeff2 = " << COEFFICIENT(mu2Coeff2) << '\n';
   out << "mu2Coeff3 = " << COEFFICIENT(mu2Coeff3) << '\n';
   out << "mu2Coeff4 = " << COEFFICIENT(mu2Coeff4) << '\n';

}

void CLASSNAME::set_precision(double p)
{
   CNMSSM_mass_eigenstates::set_precision(p);
}

const CNMSSM_semi_analytic_solutions& CLASSNAME::get_semi_analytic_solutions() const
{
   return solutions;
}

CNMSSM_semi_analytic_solutions& CLASSNAME::get_semi_analytic_solutions()
{
   return solutions;
}

/**
 * Calculates the semi-analytic coefficients at the current scale,
 * assuming that the boundary conditions are imposed at the scale
 * \c input_scale.  The soft parameters are then set to the values
 * obtained using these coefficients and the current values of the
 * boundary value parameters.
 */
void CLASSNAME::calculate_semi_analytic_solutions(double input_scale)
{
   solutions.set_input_scale(input_scale);
   solutions.set_output_scale(get_scale());
   solutions.calculate_coefficients(*this);
   solutions.evaluate_solutions(*this);
}

std::ostream& operator<<(std::ostream& ostr, const CNMSSM<Semi_analytic>& model)
{
   model.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
