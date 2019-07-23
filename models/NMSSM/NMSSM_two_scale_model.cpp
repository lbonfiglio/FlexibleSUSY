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

// File generated at Tue 23 Jul 2019 15:16:39

/**
 * @file NMSSM_two_scale_model.cpp
 * @brief implementation of the NMSSM model class
 *
 * Contains the definition of the NMSSM model class methods
 * which solve EWSB and calculate pole masses and mixings from DRbar
 * parameters.
 *
 * This file was generated at Tue 23 Jul 2019 15:16:39 with FlexibleSUSY
 * 2.3.0 (git commit: b132b3936f5e72e9f34d3f552dab7c81fb2c369e) and SARAH 4.12.3 .
 */

#include "NMSSM_two_scale_model.hpp"

namespace flexiblesusy {

#define CLASSNAME NMSSM<Two_scale>

CLASSNAME::NMSSM(const NMSSM_input_parameters& input_)
   : NMSSM_mass_eigenstates(input_)
{
}

void CLASSNAME::calculate_spectrum()
{
   NMSSM_mass_eigenstates::calculate_spectrum();
}

void CLASSNAME::clear_problems()
{
   NMSSM_mass_eigenstates::clear_problems();
}

std::string CLASSNAME::name() const
{
   return NMSSM_mass_eigenstates::name();
}

void CLASSNAME::run_to(double scale, double eps)
{
   NMSSM_mass_eigenstates::run_to(scale, eps);
}

void CLASSNAME::print(std::ostream& out) const
{
   NMSSM_mass_eigenstates::print(out);
}

void CLASSNAME::set_precision(double p)
{
   NMSSM_mass_eigenstates::set_precision(p);
}

std::ostream& operator<<(std::ostream& ostr, const NMSSM<Two_scale>& model)
{
   model.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
