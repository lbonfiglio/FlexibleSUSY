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

// File generated at Thu 21 Mar 2019 14:13:30

#include "CNMSSMEFTHiggs_input_parameters.hpp"
#include "wrappers.hpp"

#define INPUT(p) input.p

namespace flexiblesusy {

Eigen::ArrayXd CNMSSMEFTHiggs_input_parameters::get() const
{
   Eigen::ArrayXd pars(5);

   pars(0) = m12;
   pars(1) = TanBeta;
   pars(2) = SignvS;
   pars(3) = Azero;
   pars(4) = LambdaInput;

   return pars;
}

void CNMSSMEFTHiggs_input_parameters::set(const Eigen::ArrayXd& pars)
{
   m12 = pars(0);
   TanBeta = pars(1);
   SignvS = pars(2);
   Azero = pars(3);
   LambdaInput = pars(4);

}

std::ostream& operator<<(std::ostream& ostr, const CNMSSMEFTHiggs_input_parameters& input)
{
   ostr << "m12 = " << INPUT(m12) << ", ";
   ostr << "TanBeta = " << INPUT(TanBeta) << ", ";
   ostr << "SignvS = " << INPUT(SignvS) << ", ";
   ostr << "Azero = " << INPUT(Azero) << ", ";
   ostr << "LambdaInput = " << INPUT(LambdaInput) << ", ";

   return ostr;
}

} // namespace flexiblesusy
