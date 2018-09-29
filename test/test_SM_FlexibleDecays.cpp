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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_SM_FlexibleSUSY

#include <boost/test/unit_test.hpp>

//#include "test_SM.hpp"

//#include "SM_two_scale_model.hpp"
//#include "SM_a_muon.hpp"

//#include "softsusy.h"
#include "SM_decays.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_zero2 )
{
   /*
   SM_input_parameters input;
   input.LambdaIN = 0.25;
   input.Qin = 1000;

   SM<Two_scale> sm;
   setup_SM_const(sm, input);
   sm.calculate_DRbar_masses();

   //const double amu = SM_a_muon::calculate_a_muon(sm);

   softsusy::QedQcd qedqcd;
   const auto decays = SM_decays(sm, qedqcd, true);
   */
   BOOST_CHECK_SMALL(1.0, 1e-15);
}
