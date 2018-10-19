
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_CMSSMNoFV_tree_level_spectrum

#include <boost/test/unit_test.hpp>

#include "test_CMSSMNoFV.hpp"
#include "CMSSMNoFV_two_scale_model.hpp"
#include "CMSSMNoFV_decays.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_CMSSMNoFV_FlexibleDecays )
{
   CMSSMNoFV_input_parameters input;
   input.TanBeta = 10.;
   input.m0 = 125.;
   input.m12 = 200.;
   input.SignMu = 1;
   input.Azero = 0.;

   CMSSMNoFV<Two_scale> m1;

   setup_CMSSM_const(m1, input);

   BOOST_REQUIRE(gErrors == 0);
   if (gErrors) {
      BOOST_FAIL("spectra are not equal");
      gErrors = 0;
   }

   softsusy::QedQcd qedqcd;
   SM_decays decays(m, qedqcd, true);

   // h -> b bbar
   BOOST_CHECK_CLOSE_FRACTION(decays.partial_width_hh_to_barFdFd(m, 2, 2),
                              2.6118180765322455E-003, 2e-15);
   // h -> tau+ tau-
   BOOST_CHECK_CLOSE_FRACTION(decays.partial_width_hh_to_barFeFe(m, 2, 2),
                              2.6800741077127161E-004, 1e-15);
   // h -> W+ W-
   BOOST_CHECK_CLOSE_FRACTION(decays.partial_width_hh_to_VWpconjVWp(m),
                              8.4705126919250480E-004, 1e-14);
   // h -> Z Z
   BOOST_CHECK_CLOSE_FRACTION(decays.partial_width_hh_to_VZVZ(m),
                              9.4231401208556083E-005, 1e-14);
   // h -> gluon gluon
   // BOOST_CHECK_CLOSE_FRACTION(decays.partial_width_hh_to_VGVG(m), , 1e-15);
   // h -> gamma gamma
   // BOOST_CHECK_CLOSE_FRACTION(decays.partial_width_hh_to_VPVP(m), , 1e-15);
   // h -> Z gamma
   // BOOST_CHECK_CLOSE_FRACTION(decays.partial_width_hh_to_VZVP(m), , 1e-15);
}
