
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_CMSSM_FlexibleDecays

#include <boost/test/unit_test.hpp>

#include "test_CMSSM.hpp"
#include "CMSSM_two_scale_model.hpp"
#include "CMSSM_decays.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_CMSSM_FlexibleDecays )
{
   CMSSM_input_parameters input;
   input.TanBeta = 10.;
   input.m0 = 125.;
   input.m12 = 200.;
   input.SignMu = 1;
   input.Azero = 0.;

   CMSSM<Two_scale> m;

   setup_CMSSM_const(m, input);
   m.calculate_DRbar_masses();

   m.set_pole_mass_loop_order(1);
   m.do_calculate_sm_pole_masses(true);
   m.solve_ewsb_one_loop();
   m.calculate_pole_masses();

   if (m.get_problems().have_problem()) {
      std::ostringstream ostr;
      m.get_problems().print_problems(ostr);
      BOOST_FAIL(ostr.str());
   }

   softsusy::QedQcd qedqcd;
   CMSSM_decays decays(m, qedqcd, true);
   std::cout << m.get_physical().Mhh << '\n';

   // h -> b bbar
   BOOST_CHECK_CLOSE_FRACTION(decays.partial_width_hh_to_barFdFd(m, 0, 2, 2),
                              2.6118180765322455E-003, 2e-15);
   // h -> tau+ tau-
   BOOST_CHECK_CLOSE_FRACTION(decays.partial_width_hh_to_barFeFe(m, 0, 2, 2),
                              2.6800741077127161E-004, 1e-15);
   // h -> W+ W-
   BOOST_CHECK_CLOSE_FRACTION(decays.partial_width_hh_to_VWmconjVWm(m, 0),
                              8.4705126919250480E-004, 1e-14);
   // h -> Z Z
   BOOST_CHECK_CLOSE_FRACTION(decays.partial_width_hh_to_VZVZ(m, 0),
                              9.4231401208556083E-005, 1e-14);
   // h -> gluon gluon
   // BOOST_CHECK_CLOSE_FRACTION(decays.partial_width_hh_to_VGVG(m), , 1e-15);
   // h -> gamma gamma
   // BOOST_CHECK_CLOSE_FRACTION(decays.partial_width_hh_to_VPVP(m), , 1e-15);
   // h -> Z gamma
   // BOOST_CHECK_CLOSE_FRACTION(decays.partial_width_hh_to_VZVP(m), , 1e-15);
}
