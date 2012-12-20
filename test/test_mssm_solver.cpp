#include "mssm_solver.h"
#include "mssm_two_scale.hpp"
#include "mssm_two_scale_initial_guesser.hpp"
#include "mssm_two_scale_sugra_constraint.hpp"
#include "mssm_two_scale_msusy_constraint.hpp"
#include "mssm_two_scale_mz_constraint.hpp"
#include "mssm_two_scale_convergence_tester.hpp"
#include "softsusy.h"
#include "two_scale_solver.hpp"
#include "logger.hpp"
#include "stopwatch.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_mssm_solver

#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE( test_softsusy_mssm_solver )
{
   DoubleVector highScaleSoftPars(3);
   const double m12 = 500., a0 = 0., m0 = 125.;
   highScaleSoftPars(1) = m0;
   highScaleSoftPars(2) = m12;
   highScaleSoftPars(3) = a0;

   const int signMu = 1;
   const double tanBeta = 10.0;
   const bool uni = true;
   const double mxGuess = 1.0e16;

   QedQcd oneset;
   const double alphasMZ = 0.1187, mtop = 173.4, mbmb = 4.2;
   oneset.setAlpha(ALPHAS, alphasMZ);
   oneset.setPoleMt(mtop);
   oneset.setMass(mBottom, mbmb);
   oneset.toMz();

   RGFlow<SoftSusy_t> mssmSolver;
   mssmSolver.setSoftHighScalePars(highScaleSoftPars);
   mssmSolver.setSignMu(signMu);
   mssmSolver.setTanBeta(tanBeta);
   mssmSolver.setLowScaleBoundaryConditions(oneset);
   mssmSolver.setGaugeUnification(uni);
   mssmSolver.setMxGuess(mxGuess);
   mssmSolver.setHighScaleBoundaryCondition(sugraBcs);
   mssmSolver.solve();

   MssmSoftsusy softSusy;
   softSusy.lowOrg(sugraBcs, mxGuess, highScaleSoftPars, signMu, tanBeta, oneset, uni);

   BOOST_CHECK_EQUAL(mssmSolver.displayPhys(), softSusy.displayPhys());
}

struct Mssm_parameter_point {
   Mssm_parameter_point()
      : m0(125.0)
      , m12(500.0)
      , a0(0.0)
      , mxGuess(1.0e16)
      , signMu(1)
      , tanBeta(10.0)
      , oneset()
   {
      const double alphasMZ = 0.1187, mtop = 173.4, mbmb = 4.2;
      oneset.setAlpha(ALPHAS, alphasMZ);
      oneset.setPoleMt(mtop);
      oneset.setMass(mBottom, mbmb);
      oneset.toMz();
   }
   DoubleVector get_soft_pars() const {
      DoubleVector highScaleSoftPars(3);
      highScaleSoftPars(1) = m0;
      highScaleSoftPars(2) = m12;
      highScaleSoftPars(3) = a0;
      return highScaleSoftPars;
   }
   friend std::ostream& operator<<(std::ostream& os, const Mssm_parameter_point& pp) {
      os << "CMSSM parameter point:"
         << " m0=" << pp.m0
         << ", m12=" << pp.m12
         << ", a0=" << pp.a0
         << ", mxGuess=" << pp.mxGuess
         << ", signMu=" << pp.signMu
         << ", tanBeta=" << pp.tanBeta
         << ", QedQcd=" << pp.oneset
         << '\n';
      return os;
   }

   double m0, m12, a0, mxGuess, signMu, tanBeta;
   QedQcd oneset;
};

/**
 * Tests if our two scale algorithm calculates the same spectrum as
 * SoftSusy
 *
 * @param pp CMSSM parameter point
 */
void test_point(const Mssm_parameter_point& pp)
{
   // record time for the two scale method to solve the MSSM
   Stopwatch stopwatch;
   stopwatch.start();

   // setup the MSSM with the two scale method
   Mssm<Two_scale> mssm;
   Mssm_sugra_constraint mssm_sugra_constraint(&mssm, pp.mxGuess, pp.m0, pp.m12, pp.a0, pp.signMu);
   Mssm_mz_constraint mssm_mz_constraint(&mssm, pp.tanBeta);
   Mssm_msusy_constraint mssm_msusy_constraint(&mssm, pp.get_soft_pars(), 1000.0, pp.signMu);
   Mssm_convergence_tester mssm_convergence_tester(&mssm, 1.0e-4);
   Mssm_initial_guesser initial_guesser(&mssm, pp.oneset, pp.mxGuess, pp.tanBeta, pp.signMu, pp.get_soft_pars(), false);

   std::vector<Constraint<Two_scale>*> mssm_upward_constraints;
   mssm_upward_constraints.push_back(&mssm_mz_constraint);
   mssm_upward_constraints.push_back(&mssm_sugra_constraint);

   std::vector<Constraint<Two_scale>*> mssm_downward_constraints;
   mssm_downward_constraints.push_back(&mssm_sugra_constraint);
   mssm_downward_constraints.push_back(&mssm_msusy_constraint);
   mssm_downward_constraints.push_back(&mssm_mz_constraint);

   RGFlow<Two_scale> solver;
   solver.set_max_iterations(10);
   solver.set_convergence_tester(&mssm_convergence_tester);
   solver.set_increasing_running_precision(true);
   solver.set_initial_guesser(&initial_guesser);
   solver.add_model(&mssm, mssm_upward_constraints, mssm_downward_constraints);
   try {
      solver.solve();
   } catch (RGFlow<Two_scale>::Error& e) {
      BOOST_ERROR(e.what());
   }
   mssm.calculate_spectrum();

   stopwatch.stop();
   VERBOSE_MSG("Mssm<Two_scale> solved in " << stopwatch.get_time_in_seconds()
               << " seconds (" << stopwatch.get_clicks() << " clicks)");

   // record time for SoftSusy to solve the MSSM
   stopwatch.start();

   // run softsusy
   softsusy::TOLERANCE = 1.0e-4;
   MssmSoftsusy softSusy;
   const double mxSoftSusy
      = softSusy.lowOrg(sugraBcs, pp.mxGuess, pp.get_soft_pars(), pp.signMu, pp.tanBeta, pp.oneset, true);

   stopwatch.stop();
   VERBOSE_MSG("MssmSoftsusy solved in " << stopwatch.get_time_in_seconds()
               << " seconds (" << stopwatch.get_clicks() << " clicks)");

   // check equality of physical parameters
   const sPhysical softSusyPhys(softSusy.displayPhys()), mssmPhys(mssm.displayPhys());

   BOOST_CHECK_CLOSE(softSusyPhys.mh0         , mssmPhys.mh0         , 0.1);
   BOOST_CHECK_CLOSE(softSusyPhys.mA0         , mssmPhys.mA0         , 0.1);
   BOOST_CHECK_CLOSE(softSusyPhys.mH0         , mssmPhys.mH0         , 0.1);
   BOOST_CHECK_CLOSE(softSusyPhys.mHpm        , mssmPhys.mHpm        , 0.1);
   BOOST_CHECK_CLOSE(softSusyPhys.mGluino     , mssmPhys.mGluino     , 0.1);
   BOOST_CHECK_CLOSE(softSusyPhys.thetaL      , mssmPhys.thetaL      , 0.1);
   BOOST_CHECK_CLOSE(softSusyPhys.thetaR      , mssmPhys.thetaR      , 0.1);
   BOOST_CHECK_CLOSE(softSusyPhys.thetat      , mssmPhys.thetat      , 0.1);
   BOOST_CHECK_CLOSE(softSusyPhys.thetab      , mssmPhys.thetab      , 0.1);
   BOOST_CHECK_CLOSE(softSusyPhys.thetatau    , mssmPhys.thetatau    , 0.1);
   BOOST_CHECK_CLOSE(softSusyPhys.thetaH      , mssmPhys.thetaH      , 0.1);
   BOOST_CHECK_CLOSE(softSusyPhys.t1OV1Ms     , mssmPhys.t1OV1Ms     , 0.1);
   BOOST_CHECK_CLOSE(softSusyPhys.t2OV2Ms     , mssmPhys.t2OV2Ms     , 0.1);
   BOOST_CHECK_CLOSE(softSusyPhys.t1OV1Ms1loop, mssmPhys.t1OV1Ms1loop, 0.1);
   BOOST_CHECK_CLOSE(softSusyPhys.t2OV2Ms1loop, mssmPhys.t2OV2Ms1loop, 0.1);

   // sneutrino masses
   for (int i = softSusyPhys.msnu.displayStart();
        i < softSusyPhys.msnu.displayEnd(); ++i)
      BOOST_CHECK_CLOSE(softSusyPhys.msnu(i), mssmPhys.msnu(i), 0.1);

   // chargino masses
   for (int i = softSusyPhys.mch.displayStart();
        i < softSusyPhys.mch.displayEnd(); ++i)
      BOOST_CHECK_CLOSE(softSusyPhys.mch(i), mssmPhys.mch(i), 0.1);

   // neutralino masses
   for (int i = softSusyPhys.mneut.displayStart();
        i < softSusyPhys.mneut.displayEnd(); ++i)
      BOOST_CHECK_CLOSE(softSusyPhys.mneut(i), mssmPhys.mneut(i), 0.1);

   // neuralino mixing matrix
   for (int i = 1; i < softSusyPhys.mixNeut.displayRows(); ++i)
      for (int k = 1; k < softSusyPhys.mixNeut.displayCols(); ++k)
         BOOST_CHECK_CLOSE(softSusyPhys.mixNeut(i,k), mssmPhys.mixNeut(i,k), 0.1);

   // up squarks
   for (int i = 1; i < softSusyPhys.mu.displayRows(); ++i)
      for (int k = 1; k < softSusyPhys.mu.displayCols(); ++k)
         BOOST_CHECK_CLOSE(softSusyPhys.mu(i,k), mssmPhys.mu(i,k), 0.1);

   // down squarks
   for (int i = 1; i < softSusyPhys.md.displayRows(); ++i)
      for (int k = 1; k < softSusyPhys.md.displayCols(); ++k)
         BOOST_CHECK_CLOSE(softSusyPhys.md(i,k), mssmPhys.md(i,k), 0.1);

   // down sleptons
   for (int i = 1; i < softSusyPhys.me.displayRows(); ++i)
      for (int k = 1; k < softSusyPhys.me.displayCols(); ++k)
         BOOST_CHECK_CLOSE(softSusyPhys.me(i,k), mssmPhys.me(i,k), 0.1);

   BOOST_CHECK_CLOSE(mxSoftSusy, mssm_sugra_constraint.get_scale(), 0.1);
}

BOOST_AUTO_TEST_CASE( test_default_cmssm_parameter_point )
{
   Mssm_parameter_point pp;
   BOOST_MESSAGE("testing parameter point " << pp);
   test_point(pp);
}

BOOST_AUTO_TEST_CASE( test_slow_convergence_point )
{
   // slow convergence point taken from arXiv:1211.3231 Fig. 7
   Mssm_parameter_point pp;
   pp.oneset.setPoleMt(173.5);
   pp.tanBeta = 10.0;
   pp.a0 = 0.0;
   pp.m12 = 337.5;
   pp.m0 = 3400.0;

   BOOST_MESSAGE("testing slow convergent parameter point " << pp);
   test_point(pp);
}
