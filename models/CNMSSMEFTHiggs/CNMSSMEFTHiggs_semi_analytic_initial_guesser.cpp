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

// File generated at Thu 21 Mar 2019 14:21:43

#include "CNMSSMEFTHiggs_semi_analytic_initial_guesser.hpp"
#include "CNMSSMEFTHiggs_semi_analytic_model.hpp"
#include "CNMSSMEFTHiggs_semi_analytic_susy_convergence_tester.hpp"
#include "standard_model_two_scale_model.hpp"
#include "lowe.h"
#include "error.hpp"
#include "ew_input.hpp"
#include "two_scale_running_precision.hpp"
#include "two_scale_solver.hpp"
#include "wrappers.hpp"

#include <Eigen/Core>

namespace flexiblesusy {

#define INPUTPARAMETER(p) model->get_input().p
#define MODELPARAMETER(p) model->get_##p()
#define SMPARAMETER(p) eft->get_##p()
#define PHASE(p) model->get_##p()
#define LowEnergyConstant(p) Electroweak_constants::p
#define MODEL model

CNMSSMEFTHiggs_standard_model_initial_guesser<Semi_analytic>::CNMSSMEFTHiggs_standard_model_initial_guesser(
   CNMSSMEFTHiggs<Semi_analytic>* model_,
   standard_model::StandardModel<Two_scale>* eft_,
   const softsusy::QedQcd& qedqcd_,
   const standard_model::Standard_model_low_scale_constraint<Two_scale>& low_constraint_,
   CNMSSMEFTHiggs_susy_scale_constraint<Semi_analytic>& susy_constraint_,
   CNMSSMEFTHiggs_high_scale_constraint<Semi_analytic>& high_constraint_,
   CNMSSMEFTHiggs_matching_constraint<Semi_analytic>& matching_constraint_,
   CNMSSMEFTHiggs_standard_model_matching_up<Semi_analytic>& matching_up_
)
   : model(model_)
   , eft(eft_)
   , qedqcd(qedqcd_)
   , low_constraint(low_constraint_)
   , susy_constraint(susy_constraint_)
   , high_constraint(high_constraint_)
   , matching_constraint(matching_constraint_)
   , matching_up(matching_up_)
{
   if (!model)
      throw SetupError("CNMSSMEFTHiggs_initial_guesser: Error: pointer to model"
                       " CNMSSMEFTHiggs<Semi_analytic> must not be zero");
}

/**
 * Guesses the DR-bar model parameters by calling
 * guess_susy_parameters() and guess_soft_parameters() .
 */
void CNMSSMEFTHiggs_standard_model_initial_guesser<Semi_analytic>::guess()
{
   guess_eft_parameters();
   guess_model_parameters();
}

/**
 * Guesses the SUSY parameters (gauge, Yukawa couplings) at
 * \f$m_\text{top}^\text{pole}\f$ from the Standard Model gauge
 * couplings and fermion masses.  Threshold corrections are ignored.
 * The user-defined initial guess at the low-scale
 * (InitialGuessAtLowScale) is applied here:
 *
 * \code{.cpp}
   

 * \endcode
 */
void CNMSSMEFTHiggs_standard_model_initial_guesser<Semi_analytic>::guess_eft_parameters()
{
   softsusy::QedQcd leAtMt(qedqcd);
   const double mtpole = leAtMt.displayPoleMt();

   mu_guess = leAtMt.displayMass(softsusy::mUp);
   mc_guess = leAtMt.displayMass(softsusy::mCharm);
   mt_guess = model->get_thresholds() > 0 && model->get_threshold_corrections().mt > 0 ?
      leAtMt.displayMass(softsusy::mTop) - 30.0 :
      leAtMt.displayPoleMt();
   md_guess = leAtMt.displayMass(softsusy::mDown);
   ms_guess = leAtMt.displayMass(softsusy::mStrange);
   mb_guess = leAtMt.displayMass(softsusy::mBottom);
   me_guess = model->get_thresholds() > 0 ?
      leAtMt.displayMass(softsusy::mElectron) :
      leAtMt.displayPoleMel();
   mm_guess = model->get_thresholds() > 0 ?
      leAtMt.displayMass(softsusy::mMuon) :
      leAtMt.displayPoleMmuon();
   mtau_guess = leAtMt.displayMass(softsusy::mTau);

   // guess gauge couplings at mt
   const auto alpha_sm(leAtMt.guess_alpha_SM5(mtpole));

   eft->set_g1(sqrt(4.0 * M_PI * alpha_sm(0)));
   eft->set_g2(sqrt(4.0 * M_PI * alpha_sm(1)));
   eft->set_g3(sqrt(4.0 * M_PI * alpha_sm(2)));
   eft->set_scale(mtpole);

   eft->set_v(Electroweak_constants::vev);
   eft->set_Yu(ZEROMATRIX(3,3));
   eft->set_Yd(ZEROMATRIX(3,3));
   eft->set_Ye(ZEROMATRIX(3,3));

   eft->set_Yu(0, 0, -Sqrt(2.)* mu_guess/ eft->get_v());
   eft->set_Yu(1, 1, -Sqrt(2.)* mc_guess/ eft->get_v());
   eft->set_Yu(2, 2, -Sqrt(2.)* mt_guess/ eft->get_v());

   eft->set_Yd(0, 0, Sqrt(2.)* md_guess/ eft->get_v());
   eft->set_Yd(1, 1, Sqrt(2.)* ms_guess/ eft->get_v());
   eft->set_Yd(2, 2, Sqrt(2.)* mb_guess/ eft->get_v());

   eft->set_Ye(0, 0, Sqrt(2.)* me_guess/ eft->get_v());
   eft->set_Ye(1, 1, Sqrt(2.)* mm_guess/ eft->get_v());
   eft->set_Ye(2, 2, Sqrt(2.)* mtau_guess/ eft->get_v());

   eft->set_Lambdax(0.12604);
   eft->solve_ewsb_tree_level();
}

void CNMSSMEFTHiggs_standard_model_initial_guesser<Semi_analytic>::calculate_DRbar_yukawa_couplings()
{
   calculate_Yu_DRbar();
   calculate_Yd_DRbar();
   calculate_Ye_DRbar();
}

/**
 * Calculates the Yukawa couplings Yu of the up-type quarks
 * from the Standard Model up-type quark masses (ignoring threshold
 * corrections).
 */
void CNMSSMEFTHiggs_standard_model_initial_guesser<Semi_analytic>::calculate_Yu_DRbar()
{
   Eigen::Matrix<std::complex<double>,3,3> upQuarksDRbar(ZEROMATRIXCOMPLEX(3,3));
   upQuarksDRbar(0,0) = mu_guess;
   upQuarksDRbar(1,1) = mc_guess;
   upQuarksDRbar(2,2) = mt_guess;


}

/**
 * Calculates the Yukawa couplings Yd of the down-type
 * quarks from the Standard Model down-type quark masses (ignoring
 * threshold corrections).
 */
void CNMSSMEFTHiggs_standard_model_initial_guesser<Semi_analytic>::calculate_Yd_DRbar()
{
   Eigen::Matrix<std::complex<double>,3,3> downQuarksDRbar(ZEROMATRIXCOMPLEX(3,3));
   downQuarksDRbar(0,0) = md_guess;
   downQuarksDRbar(1,1) = ms_guess;
   downQuarksDRbar(2,2) = mb_guess;


}

/**
 * Calculates the Yukawa couplings Ye of the leptons
 * from the Standard Model down-type lepton masses (ignoring threshold
 * corrections).
 */
void CNMSSMEFTHiggs_standard_model_initial_guesser<Semi_analytic>::calculate_Ye_DRbar()
{
   Eigen::Matrix<std::complex<double>,3,3> downLeptonsDRbar(ZEROMATRIXCOMPLEX(3,3));
   downLeptonsDRbar(0,0) = me_guess;
   downLeptonsDRbar(1,1) = mm_guess;
   downLeptonsDRbar(2,2) = mtau_guess;


}

void CNMSSMEFTHiggs_standard_model_initial_guesser<Semi_analytic>::guess_model_parameters()
{
   initial_guess_susy_scale_parameters();
   initial_guess_high_scale_parameters();
   solve_susy_parameters();
   guess_soft_parameters();
}

void CNMSSMEFTHiggs_standard_model_initial_guesser<Semi_analytic>::initial_guess_susy_scale_parameters()
{
   const double susy_scale_guess = susy_constraint.get_initial_scale_guess();

   model->set_scale(susy_scale_guess);

   // apply susy-scale first guess
   {
   const auto TanBeta = INPUTPARAMETER(TanBeta);

   MODEL->set_vu(Re((TanBeta*LowEnergyConstant(vev))/Sqrt(1 + Sqr(TanBeta))));
   MODEL->set_vd(Re(LowEnergyConstant(vev)/Sqrt(1 + Sqr(TanBeta))));

   }

   eft->run_to(susy_scale_guess, running_precision);
   eft->calculate_DRbar_masses();

   // get gauge and Yukawa couplings from effective theory
   matching_up.match_tree_level();

   model->run_to(susy_scale_guess, running_precision);

   // apply susy-scale constraint
   susy_constraint.apply();
}

void CNMSSMEFTHiggs_standard_model_initial_guesser<Semi_analytic>::initial_guess_high_scale_parameters()
{
   const double high_scale_guess = high_constraint.get_initial_scale_guess();

   // run to high scale
   model->run_to(high_scale_guess, running_precision);

   // apply user-defined initial guess at the high scale
   {
   const auto Azero = INPUTPARAMETER(Azero);
   const auto m12 = INPUTPARAMETER(m12);
   const auto LambdaInput = INPUTPARAMETER(LambdaInput);
   const auto Ye = MODELPARAMETER(Ye);
   const auto Yd = MODELPARAMETER(Yd);
   const auto Yu = MODELPARAMETER(Yu);
   const auto Kappa = MODELPARAMETER(Kappa);

   MODEL->set_TYe((Azero*Ye).real());
   MODEL->set_TYd((Azero*Yd).real());
   MODEL->set_TYu((Azero*Yu).real());
   MODEL->set_mq2((Sqr(m12)*UNITMATRIX(3)).real());
   MODEL->set_ml2((Sqr(m12)*UNITMATRIX(3)).real());
   MODEL->set_md2((Sqr(m12)*UNITMATRIX(3)).real());
   MODEL->set_mu2((Sqr(m12)*UNITMATRIX(3)).real());
   MODEL->set_me2((Sqr(m12)*UNITMATRIX(3)).real());
   MODEL->set_mHu2(Re(Sqr(m12)));
   MODEL->set_mHd2(Re(Sqr(m12)));
   MODEL->set_ms2(Re(Sqr(m12)));
   MODEL->set_Lambdax(Re(LambdaInput));
   MODEL->set_TKappa(Re(Azero*Kappa));
   MODEL->set_TLambdax(Re(Azero*LambdaInput));
   MODEL->set_MassB(Re(m12));
   MODEL->set_MassWB(Re(m12));
   MODEL->set_MassG(Re(m12));

   }

   // apply high-scale constraint
   high_constraint.apply();
}

/**
 * Performs an initial iteration for the SUSY parameters, ignoring
 * threshold corrections.
 */
void CNMSSMEFTHiggs_standard_model_initial_guesser<Semi_analytic>::solve_susy_parameters()
{
   const double susy_scale_guess = susy_constraint.get_initial_scale_guess();
   model->run_to(susy_scale_guess, running_precision);

   CNMSSMEFTHiggs_susy_convergence_tester<Semi_analytic> convergence_tester(
      model, running_precision);
   Two_scale_increasing_precision precision(10.0, running_precision);

   RGFlow<Two_scale> initial_solver;
   initial_solver.set_convergence_tester(&convergence_tester);
   initial_solver.set_running_precision(&precision);
   initial_solver.add(&matching_constraint, model);
   initial_solver.add(&high_constraint, model);
   initial_solver.solve();
}

void CNMSSMEFTHiggs_standard_model_initial_guesser<Semi_analytic>::guess_soft_parameters()
{
   const double susy_scale_guess = matching_constraint.get_scale();
   const double input_scale = high_constraint.get_scale();

   model->run_to(susy_scale_guess, running_precision);

   // apply EWSB constraint
   model->calculate_semi_analytic_solutions(input_scale);
   model->solve_ewsb_tree_level();

   // calculate tree-level spectrum
   model->calculate_DRbar_masses();
}

} // namespace flexiblesusy
