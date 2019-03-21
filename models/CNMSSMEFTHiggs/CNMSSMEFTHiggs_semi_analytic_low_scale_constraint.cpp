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

// File generated at Thu 21 Mar 2019 14:21:42

#include "CNMSSMEFTHiggs_semi_analytic_low_scale_constraint.hpp"
#include "CNMSSMEFTHiggs_semi_analytic_model.hpp"
#include "CNMSSMEFTHiggs_info.hpp"
#include "CNMSSMEFTHiggs_weinberg_angle.hpp"
#include "wrappers.hpp"
#include "logger.hpp"
#include "ew_input.hpp"
#include "gsl_utils.hpp"
#include "minimizer.hpp"
#include "root_finder.hpp"
#include "threshold_loop_functions.hpp"


#include <cmath>
#include <limits>

namespace flexiblesusy {

#define DERIVEDPARAMETER(p) model->p()
#define EXTRAPARAMETER(p) model->get_##p()
#define INPUTPARAMETER(p) model->get_input().p
#define MODELPARAMETER(p) model->get_##p()
#define PHASE(p) model->get_##p()
#define BETAPARAMETER(p) beta_functions.get_##p()
#define BETAPARAMETER1(l,p) beta_functions_##l##L.get_##p()
#define BETA(p) beta_##p
#define BETA1(l,p) beta_##l##L_##p
#define LowEnergyGaugeCoupling(i) new_g##i
#define LowEnergyConstant(p) Electroweak_constants::p
#define MZPole qedqcd.displayPoleMZ()
#define STANDARDDEVIATION(p) Electroweak_constants::Error_##p
#define Pole(p) model->get_physical().p
#define SCALE model->get_scale()
#define MODEL model
#define MODELCLASSNAME CNMSSMEFTHiggs<Semi_analytic>
#define MWMSbar mW_run
#define MWDRbar mW_run
#define MZMSbar mZ_run
#define MZDRbar mZ_run
#define EDRbar e_run
#define EMSbar e_run
#define CKM ckm
#define PMNS pmns
#define THETAW theta_w
#define THRESHOLD static_cast<int>(model->get_thresholds())
#define ALPHA_EM_DRBAR alpha_em_drbar
#define CALCULATE_DRBAR_MASSES() model->calculate_DRbar_masses()

CNMSSMEFTHiggs_low_scale_constraint<Semi_analytic>::CNMSSMEFTHiggs_low_scale_constraint(
   CNMSSMEFTHiggs<Semi_analytic>* model_, const softsusy::QedQcd& qedqcd_)
   : model(model_)
   , qedqcd(qedqcd_)
{
   initialize();
}

double CNMSSMEFTHiggs_low_scale_constraint<Semi_analytic>::get_scale() const
{
   return scale;
}

void CNMSSMEFTHiggs_low_scale_constraint<Semi_analytic>::apply()
{
   check_model_ptr();

   


   model->calculate_DRbar_masses();
   update_scale();
   qedqcd.run_to(scale, 1.0e-5);
   calculate_DRbar_gauge_couplings();
   calculate_running_SM_masses();

   
   MODEL->set_g1(new_g1);
   MODEL->set_g2(new_g2);
   MODEL->set_g3(new_g3);
   if (is_initial_guess) {
      const auto TanBeta = INPUTPARAMETER(TanBeta);
      const auto vd = MODELPARAMETER(vd);
      const auto vu = MODELPARAMETER(vu);

      MODEL->set_vu(Re(TanBeta*Sqrt((Sqr(vd) + Sqr(vu))/(1 + Sqr(TanBeta)))));
      MODEL->set_vd(Re(Sqrt((Sqr(vd) + Sqr(vu))/(1 + Sqr(TanBeta)))));
   }




}

void CNMSSMEFTHiggs_low_scale_constraint<Semi_analytic>::set_model(Model* model_)
{
   model = cast_model<CNMSSMEFTHiggs<Semi_analytic>*>(model_);
}

void CNMSSMEFTHiggs_low_scale_constraint<Semi_analytic>::clear()
{
   scale = 0.;
   initial_scale_guess = 0.;
   model = nullptr;
   qedqcd = softsusy::QedQcd();
   ckm.setIdentity();
   pmns.setIdentity();
   upQuarksDRbar.setZero();
   downQuarksDRbar.setZero();
   downLeptonsDRbar.setZero();
   neutrinoDRbar.setZero();
   mW_run = 0.;
   mZ_run = 0.;
   AlphaS = 0.;
   e_run = 0.;
   ThetaWDRbar = 0.;
   new_g1 = 0.;
   new_g2 = 0.;
   new_g3 = 0.;
   is_initial_guess = false;
}

void CNMSSMEFTHiggs_low_scale_constraint<Semi_analytic>::initialize()
{
   check_model_ptr();

   initial_scale_guess = qedqcd.displayPoleMZ();

   scale = initial_scale_guess;

   ckm = qedqcd.get_complex_ckm();
   pmns = qedqcd.get_complex_pmns();
   upQuarksDRbar.setZero();
   downQuarksDRbar.setZero();
   downLeptonsDRbar.setZero();
   neutrinoDRbar.setZero();
   mW_run = 0.;
   mZ_run = 0.;
   AlphaS = 0.;
   e_run = 0.;
   ThetaWDRbar = 0.;
   new_g1 = 0.;
   new_g2 = 0.;
   new_g3 = 0.;
   is_initial_guess = false;
}

void CNMSSMEFTHiggs_low_scale_constraint<Semi_analytic>::update_scale()
{
   check_model_ptr();

   scale = qedqcd.displayPoleMZ();


}

void CNMSSMEFTHiggs_low_scale_constraint<Semi_analytic>::calculate_threshold_corrections()
{
   check_model_ptr();

   if (qedqcd.get_scale() != get_scale())
      throw SetupError("Error: low-energy data"
          " set is not defined at the same scale as the low-energy"
          " constraint.  You need to run the low-energy data set to this"
          " scale!");

   const double alpha_em = qedqcd.displayAlpha(softsusy::ALPHA);
   const double alpha_s  = qedqcd.displayAlpha(softsusy::ALPHAS);
   const double mw_pole  = qedqcd.displayPoleMW();
   const double mz_pole  = qedqcd.displayPoleMZ();

   double delta_alpha_em = 0.;
   double delta_alpha_s  = 0.;

   if (model->get_thresholds() && model->get_threshold_corrections().alpha_em > 0)
      delta_alpha_em = calculate_delta_alpha_em(alpha_em);

   if (model->get_thresholds() && model->get_threshold_corrections().alpha_s > 0)
      delta_alpha_s  = calculate_delta_alpha_s(alpha_s);

   const double alpha_em_drbar = alpha_em / (1.0 - delta_alpha_em);
   const double alpha_s_drbar  = alpha_s  / (1.0 - delta_alpha_s);
   const double e_drbar        = Sqrt(4.0 * Pi * alpha_em_drbar);

   // interface variables
   mZ_run = mz_pole;
   mW_run = mw_pole;

   if (model->get_thresholds() && model->get_threshold_corrections().mz > 0)
      mZ_run = model->calculate_MVZ_DRbar(mz_pole);

   if (model->get_thresholds() && model->get_threshold_corrections().mw > 0)
      mW_run = model->calculate_MVWm_DRbar(mw_pole);

   AlphaS = alpha_s_drbar;
   e_run = e_drbar;
   ThetaWDRbar = calculate_theta_w();
}

double CNMSSMEFTHiggs_low_scale_constraint<Semi_analytic>::calculate_theta_w()
{
   check_model_ptr();

   double theta_w = std::asin(Electroweak_constants::sinThetaW);

   CNMSSMEFTHiggs_weinberg_angle::Sm_parameters sm_pars;
   sm_pars.fermi_constant = qedqcd.displayFermiConstant();
   sm_pars.mw_pole = qedqcd.displayPoleMW();
   sm_pars.mz_pole = qedqcd.displayPoleMZ();
   sm_pars.mt_pole = qedqcd.displayPoleMt();
   sm_pars.alpha_s = calculate_alpha_s_SM5_at(qedqcd, qedqcd.displayPoleMt());

   const int number_of_iterations =
       std::max(20, static_cast<int>(std::abs(-log10(MODEL->get_precision()) * 10)
          ));

   CNMSSMEFTHiggs_weinberg_angle weinberg(MODEL, sm_pars);
   weinberg.set_number_of_loops(MODEL->get_threshold_corrections().sin_theta_w);
   weinberg.set_number_of_iterations(number_of_iterations);

   try {
      const auto result = weinberg.calculate();
      THETAW = ArcSin(result.first);

      if (MODEL->get_thresholds() && MODEL->get_threshold_corrections().
         sin_theta_w > 0)
         qedqcd.setPoleMW(result.second);

      MODEL->get_problems().unflag_no_sinThetaW_convergence();
   } catch (const Error& e) {
      VERBOSE_MSG(e.what());
      MODEL->get_problems().flag_no_sinThetaW_convergence();
   }

   return theta_w;
}

void CNMSSMEFTHiggs_low_scale_constraint<Semi_analytic>::calculate_DRbar_gauge_couplings()
{
   check_model_ptr();
   calculate_threshold_corrections();

   new_g1 = 1.2909944487358056*EDRbar*Sec(ThetaWDRbar);
   new_g2 = EDRbar*Csc(ThetaWDRbar);
   new_g3 = 3.5449077018110318*Sqrt(AlphaS);

   if (IsFinite(new_g1)) {
      model->get_problems().unflag_non_perturbative_parameter(CNMSSMEFTHiggs_info
         ::g1);
   } else {
      model->get_problems().flag_non_perturbative_parameter(
         CNMSSMEFTHiggs_info::g1, new_g1, get_scale());
      new_g1 = Electroweak_constants::g1;
   }

   if (IsFinite(new_g2)) {
      model->get_problems().unflag_non_perturbative_parameter(CNMSSMEFTHiggs_info
         ::g2);
   } else {
      model->get_problems().flag_non_perturbative_parameter(
         CNMSSMEFTHiggs_info::g2, new_g2, get_scale());
      new_g2 = Electroweak_constants::g2;
   }
}

double CNMSSMEFTHiggs_low_scale_constraint<Semi_analytic>::calculate_delta_alpha_em(double alphaEm) const
{
   check_model_ptr();

   const double currentScale = model->get_scale();
   const auto MCha = MODELPARAMETER(MCha);
   const auto MHpm = MODELPARAMETER(MHpm);
   const auto MSd = MODELPARAMETER(MSd);
   const auto MSe = MODELPARAMETER(MSe);
   const auto MSu = MODELPARAMETER(MSu);
   const auto MFu = MODELPARAMETER(MFu);

   const double delta_alpha_em_SM = -0.28294212105225836*alphaEm*FiniteLog(Abs(MFu
      (2)/currentScale));

   const double delta_alpha_em = 0.15915494309189535*alphaEm*(0.3333333333333333 -
      1.3333333333333333*FiniteLog(Abs(MCha(0)/currentScale)) - 1.3333333333333333
      *FiniteLog(Abs(MCha(1)/currentScale)) - 0.3333333333333333*FiniteLog(Abs(
      MHpm(1)/currentScale)) - 0.1111111111111111*FiniteLog(Abs(MSd(0)/
      currentScale)) - 0.1111111111111111*FiniteLog(Abs(MSd(1)/currentScale)) -
      0.1111111111111111*FiniteLog(Abs(MSd(2)/currentScale)) - 0.1111111111111111*
      FiniteLog(Abs(MSd(3)/currentScale)) - 0.1111111111111111*FiniteLog(Abs(MSd(4
      )/currentScale)) - 0.1111111111111111*FiniteLog(Abs(MSd(5)/currentScale)) -
      0.3333333333333333*FiniteLog(Abs(MSe(0)/currentScale)) - 0.3333333333333333*
      FiniteLog(Abs(MSe(1)/currentScale)) - 0.3333333333333333*FiniteLog(Abs(MSe(2
      )/currentScale)) - 0.3333333333333333*FiniteLog(Abs(MSe(3)/currentScale)) -
      0.3333333333333333*FiniteLog(Abs(MSe(4)/currentScale)) - 0.3333333333333333*
      FiniteLog(Abs(MSe(5)/currentScale)) - 0.4444444444444444*FiniteLog(Abs(MSu(0
      )/currentScale)) - 0.4444444444444444*FiniteLog(Abs(MSu(1)/currentScale)) -
      0.4444444444444444*FiniteLog(Abs(MSu(2)/currentScale)) - 0.4444444444444444*
      FiniteLog(Abs(MSu(3)/currentScale)) - 0.4444444444444444*FiniteLog(Abs(MSu(4
      )/currentScale)) - 0.4444444444444444*FiniteLog(Abs(MSu(5)/currentScale)));

   return delta_alpha_em + delta_alpha_em_SM;

}

double CNMSSMEFTHiggs_low_scale_constraint<Semi_analytic>::calculate_delta_alpha_s(double alphaS) const
{
   check_model_ptr();

   const double currentScale = model->get_scale();
   const auto MSd = MODELPARAMETER(MSd);
   const auto MSu = MODELPARAMETER(MSu);
   const auto MFu = MODELPARAMETER(MFu);
   const auto MGlu = MODELPARAMETER(MGlu);

   const double delta_alpha_s_SM = -0.1061032953945969*alphaS*FiniteLog(Abs(MFu(2)
      /currentScale));

   const double delta_alpha_s = 0.15915494309189535*alphaS*(0.5 - 2*FiniteLog(Abs(
      MGlu/currentScale)) - 0.16666666666666666*FiniteLog(Abs(MSd(0)/currentScale)
      ) - 0.16666666666666666*FiniteLog(Abs(MSd(1)/currentScale)) -
      0.16666666666666666*FiniteLog(Abs(MSd(2)/currentScale)) -
      0.16666666666666666*FiniteLog(Abs(MSd(3)/currentScale)) -
      0.16666666666666666*FiniteLog(Abs(MSd(4)/currentScale)) -
      0.16666666666666666*FiniteLog(Abs(MSd(5)/currentScale)) -
      0.16666666666666666*FiniteLog(Abs(MSu(0)/currentScale)) -
      0.16666666666666666*FiniteLog(Abs(MSu(1)/currentScale)) -
      0.16666666666666666*FiniteLog(Abs(MSu(2)/currentScale)) -
      0.16666666666666666*FiniteLog(Abs(MSu(3)/currentScale)) -
      0.16666666666666666*FiniteLog(Abs(MSu(4)/currentScale)) -
      0.16666666666666666*FiniteLog(Abs(MSu(5)/currentScale)));

   const double delta_alpha_s_1loop = delta_alpha_s + delta_alpha_s_SM;
   double delta_alpha_s_2loop = 0.;
   double delta_alpha_s_3loop = 0.;

   return delta_alpha_s_1loop + delta_alpha_s_2loop + delta_alpha_s_3loop;

}

double CNMSSMEFTHiggs_low_scale_constraint<Semi_analytic>::calculate_alpha_s_SM5_at(
   softsusy::QedQcd qedqcd_tmp, double scale) const
{
   qedqcd_tmp.run_to(scale); // running in SM(5)
   return qedqcd_tmp.displayAlpha(softsusy::ALPHAS);
}

void CNMSSMEFTHiggs_low_scale_constraint<Semi_analytic>::calculate_running_SM_masses()
{
   check_model_ptr();

   upQuarksDRbar.setZero();
   upQuarksDRbar(0,0) = qedqcd.displayMass(softsusy::mUp);
   upQuarksDRbar(1,1) = qedqcd.displayMass(softsusy::mCharm);
   upQuarksDRbar(2,2) = qedqcd.displayPoleMt();

   downQuarksDRbar.setZero();
   downQuarksDRbar(0,0) = qedqcd.displayMass(softsusy::mDown);
   downQuarksDRbar(1,1) = qedqcd.displayMass(softsusy::mStrange);
   downQuarksDRbar(2,2) = qedqcd.displayMass(softsusy::mBottom);

   downLeptonsDRbar.setZero();
   downLeptonsDRbar(0,0) = qedqcd.displayPoleMel();
   downLeptonsDRbar(1,1) = qedqcd.displayPoleMmuon();
   downLeptonsDRbar(2,2) = qedqcd.displayPoleMtau();

   neutrinoDRbar.setZero();
   neutrinoDRbar(0,0) = qedqcd.displayNeutrinoPoleMass(1);
   neutrinoDRbar(1,1) = qedqcd.displayNeutrinoPoleMass(2);
   neutrinoDRbar(2,2) = qedqcd.displayNeutrinoPoleMass(3);

   if (model->get_thresholds() && model->get_threshold_corrections().mt > 0) {
      upQuarksDRbar(2,2) = MODEL->calculate_MFu_DRbar(qedqcd.displayPoleMt(), 2);
   }

   if (model->get_thresholds() && model->get_threshold_corrections().mb > 0) {
      downQuarksDRbar(2,2) = MODEL->calculate_MFd_DRbar(qedqcd.displayMass(softsusy::mBottom), 2);
   }

   if (model->get_thresholds()) {
      downLeptonsDRbar(0,0) = MODEL->calculate_MFe_DRbar(qedqcd.displayMass(softsusy::mElectron), 0);
      downLeptonsDRbar(1,1) = MODEL->calculate_MFe_DRbar(qedqcd.displayMass(softsusy::mMuon), 1);
   }

   if (model->get_thresholds() && model->get_threshold_corrections().mtau > 0) {
      downLeptonsDRbar(2,2) = MODEL->calculate_MFe_DRbar(qedqcd.displayMass(softsusy::mTau), 2);
   }
}

void CNMSSMEFTHiggs_low_scale_constraint<Semi_analytic>::calculate_DRbar_yukawa_couplings()
{
   calculate_running_SM_masses();
   calculate_Yu_DRbar();
   calculate_Yd_DRbar();
   calculate_Ye_DRbar();
}

void CNMSSMEFTHiggs_low_scale_constraint<Semi_analytic>::calculate_Yu_DRbar()
{
   check_model_ptr();


}

void CNMSSMEFTHiggs_low_scale_constraint<Semi_analytic>::calculate_Yd_DRbar()
{
   check_model_ptr();


}

void CNMSSMEFTHiggs_low_scale_constraint<Semi_analytic>::calculate_Ye_DRbar()
{
   check_model_ptr();


}

void CNMSSMEFTHiggs_low_scale_constraint<Semi_analytic>::check_model_ptr() const
{
   if (!model)
      throw SetupError("CNMSSMEFTHiggs_low_scale_constraint<Semi_analytic>: "
                       "model pointer is zero!");
}

} // namespace flexiblesusy
