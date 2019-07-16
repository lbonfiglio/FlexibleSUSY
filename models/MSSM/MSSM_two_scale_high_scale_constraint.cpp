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

// File generated at Tue 16 Jul 2019 15:02:20

#include "MSSM_two_scale_high_scale_constraint.hpp"
#include "MSSM_two_scale_model.hpp"
#include "MSSM_info.hpp"
#include "config.h"
#include "wrappers.hpp"
#include "logger.hpp"
#include "ew_input.hpp"
#include "error.hpp"
#include "gsl_utils.hpp"
#include "minimizer.hpp"
#include "raii.hpp"
#include "root_finder.hpp"
#include "threshold_loop_functions.hpp"
#include "numerics2.hpp"

#ifdef ENABLE_HIMALAYA
#include "HierarchyCalculator.hpp"
#include "version.hpp"
#endif

#include <cmath>
#include <cerrno>
#include <cstring>

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
#define LowEnergyConstant(p) Electroweak_constants::p
#define MZPole Electroweak_constants::MZ
#define STANDARDDEVIATION(p) Electroweak_constants::Error_##p
#define Pole(p) model->get_physical().p
#define SCALE model->get_scale()
#define THRESHOLD static_cast<int>(model->get_thresholds())
#define MODEL model
#define MODELCLASSNAME MSSM<Two_scale>

#if defined(ENABLE_HIMALAYA) && Himalaya_VERSION_MAJOR >= 2
#define FSHimalayaMh23L [&] () {                                        \
      MODEL->calculate_DRbar_masses();                                  \
                                                                        \
      himalaya::Parameters pars;                                        \
      const auto TYu = MODELPARAMETER(TYu); \
      const auto TYd = MODELPARAMETER(TYd); \
      const auto TYe = MODELPARAMETER(TYe); \
      const auto Mu = MODELPARAMETER(Mu); \
      const auto g1 = MODELPARAMETER(g1); \
      const auto g2 = MODELPARAMETER(g2); \
      const auto g3 = MODELPARAMETER(g3); \
      const auto vu = MODELPARAMETER(vu); \
      const auto vd = MODELPARAMETER(vd); \
      const auto Yu = MODELPARAMETER(Yu); \
      const auto Yd = MODELPARAMETER(Yd); \
      const auto Ye = MODELPARAMETER(Ye); \
      const auto ZUL = MODELPARAMETER(ZUL); \
      const auto ZDR = MODELPARAMETER(ZDR); \
      const auto ZUR = MODELPARAMETER(ZUR); \
      const auto ZEL = MODELPARAMETER(ZEL); \
      const auto ZER = MODELPARAMETER(ZER); \
      const auto MGlu = MODELPARAMETER(MGlu); \
      const auto MAh = MODELPARAMETER(MAh); \
       \
       \
      pars.scale = MODELPARAMETER(scale); \
      pars.mu = Re(Mu); \
      pars.g1 = Re(g1); \
      pars.g2 = Re(g2); \
      pars.g3 = Re(g3); \
      pars.vd = Re(vd); \
      pars.vu = Re(vu); \
      pars.mq2 = Re(ZUL); \
      pars.md2 = Re(ZDR); \
      pars.mu2 = Re(ZUR); \
      pars.ml2 = Re(ZEL); \
      pars.me2 = Re(ZER); \
      pars.Au = Re(TYu); \
      pars.Ad = Re(TYd); \
      pars.Ae = Re(TYe); \
      pars.Yu = Re(Yu); \
      pars.Yd = Re(Yd); \
      pars.Ye = Re(Ye); \
      pars.M1 = 0; \
      pars.M2 = 0; \
      pars.MG = MGlu; \
      pars.MA = MAh; \
       \
      const double msbar_scheme = 1; \
      const double lambda_3L_eft = 1; \
      const double lambda_3L_uncertainty = 0; \
       \
                                                                        \
      double lambda_3L = 0.;                                            \
                                                                        \
      try {                                                             \
         const bool verbose = false;                                    \
         himalaya::HierarchyCalculator hc(pars, verbose);               \
                                                                        \
         const auto ho = hc.calculateDMh3L(false);                      \
                                                                        \
         lambda_3L =                                                    \
            lambda_3L_eft * (                                           \
               ho.getDLambda(3)                                         \
               + msbar_scheme*ho.getDLambdaDRbarPrimeToMSbarShift(3)    \
               + lambda_3L_uncertainty*ho.getDLambdaUncertainty(3)      \
            );                                                          \
                                                                        \
         VERBOSE_MSG("Himalaya top (hierarchy, Dlambda_3L) = ("         \
                     << ho.getSuitableHierarchy() << ", "               \
                     << lambda_3L <<")");                               \
      } catch (const std::exception& e) {                               \
         model->get_problems().flag_bad_mass(MSSM_info::hh); \
         WARNING(e.what());                                             \
         VERBOSE_MSG(pars);                                             \
      }                                                                 \
                                                                        \
      return lambda_3L;                                                 \
   }()
#else
#define FSHimalayaMh23L [] () {                                         \
      throw HimalayaError("The 3-loop corrections to lambda "           \
                          "require Himalaya 2.0.0 (or higher), but "    \
                          "FlexibleSUSY has not been configured with "  \
                          "this Himalaya version!");                    \
      return 0.;                                                        \
   }()
#endif

MSSM_high_scale_constraint<Two_scale>::MSSM_high_scale_constraint(
   MSSM<Two_scale>* model_)
   : model(model_)
{
   initialize();
}

void MSSM_high_scale_constraint<Two_scale>::apply()
{
   check_model_ptr();

   


   update_scale();

   const auto Aeij = INPUTPARAMETER(Aeij);
   const auto Adij = INPUTPARAMETER(Adij);
   const auto Auij = INPUTPARAMETER(Auij);
   const auto mHd2IN = INPUTPARAMETER(mHd2IN);
   const auto mHu2IN = INPUTPARAMETER(mHu2IN);
   const auto mq2Input = INPUTPARAMETER(mq2Input);
   const auto ml2Input = INPUTPARAMETER(ml2Input);
   const auto md2Input = INPUTPARAMETER(md2Input);
   const auto mu2Input = INPUTPARAMETER(mu2Input);
   const auto me2Input = INPUTPARAMETER(me2Input);
   const auto MassBInput = INPUTPARAMETER(MassBInput);
   const auto MassWBInput = INPUTPARAMETER(MassWBInput);
   const auto MassGInput = INPUTPARAMETER(MassGInput);
   const auto Ye = MODELPARAMETER(Ye);
   const auto Yd = MODELPARAMETER(Yd);
   const auto Yu = MODELPARAMETER(Yu);

   MODEL->set_TYe(((Aeij).cwiseProduct(Ye)).real());
   MODEL->set_TYd(((Adij).cwiseProduct(Yd)).real());
   MODEL->set_TYu(((Auij).cwiseProduct(Yu)).real());
   MODEL->set_mHd2(Re(mHd2IN));
   MODEL->set_mHu2(Re(mHu2IN));
   MODEL->set_mq2((mq2Input).real());
   MODEL->set_ml2((ml2Input).real());
   MODEL->set_md2((md2Input).real());
   MODEL->set_mu2((mu2Input).real());
   MODEL->set_me2((me2Input).real());
   MODEL->set_MassB(Re(MassBInput));
   MODEL->set_MassWB(Re(MassWBInput));
   MODEL->set_MassG(Re(MassGInput));


   check_non_perturbative();
}

bool MSSM_high_scale_constraint<Two_scale>::check_non_perturbative()
{
   bool problem = false;

   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto g3 = MODELPARAMETER(g3);
   const auto Yd = MODELPARAMETER(Yd);
   const auto Ye = MODELPARAMETER(Ye);
   const auto Yu = MODELPARAMETER(Yu);

   if (MaxAbsValue(g1) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(MSSM_info::g1, MaxAbsValue(g1), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(MSSM_info::g1);
   }
   if (MaxAbsValue(g2) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(MSSM_info::g2, MaxAbsValue(g2), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(MSSM_info::g2);
   }
   if (MaxAbsValue(g3) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(MSSM_info::g3, MaxAbsValue(g3), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(MSSM_info::g3);
   }
   if (MaxAbsValue(Yd(0,0)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(MSSM_info::Yd0_0, MaxAbsValue(Yd(0,0)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(MSSM_info::Yd0_0);
   }

   if (MaxAbsValue(Yd(0,1)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(MSSM_info::Yd0_1, MaxAbsValue(Yd(0,1)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(MSSM_info::Yd0_1);
   }

   if (MaxAbsValue(Yd(0,2)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(MSSM_info::Yd0_2, MaxAbsValue(Yd(0,2)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(MSSM_info::Yd0_2);
   }

   if (MaxAbsValue(Yd(1,0)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(MSSM_info::Yd1_0, MaxAbsValue(Yd(1,0)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(MSSM_info::Yd1_0);
   }

   if (MaxAbsValue(Yd(1,1)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(MSSM_info::Yd1_1, MaxAbsValue(Yd(1,1)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(MSSM_info::Yd1_1);
   }

   if (MaxAbsValue(Yd(1,2)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(MSSM_info::Yd1_2, MaxAbsValue(Yd(1,2)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(MSSM_info::Yd1_2);
   }

   if (MaxAbsValue(Yd(2,0)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(MSSM_info::Yd2_0, MaxAbsValue(Yd(2,0)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(MSSM_info::Yd2_0);
   }

   if (MaxAbsValue(Yd(2,1)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(MSSM_info::Yd2_1, MaxAbsValue(Yd(2,1)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(MSSM_info::Yd2_1);
   }

   if (MaxAbsValue(Yd(2,2)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(MSSM_info::Yd2_2, MaxAbsValue(Yd(2,2)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(MSSM_info::Yd2_2);
   }
   if (MaxAbsValue(Ye(0,0)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(MSSM_info::Ye0_0, MaxAbsValue(Ye(0,0)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(MSSM_info::Ye0_0);
   }

   if (MaxAbsValue(Ye(0,1)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(MSSM_info::Ye0_1, MaxAbsValue(Ye(0,1)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(MSSM_info::Ye0_1);
   }

   if (MaxAbsValue(Ye(0,2)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(MSSM_info::Ye0_2, MaxAbsValue(Ye(0,2)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(MSSM_info::Ye0_2);
   }

   if (MaxAbsValue(Ye(1,0)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(MSSM_info::Ye1_0, MaxAbsValue(Ye(1,0)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(MSSM_info::Ye1_0);
   }

   if (MaxAbsValue(Ye(1,1)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(MSSM_info::Ye1_1, MaxAbsValue(Ye(1,1)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(MSSM_info::Ye1_1);
   }

   if (MaxAbsValue(Ye(1,2)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(MSSM_info::Ye1_2, MaxAbsValue(Ye(1,2)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(MSSM_info::Ye1_2);
   }

   if (MaxAbsValue(Ye(2,0)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(MSSM_info::Ye2_0, MaxAbsValue(Ye(2,0)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(MSSM_info::Ye2_0);
   }

   if (MaxAbsValue(Ye(2,1)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(MSSM_info::Ye2_1, MaxAbsValue(Ye(2,1)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(MSSM_info::Ye2_1);
   }

   if (MaxAbsValue(Ye(2,2)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(MSSM_info::Ye2_2, MaxAbsValue(Ye(2,2)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(MSSM_info::Ye2_2);
   }
   if (MaxAbsValue(Yu(0,0)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(MSSM_info::Yu0_0, MaxAbsValue(Yu(0,0)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(MSSM_info::Yu0_0);
   }

   if (MaxAbsValue(Yu(0,1)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(MSSM_info::Yu0_1, MaxAbsValue(Yu(0,1)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(MSSM_info::Yu0_1);
   }

   if (MaxAbsValue(Yu(0,2)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(MSSM_info::Yu0_2, MaxAbsValue(Yu(0,2)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(MSSM_info::Yu0_2);
   }

   if (MaxAbsValue(Yu(1,0)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(MSSM_info::Yu1_0, MaxAbsValue(Yu(1,0)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(MSSM_info::Yu1_0);
   }

   if (MaxAbsValue(Yu(1,1)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(MSSM_info::Yu1_1, MaxAbsValue(Yu(1,1)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(MSSM_info::Yu1_1);
   }

   if (MaxAbsValue(Yu(1,2)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(MSSM_info::Yu1_2, MaxAbsValue(Yu(1,2)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(MSSM_info::Yu1_2);
   }

   if (MaxAbsValue(Yu(2,0)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(MSSM_info::Yu2_0, MaxAbsValue(Yu(2,0)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(MSSM_info::Yu2_0);
   }

   if (MaxAbsValue(Yu(2,1)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(MSSM_info::Yu2_1, MaxAbsValue(Yu(2,1)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(MSSM_info::Yu2_1);
   }

   if (MaxAbsValue(Yu(2,2)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(MSSM_info::Yu2_2, MaxAbsValue(Yu(2,2)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(MSSM_info::Yu2_2);
   }


   return problem;
}

double MSSM_high_scale_constraint<Two_scale>::get_scale() const
{
   return scale;
}

double MSSM_high_scale_constraint<Two_scale>::get_initial_scale_guess() const
{
   return initial_scale_guess;
}

const MSSM_input_parameters& MSSM_high_scale_constraint<Two_scale>::get_input_parameters() const
{
   return model->get_input();
}

MSSM<Two_scale>* MSSM_high_scale_constraint<Two_scale>::get_model() const
{
   return model;
}

void MSSM_high_scale_constraint<Two_scale>::set_model(Model* model_)
{
   model = cast_model<MSSM<Two_scale>*>(model_);
}

void MSSM_high_scale_constraint<Two_scale>::set_scale(double s)
{
   scale = s;
}

void MSSM_high_scale_constraint<Two_scale>::clear()
{
   scale = 0.;
   initial_scale_guess = 0.;
   model = nullptr;
}

void MSSM_high_scale_constraint<Two_scale>::initialize()
{
   check_model_ptr();

   const auto Qin = INPUTPARAMETER(Qin);

   initial_scale_guess = Qin;

   scale = initial_scale_guess;
}

void MSSM_high_scale_constraint<Two_scale>::update_scale()
{
   check_model_ptr();

   const auto Qin = INPUTPARAMETER(Qin);

   scale = Qin;


}

void MSSM_high_scale_constraint<Two_scale>::check_model_ptr() const
{
   if (!model)
      throw SetupError("MSSM_high_scale_constraint<Two_scale>: "
                       "model pointer is zero!");
}

} // namespace flexiblesusy
