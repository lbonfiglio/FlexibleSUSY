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

// File generated at Wed 10 Jul 2019 12:34:50

/**
 * @file NMSSM_two_scale_ewsb_solver.cpp
 *
 * @brief implementation of EWSB solver for two-scale iteration
 *
 * This file was generated at Wed 10 Jul 2019 12:34:50 with FlexibleSUSY
 * 2.3.0 (git commit: dd385b4d85d1158da96abc8a66b31791e0c651c9) and SARAH 4.12.3 .
 */

#include "NMSSM_two_scale_ewsb_solver.hpp"
#include "NMSSM_mass_eigenstates.hpp"
#include "logger.hpp"
#include "root_finder.hpp"
#include "fixed_point_iterator.hpp"
#include "raii.hpp"

#include <memory>

namespace flexiblesusy {

#define INPUT(parameter) model.get_input().parameter
#define LOCALINPUT(parameter) INPUT(parameter)
#define INPUTPARAMETER(parameter) LOCALINPUT(parameter)
#define MODELPARAMETER(parameter) model.get_##parameter()
#define EXTRAPARAMETER(parameter) model.get_##parameter()
#define DERIVEDPARAMETER(parameter) model.parameter()
#define PHASE(p) MODELPARAMETER(p)
#define LowEnergyConstant(p) Electroweak_constants::p
#define CLASSNAME NMSSM_ewsb_solver<Two_scale>

/**
 * This method solves the EWSB conditions iteratively, trying several
 * root finding methods until a solution is found.
 */
int CLASSNAME::solve_iteratively(NMSSM_mass_eigenstates& model_to_solve)
{
   auto model = model_to_solve;
   model.set_ewsb_loop_order(loop_order);

   auto ewsb_stepper = [this, model](const EWSB_vector_t& ewsb_pars) mutable -> EWSB_vector_t {
      const double Kappa = ewsb_pars(0);
      const double vS = INPUT(SignvS)*Abs(ewsb_pars(1));
      const double ms2 = ewsb_pars(2);

      model.set_Kappa(Kappa);
      model.set_vS(vS);
      model.set_ms2(ms2);


      if (this->loop_order > 0)
         model.calculate_DRbar_masses();

      return this->ewsb_step(model);
   };

   auto tadpole_stepper = [this, model](const EWSB_vector_t& ewsb_pars) mutable -> EWSB_vector_t {
      const double Kappa = ewsb_pars(0);
      const double vS = INPUT(SignvS)*Abs(ewsb_pars(1));
      const double ms2 = ewsb_pars(2);

      model.set_Kappa(Kappa);
      model.set_vS(vS);
      model.set_ms2(ms2);


      if (this->loop_order > 0)
         model.calculate_DRbar_masses();

      return this->tadpole_equations(model);
   };

   std::unique_ptr<EWSB_solver> solvers[] = {
      std::unique_ptr<EWSB_solver>(new Fixed_point_iterator<number_of_ewsb_equations, fixed_point_iterator::Convergence_tester_relative>(ewsb_stepper, number_of_iterations, fixed_point_iterator::Convergence_tester_relative(precision))),
      std::unique_ptr<EWSB_solver>(new Root_finder<number_of_ewsb_equations>(tadpole_stepper, number_of_iterations, precision, Root_finder<number_of_ewsb_equations>::GSLHybridS)),
      std::unique_ptr<EWSB_solver>(new Root_finder<number_of_ewsb_equations>(tadpole_stepper, number_of_iterations, precision, Root_finder<number_of_ewsb_equations>::GSLBroyden))
   };

   const auto x_init(initial_guess(model_to_solve));

   VERBOSE_MSG("\t\tSolving EWSB equations ...");
   VERBOSE_MSG("\t\tInitial guess: x_init = " << x_init.transpose());

   int status;
   for (auto& solver: solvers) {
      VERBOSE_MSG("\t\t\tStarting EWSB iteration using " << solver->name());
      status = solve_iteratively_with(model_to_solve, solver.get(), x_init);
      if (status == EWSB_solver::SUCCESS) {
         VERBOSE_MSG("\t\t\t" << solver->name() << " finished successfully!");
         break;
      }
#ifdef ENABLE_VERBOSE
      else {
         WARNING("\t\t\t" << solver->name() << " could not find a solution!"
                 " (requested precision: " << precision << ")");
      }
#endif
   }

   if (status == EWSB_solver::SUCCESS) {
      model_to_solve.get_problems().unflag_no_ewsb();
   } else {
      set_best_ewsb_solution(model_to_solve, std::begin(solvers), std::end(solvers));
      model_to_solve.get_problems().flag_no_ewsb();
#ifdef ENABLE_VERBOSE
      WARNING("\t\tCould not find a solution to the EWSB equations!"
              " (requested precision: " << precision << ")");
#endif
   }

   return status;
}

/**
 * Solves EWSB equations with given EWSB solver
 *
 * @param model_to_solve model to solve EWSB conditions for
 * @param solver EWSB solver
 * @param x_init initial values
 *
 * @return status of the EWSB solver
 */
int CLASSNAME::solve_iteratively_with(
   NMSSM_mass_eigenstates& model_to_solve, EWSB_solver* solver, const EWSB_vector_t& x_init)
{
   const int status = solver->solve(x_init);

   if (status == EWSB_solver::SUCCESS)
      set_ewsb_solution(model_to_solve, solver);

   return status;
}

/**
 * Sets EWSB output parameters from given solver.
 *
 * @param model model to set EWSB output parameters in
 * @param solver solver
 */
void CLASSNAME::set_ewsb_solution(NMSSM_mass_eigenstates& model, const EWSB_solver* solver)
{
   const auto solution = solver->get_solution();

   const double Kappa = solution(0);
   const double vS = INPUT(SignvS)*Abs(solution(1));
   const double ms2 = solution(2);
   model.set_Kappa(Kappa);
   model.set_vS(vS);
   model.set_ms2(ms2);


   model.calculate_DRbar_masses();
}

/**
 * Sets EWSB output parameters from the solver from the range [first,
 * last), which minimizes the tadpole equations at most.
 *
 * @param model model to set EWSB output parameters in
 * @param first iterator to first solver
 * @param last iterator to last solver
 */
template <typename It>
void CLASSNAME::set_best_ewsb_solution(NMSSM_mass_eigenstates& model, It first, It last)
{
   auto ma(model), mb(model);

   const auto best_solver =
      std::min_element(first, last,
                       [this, &ma, &mb](const std::unique_ptr<EWSB_solver>& a, const std::unique_ptr<EWSB_solver>& b) {
                          this->set_ewsb_solution(ma, a.get());
                          this->set_ewsb_solution(mb, b.get());
                          return Total(Abs(Re(ma.tadpole_equations()))) < Total(Abs(Re(mb.tadpole_equations())));
                       });

   VERBOSE_MSG("\t\tUsing best solution from " << (*best_solver)->name());

   set_ewsb_solution(model, best_solver->get());
}

int CLASSNAME::solve_iteratively_at(NMSSM_mass_eigenstates& model_to_solve, int l)
{
   // temporarily set `ewsb_loop_order' to `loop_order' and do
   // iteration
   const auto save_loop_order_raii = make_raii_save(loop_order);
   loop_order = l;

   return solve_iteratively(model_to_solve);
}

int CLASSNAME::solve(NMSSM_mass_eigenstates& model_to_solve)
{
   if (loop_order == 0) {
      return solve_tree_level(model_to_solve);
   }
   return solve_iteratively_at(model_to_solve, loop_order);
}

int CLASSNAME::solve_tree_level(NMSSM_mass_eigenstates& model)
{
   int error = EWSB_solver::SUCCESS;

   const auto TLambdax = MODELPARAMETER(TLambdax);
   const auto TKappa = MODELPARAMETER(TKappa);
   const auto mHd2 = MODELPARAMETER(mHd2);
   const auto vd = MODELPARAMETER(vd);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto mHu2 = MODELPARAMETER(mHu2);
   const auto vu = MODELPARAMETER(vu);
   const auto Lambdax = MODELPARAMETER(Lambdax);

   double vS;
   double Kappa;
   double ms2;

   vS = Re(0.22360679774997896*LOCALINPUT(SignvS)*Sqrt((-3*Quad(vd)*Sqr(g1) + 3*
      Quad(vu)*Sqr(g1) - 5*Quad(vd)*Sqr(g2) + 5*Quad(vu)*Sqr(g2) - 40*mHd2*Sqr(vd)
      + 40*mHu2*Sqr(vu))/(AbsSqr(Lambdax)*Sqr(vd) - AbsSqr(Lambdax)*Sqr(vu))));
   Kappa = Re((0.1*(40*mHd2*vd - 14.142135623730951*vS*vu*Conj(TLambdax) + 3*Cube(
      vd)*Sqr(g1) + 5*Cube(vd)*Sqr(g2) + 20*vd*AbsSqr(Lambdax)*Sqr(vS) + 20*vd*
      AbsSqr(Lambdax)*Sqr(vu) - 3*vd*Sqr(g1)*Sqr(vu) - 5*vd*Sqr(g2)*Sqr(vu) -
      14.142135623730951*vS*vu*TLambdax))/(vu*(Conj(Lambdax) + Lambdax)*Sqr(vS)));
   ms2 = Re((0.25*(2*Kappa*vd*vS*vu*Conj(Lambdax) + 1.4142135623730951*vd*vu*Conj(
      TLambdax) + 2*Kappa*vd*vS*vu*Lambdax - 4*Cube(vS)*Sqr(Kappa) - 2*vS*AbsSqr(
      Lambdax)*Sqr(vd) - 1.4142135623730951*Conj(TKappa)*Sqr(vS) - 2*vS*AbsSqr(
      Lambdax)*Sqr(vu) - 1.4142135623730951*Sqr(vS)*TKappa + 1.4142135623730951*vd
      *vu*TLambdax))/vS);

   
   const bool is_finite = IsFinite(vS) && IsFinite(Kappa) && IsFinite(ms2);

   if (is_finite) {
      model.set_vS(vS);
      model.set_Kappa(Kappa);
      model.set_ms2(ms2);
      model.get_problems().unflag_no_ewsb_tree_level();
   } else {
      error = EWSB_solver::FAIL;
      model.get_problems().flag_no_ewsb_tree_level();
   }
   return error;
}

CLASSNAME::EWSB_vector_t CLASSNAME::initial_guess(const NMSSM_mass_eigenstates& model) const
{
   EWSB_vector_t x_init(EWSB_vector_t::Zero());

   const auto Kappa = MODELPARAMETER(Kappa);
   const auto vS = MODELPARAMETER(vS);
   const auto ms2 = MODELPARAMETER(ms2);
   x_init[0] = Kappa;
   x_init[1] = vS;
   x_init[2] = ms2;


   return x_init;
}

CLASSNAME::EWSB_vector_t CLASSNAME::tadpole_equations(const NMSSM_mass_eigenstates& model) const
{
   return model.tadpole_equations();
}

/**
 * Calculates EWSB output parameters including loop-corrections.
 *
 * Throws exception of type EEWSBStepFailed if new EWSB parameters are
 * inf or nan.
 *
 * @param model model to use for calculating the EWSB output parameters
 *
 * @return new set of EWSB output parameters
 */
CLASSNAME::EWSB_vector_t CLASSNAME::ewsb_step(const NMSSM_mass_eigenstates& model) const
{
   std::array<double, number_of_ewsb_equations> tadpole{};
   EWSB_vector_t ewsb_parameters(EWSB_vector_t::Zero());

   if (loop_order > 0) {
   tadpole[0] += Re(model.tadpole_hh_1loop(0));
   tadpole[1] += Re(model.tadpole_hh_1loop(1));
   tadpole[2] += Re(model.tadpole_hh_1loop(2));

      if (loop_order > 1) {
   const auto tadpole_2l(model.tadpole_hh_2loop());
   tadpole[0] += tadpole_2l(0);
   tadpole[1] += tadpole_2l(1);
   tadpole[2] += tadpole_2l(2);

      }
   }

   const auto TLambdax = MODELPARAMETER(TLambdax);
   const auto TKappa = MODELPARAMETER(TKappa);
   const auto vd = MODELPARAMETER(vd);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto vu = MODELPARAMETER(vu);
   const auto mHd2 = MODELPARAMETER(mHd2);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto mHu2 = MODELPARAMETER(mHu2);
   double vS;
   double Kappa;
   double ms2;

   vS = Re(0.22360679774997896*LOCALINPUT(SignvS)*Sqrt((40*vd*tadpole[0] - 40*vu*
      tadpole[1] - 3*Quad(vd)*Sqr(g1) + 3*Quad(vu)*Sqr(g1) - 5*Quad(vd)*Sqr(g2) +
      5*Quad(vu)*Sqr(g2) - 40*mHd2*Sqr(vd) + 40*mHu2*Sqr(vu))/(AbsSqr(Lambdax)*Sqr
      (vd) - AbsSqr(Lambdax)*Sqr(vu))));
   Kappa = Re((0.1*(40*mHd2*vd - 14.142135623730951*vS*vu*Conj(TLambdax) - 40*
      tadpole[0] + 3*Cube(vd)*Sqr(g1) + 5*Cube(vd)*Sqr(g2) + 20*vd*AbsSqr(Lambdax)
      *Sqr(vS) + 20*vd*AbsSqr(Lambdax)*Sqr(vu) - 3*vd*Sqr(g1)*Sqr(vu) - 5*vd*Sqr(
      g2)*Sqr(vu) - 14.142135623730951*vS*vu*TLambdax))/(vu*(Conj(Lambdax) +
      Lambdax)*Sqr(vS)));
   ms2 = Re((0.25*(2*Kappa*vd*vS*vu*Conj(Lambdax) + 1.4142135623730951*vd*vu*Conj(
      TLambdax) + 2*Kappa*vd*vS*vu*Lambdax + 4*tadpole[2] - 4*Cube(vS)*Sqr(Kappa)
      - 2*vS*AbsSqr(Lambdax)*Sqr(vd) - 1.4142135623730951*Conj(TKappa)*Sqr(vS) - 2
      *vS*AbsSqr(Lambdax)*Sqr(vu) - 1.4142135623730951*Sqr(vS)*TKappa +
      1.4142135623730951*vd*vu*TLambdax))/vS);

   const bool is_finite = IsFinite(vS) && IsFinite(Kappa) && IsFinite(ms2);


   if (!is_finite)
      throw EEWSBStepFailed();

   ewsb_parameters[0] = Kappa;
   ewsb_parameters[1] = vS;
   ewsb_parameters[2] = ms2;


   return ewsb_parameters;
}

} // namespace flexiblesusy
