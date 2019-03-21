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

// File generated at Thu 21 Mar 2019 14:21:46

/**
 * @file CNMSSMEFTHiggs_semi_analytic_ewsb_solver.cpp
 *
 * @brief implementation of EWSB solver for semi-analytic iteration
 *
 * This file was generated at Thu 21 Mar 2019 14:21:46 with FlexibleSUSY
 * 2.3.0 (git commit: 29358b5665eb4bb485434fc07d3bcacce3e489a2) and SARAH 4.12.3 .
 */

#include "CNMSSMEFTHiggs_semi_analytic_ewsb_solver.hpp"
#include "CNMSSMEFTHiggs_mass_eigenstates.hpp"
#include "CNMSSMEFTHiggs_semi_analytic_solutions.hpp"
#include "error.hpp"
#include "ew_input.hpp"
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
#define SEMIANALYTICPARAMETER(p) solutions->get_##p()
#define LowEnergyConstant(p) Electroweak_constants::p
#define CLASSNAME CNMSSMEFTHiggs_ewsb_solver<Semi_analytic>

/**
 * This method solves the EWSB conditions iteratively, trying several
 * root finding methods until a solution is found.
 */
int CLASSNAME::solve_iteratively(CNMSSMEFTHiggs_mass_eigenstates& model_to_solve)
{
   auto model = model_to_solve;
   model.set_ewsb_loop_order(loop_order);

   auto ewsb_stepper = [this, model](const EWSB_vector_t& ewsb_pars) mutable -> EWSB_vector_t {
      const double Kappa = ewsb_pars(0);
      const double vS = ewsb_pars(1);
      const double m0Sq = ewsb_pars(2);

      model.set_Kappa(Kappa);
      model.set_vS(vS);
      model.set_m0Sq(m0Sq);

      solutions->evaluate_solutions(model);

      if (this->loop_order > 0)
         model.calculate_DRbar_masses();

      return this->ewsb_step(model);
   };

   auto tadpole_stepper = [this, model](const EWSB_vector_t& ewsb_pars) mutable -> EWSB_vector_t {
      const double Kappa = ewsb_pars(0);
      const double vS = ewsb_pars(1);
      const double m0Sq = ewsb_pars(2);

      model.set_Kappa(Kappa);
      model.set_vS(vS);
      model.set_m0Sq(m0Sq);

      solutions->evaluate_solutions(model);

      if (this->loop_order > 0)
         model.calculate_DRbar_masses();

      return this->tadpole_equations(model);
   };

   std::unique_ptr<EWSB_solver> solvers[] = {
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
   CNMSSMEFTHiggs_mass_eigenstates& model_to_solve, EWSB_solver* solver, const EWSB_vector_t& x_init)
{
   const int status = solver->solve(x_init);

   if (status == EWSB_solver::SUCCESS)
      set_ewsb_solution(model_to_solve, solver);

   return status;
}

/**
 * Sets EWSB output parameters from given solver.
 *
 * @param solver solver
 */
void CLASSNAME::set_ewsb_solution(CNMSSMEFTHiggs_mass_eigenstates& model, const EWSB_solver* solver)
{
   const auto solution = solver->get_solution();

   const double Kappa = solution(0);
   const double vS = solution(1);
   const double m0Sq = solution(2);
   model.set_Kappa(Kappa);
   model.set_vS(vS);
   model.set_m0Sq(m0Sq);

   solutions->evaluate_solutions(model);

   model.calculate_DRbar_masses();
}

/**
 * Sets EWSB output parameters from the solver from the range [first,
 * last), which minimizes the tadpole equations at most.
 *
 * @param first iterator to first solver
 * @param last iterator to last solver
 */
template <typename It>
void CLASSNAME::set_best_ewsb_solution(CNMSSMEFTHiggs_mass_eigenstates& model, It first, It last)
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

int CLASSNAME::solve_iteratively_at(CNMSSMEFTHiggs_mass_eigenstates& model_to_solve, int l)
{
   // temporarily set `ewsb_loop_order' to `loop_order' and do
   // iteration
   const auto save_loop_order_raii = make_raii_save(loop_order);
   loop_order = l;

   return solve_iteratively(model_to_solve);
}

int CLASSNAME::solve(CNMSSMEFTHiggs_mass_eigenstates& model_to_solve)
{
   if (!solutions) {
      throw SetupError("CNMSSMEFTHiggs_ewsb_solver<Semi_analytic>:solve: "
                       "pointer to semi-analytic solutions is zero!");
   }

   if (loop_order == 0) {
      return solve_tree_level(model_to_solve);
   }
   return solve_iteratively_at(model_to_solve, loop_order);
}

int CLASSNAME::solve_tree_level(CNMSSMEFTHiggs_mass_eigenstates& model)
{
   int error = EWSB_solver::SUCCESS;

   error = solve_iteratively_at(model, 0);


   return error;
}

CLASSNAME::EWSB_vector_t CLASSNAME::initial_guess(
   const CNMSSMEFTHiggs_mass_eigenstates& model) const
{
   EWSB_vector_t x_init(EWSB_vector_t::Zero());

   const auto Kappa = MODELPARAMETER(Kappa);
   const auto vS = MODELPARAMETER(vS);
   const auto m0Sq = EXTRAPARAMETER(m0Sq);
   x_init[0] = Kappa;
   x_init[1] = vS;
   x_init[2] = m0Sq;


   return x_init;
}

CLASSNAME::EWSB_vector_t CLASSNAME::tadpole_equations(
   const CNMSSMEFTHiggs_mass_eigenstates& model) const
{
   return model.tadpole_equations();
}

/**
 * Calculates EWSB output parameters including loop-corrections.
 *
 * Throws exception of type EEWSBStepFailed if new EWSB parameters are
 * inf or nan.
 *
 * @return new set of EWSB output parameters
 */
CLASSNAME::EWSB_vector_t CLASSNAME::ewsb_step(
   const CNMSSMEFTHiggs_mass_eigenstates& model) const
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

   const bool is_finite = false;


   if (!is_finite)
      throw EEWSBStepFailed();



   return ewsb_parameters;
}

void CLASSNAME::set_semi_analytic_solutions(
   CNMSSMEFTHiggs_semi_analytic_solutions* s)
{
   solutions = s;
}

} // namespace flexiblesusy
