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

// File generated at Thu 21 Mar 2019 14:21:44

/**
 * @file CNMSSMEFTHiggs_semi_analytic_solutions.cpp
 * @brief contains implementation of class for computing the semi-analytic solutions
 *
 * This file was generated at Thu 21 Mar 2019 14:21:44 with FlexibleSUSY
 * 2.3.0 (git commit: 29358b5665eb4bb485434fc07d3bcacce3e489a2) and SARAH 4.12.3 .
 */

#include "CNMSSMEFTHiggs_semi_analytic_solutions.hpp"
#include "CNMSSMEFTHiggs_mass_eigenstates.hpp"

#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUTPARAMETER(parameter) model.get_input().parameter
#define MODELPARAMETER(parameter) model.get_##parameter()
#define EXTRAPARAMETER(parameter) model.get_##parameter()
#define DERIVEDPARAMETER(parameter) model.parameter()
#define PHASE(p) MODELPARAMETER(p)
#define BOUNDARYVALUE(p) boundary_values.p

CNMSSMEFTHiggs_semi_analytic_solutions::CNMSSMEFTHiggs_semi_analytic_solutions()
{
   initialize_trial_values();
}

void CNMSSMEFTHiggs_semi_analytic_solutions::initialize_trial_values()
{
   trial_data[0].boundary_values.Azero = 1;
   trial_data[0].boundary_values.m12 = 0;
   trial_data[0].boundary_values.m0Sq = 0;
   trial_data[0].basis_sets.push_back(1);
   trial_data[0].basis_sets.push_back(2);

   trial_data[1].boundary_values.Azero = 0;
   trial_data[1].boundary_values.m12 = 1;
   trial_data[1].boundary_values.m0Sq = 0;
   trial_data[1].basis_sets.push_back(1);
   trial_data[1].basis_sets.push_back(2);

   trial_data[2].boundary_values.Azero = 0;
   trial_data[2].boundary_values.m12 = 0;
   trial_data[2].boundary_values.m0Sq = 1;
   trial_data[2].basis_sets.push_back(2);

   trial_data[3].boundary_values.Azero = 1;
   trial_data[3].boundary_values.m12 = 1.5;
   trial_data[3].boundary_values.m0Sq = 0;
   trial_data[3].basis_sets.push_back(2);

}

void CNMSSMEFTHiggs_semi_analytic_solutions::calculate_trial_data(const CNMSSMEFTHiggs_mass_eigenstates& model)
{
   for (auto& point: trial_data) {
      point.model = run_to_output_scale(model, point.boundary_values);
   }
}

CNMSSMEFTHiggs_mass_eigenstates CNMSSMEFTHiggs_semi_analytic_solutions::run_to_output_scale(
   const CNMSSMEFTHiggs_mass_eigenstates& model, const Boundary_values& values) const
{
   CNMSSMEFTHiggs_mass_eigenstates result(model);
   set_to_boundary_values(result, values);

   result.run_to(output_scale);

   return result;
}

std::map<int,CNMSSMEFTHiggs_semi_analytic_solutions::Data_vector_t> CNMSSMEFTHiggs_semi_analytic_solutions::create_datasets() const
{
   std::map<int,Data_vector_t> datasets;
   for (const auto& point: trial_data) {
      for (auto idx: point.basis_sets) {
         datasets[idx].push_back(&point);
      }
   }

   for (auto& d: datasets)
      d.second.shrink_to_fit();

   return datasets;
}

/**
 * Calculates the coefficients in the semi-analytic solutions
 * by varying the boundary value parameters and running to the
 * output scale.  For each term in the semi-analytic solutions,
 * the boundary value parameters are chosen to make that term
 * non-zero.  The boundary conditions are applied at the
 * defined input scale, and the parameters run to the output scale.
 * The estimates for the soft parameters obtained in this way for
 * each term defines a linear system, which is solved for the
 * coefficients.
 */
void CNMSSMEFTHiggs_semi_analytic_solutions::calculate_coefficients(const CNMSSMEFTHiggs_mass_eigenstates& model_)
{
   CNMSSMEFTHiggs_mass_eigenstates model(model_);
   model.run_to(input_scale);

   calculate_trial_data(model);
   std::map<int,Data_vector_t> datasets(create_datasets());

   auto basis_1 = [] (const Boundary_values& boundary_values) {
      Eigen::Matrix<double,1,2> result;
      const auto Azero = BOUNDARYVALUE(Azero);
      const auto m12 = BOUNDARYVALUE(m12);

      result(0) = Azero;
      result(1) = m12;

      return result;
   };

   auto basis_2 = [] (const Boundary_values& boundary_values) {
      Eigen::Matrix<double,1,4> result;
      const auto m0Sq = BOUNDARYVALUE(m0Sq);
      const auto Azero = BOUNDARYVALUE(Azero);
      const auto m12 = BOUNDARYVALUE(m12);

      result(0) = m0Sq;
      result(1) = Sqr(Azero);
      result(2) = Azero*m12;
      result(3) = Sqr(m12);

      return result;
   };

   
   auto solver_1 = create_solver<Eigen::MatrixXd,2>(datasets[1], basis_1);
   auto solver_2 = create_solver<Eigen::MatrixXd,4>(datasets[2], basis_2);

   calculate_MassB_coefficients(solver_1, datasets[1]);
   calculate_MassG_coefficients(solver_1, datasets[1]);
   calculate_MassWB_coefficients(solver_1, datasets[1]);
   calculate_TYd_coefficients(solver_1, datasets[1]);
   calculate_TYe_coefficients(solver_1, datasets[1]);
   calculate_TYu_coefficients(solver_1, datasets[1]);
   calculate_TKappa_coefficients(solver_1, datasets[1]);
   calculate_TLambdax_coefficients(solver_1, datasets[1]);
   calculate_md2_coefficients(solver_2, datasets[2]);
   calculate_me2_coefficients(solver_2, datasets[2]);
   calculate_mHd2_coefficients(solver_2, datasets[2]);
   calculate_mHu2_coefficients(solver_2, datasets[2]);
   calculate_ml2_coefficients(solver_2, datasets[2]);
   calculate_mq2_coefficients(solver_2, datasets[2]);
   calculate_ms2_coefficients(solver_2, datasets[2]);
   calculate_mu2_coefficients(solver_2, datasets[2]);

}

void CNMSSMEFTHiggs_semi_analytic_solutions::set_to_boundary_values(
   CNMSSMEFTHiggs_mass_eigenstates& model, const Boundary_values& boundary_values) const
{
   const auto LambdaInput = INPUTPARAMETER(LambdaInput);
   const auto Ye = MODELPARAMETER(Ye);
   const auto Yd = MODELPARAMETER(Yd);
   const auto Yu = MODELPARAMETER(Yu);
   const auto Kappa = MODELPARAMETER(Kappa);
   const auto Azero = BOUNDARYVALUE(Azero);
   const auto m12 = BOUNDARYVALUE(m12);
   const auto m0Sq = BOUNDARYVALUE(m0Sq);
   model.set_TYe((Azero*Ye).real());
   model.set_TYd((Azero*Yd).real());
   model.set_TYu((Azero*Yu).real());
   model.set_mq2((m0Sq*UNITMATRIX(3)).real());
   model.set_ml2((m0Sq*UNITMATRIX(3)).real());
   model.set_md2((m0Sq*UNITMATRIX(3)).real());
   model.set_mu2((m0Sq*UNITMATRIX(3)).real());
   model.set_me2((m0Sq*UNITMATRIX(3)).real());
   model.set_mHu2(Re(m0Sq));
   model.set_mHd2(Re(m0Sq));
   model.set_ms2(Re(m0Sq));
   model.set_TKappa(Re(Azero*Kappa));
   model.set_TLambdax(Re(Azero*LambdaInput));
   model.set_MassB(Re(m12));
   model.set_MassWB(Re(m12));
   model.set_MassG(Re(m12));

}

/**
 * Sets the soft parameters in the given model to the values
 * obtained by evaluating the semi-analytic solutions with the
 * current values of the coefficients and the basis parameters
 * appearing in the boundary conditions.
 */
void CNMSSMEFTHiggs_semi_analytic_solutions::evaluate_solutions(
   CNMSSMEFTHiggs_mass_eigenstates& model) const
{
   const auto Azero = INPUTPARAMETER(Azero);
   const auto m12 = INPUTPARAMETER(m12);
   const auto m0Sq = EXTRAPARAMETER(m0Sq);

   model.set_MassB(Re(Azero*MassBCoeff1 + m12*MassBCoeff2));
   model.set_MassG(Re(Azero*MassGCoeff1 + m12*MassGCoeff2));
   model.set_MassWB(Re(Azero*MassWBCoeff1 + m12*MassWBCoeff2));
   model.set_TYd((Azero*TYdCoeff1 + m12*TYdCoeff2).real());
   model.set_TYe((Azero*TYeCoeff1 + m12*TYeCoeff2).real());
   model.set_TYu((Azero*TYuCoeff1 + m12*TYuCoeff2).real());
   model.set_TKappa(Re(Azero*TKappaCoeff1 + m12*TKappaCoeff2));
   model.set_TLambdax(Re(Azero*TLambdaxCoeff1 + m12*TLambdaxCoeff2));
   model.set_md2((m0Sq*md2Coeff1 + Azero*m12*md2Coeff3 + md2Coeff2*Sqr(Azero) +
      md2Coeff4*Sqr(m12)).real());
   model.set_me2((m0Sq*me2Coeff1 + Azero*m12*me2Coeff3 + me2Coeff2*Sqr(Azero) +
      me2Coeff4*Sqr(m12)).real());
   model.set_mHd2(Re(m0Sq*mHd2Coeff1 + Azero*m12*mHd2Coeff3 + mHd2Coeff2*Sqr(Azero
      ) + mHd2Coeff4*Sqr(m12)));
   model.set_mHu2(Re(m0Sq*mHu2Coeff1 + Azero*m12*mHu2Coeff3 + mHu2Coeff2*Sqr(Azero
      ) + mHu2Coeff4*Sqr(m12)));
   model.set_ml2((m0Sq*ml2Coeff1 + Azero*m12*ml2Coeff3 + ml2Coeff2*Sqr(Azero) +
      ml2Coeff4*Sqr(m12)).real());
   model.set_mq2((m0Sq*mq2Coeff1 + Azero*m12*mq2Coeff3 + mq2Coeff2*Sqr(Azero) +
      mq2Coeff4*Sqr(m12)).real());
   model.set_ms2(Re(m0Sq*ms2Coeff1 + Azero*m12*ms2Coeff3 + ms2Coeff2*Sqr(Azero) +
      ms2Coeff4*Sqr(m12)));
   model.set_mu2((m0Sq*mu2Coeff1 + Azero*m12*mu2Coeff3 + mu2Coeff2*Sqr(Azero) +
      mu2Coeff4*Sqr(m12)).real());

}

void CNMSSMEFTHiggs_semi_analytic_solutions::calculate_MassB_coefficients(
   const Eigen::JacobiSVD<Eigen::MatrixXd>& solver, const CNMSSMEFTHiggs_semi_analytic_solutions::Data_vector_t& data)
{
   const std::size_t n = data.size();
   Eigen::VectorXd rhs(n);
   for (std::size_t j = 0; j < n; ++j) {
      rhs(j) = data[j]->model.get_MassB();
   }
   Eigen::VectorXd solution = solver.solve(rhs);
   MassBCoeff1 = solution(0);
   MassBCoeff2 = solution(1);
}

void CNMSSMEFTHiggs_semi_analytic_solutions::calculate_MassG_coefficients(
   const Eigen::JacobiSVD<Eigen::MatrixXd>& solver, const CNMSSMEFTHiggs_semi_analytic_solutions::Data_vector_t& data)
{
   const std::size_t n = data.size();
   Eigen::VectorXd rhs(n);
   for (std::size_t j = 0; j < n; ++j) {
      rhs(j) = data[j]->model.get_MassG();
   }
   Eigen::VectorXd solution = solver.solve(rhs);
   MassGCoeff1 = solution(0);
   MassGCoeff2 = solution(1);
}

void CNMSSMEFTHiggs_semi_analytic_solutions::calculate_MassWB_coefficients(
   const Eigen::JacobiSVD<Eigen::MatrixXd>& solver, const CNMSSMEFTHiggs_semi_analytic_solutions::Data_vector_t& data)
{
   const std::size_t n = data.size();
   Eigen::VectorXd rhs(n);
   for (std::size_t j = 0; j < n; ++j) {
      rhs(j) = data[j]->model.get_MassWB();
   }
   Eigen::VectorXd solution = solver.solve(rhs);
   MassWBCoeff1 = solution(0);
   MassWBCoeff2 = solution(1);
}

void CNMSSMEFTHiggs_semi_analytic_solutions::calculate_TYd_coefficients(
   const Eigen::JacobiSVD<Eigen::MatrixXd>& solver, const CNMSSMEFTHiggs_semi_analytic_solutions::Data_vector_t& data)
{
   const std::size_t n = data.size();
   Eigen::MatrixXd rhs(n, 9);
   for (std::size_t j = 0; j < n; ++j) {
      rhs(j, 0) = data[j]->model.get_TYd(0, 0);
      rhs(j, 1) = data[j]->model.get_TYd(0, 1);
      rhs(j, 2) = data[j]->model.get_TYd(0, 2);
      rhs(j, 3) = data[j]->model.get_TYd(1, 0);
      rhs(j, 4) = data[j]->model.get_TYd(1, 1);
      rhs(j, 5) = data[j]->model.get_TYd(1, 2);
      rhs(j, 6) = data[j]->model.get_TYd(2, 0);
      rhs(j, 7) = data[j]->model.get_TYd(2, 1);
      rhs(j, 8) = data[j]->model.get_TYd(2, 2);
   }
   Eigen::MatrixXd solution = solver.solve(rhs);
   TYdCoeff1(0, 0) = solution(0, 0);
   TYdCoeff1(0, 1) = solution(0, 1);
   TYdCoeff1(0, 2) = solution(0, 2);
   TYdCoeff1(1, 0) = solution(0, 3);
   TYdCoeff1(1, 1) = solution(0, 4);
   TYdCoeff1(1, 2) = solution(0, 5);
   TYdCoeff1(2, 0) = solution(0, 6);
   TYdCoeff1(2, 1) = solution(0, 7);
   TYdCoeff1(2, 2) = solution(0, 8);
   TYdCoeff2(0, 0) = solution(1, 0);
   TYdCoeff2(0, 1) = solution(1, 1);
   TYdCoeff2(0, 2) = solution(1, 2);
   TYdCoeff2(1, 0) = solution(1, 3);
   TYdCoeff2(1, 1) = solution(1, 4);
   TYdCoeff2(1, 2) = solution(1, 5);
   TYdCoeff2(2, 0) = solution(1, 6);
   TYdCoeff2(2, 1) = solution(1, 7);
   TYdCoeff2(2, 2) = solution(1, 8);
}

void CNMSSMEFTHiggs_semi_analytic_solutions::calculate_TYe_coefficients(
   const Eigen::JacobiSVD<Eigen::MatrixXd>& solver, const CNMSSMEFTHiggs_semi_analytic_solutions::Data_vector_t& data)
{
   const std::size_t n = data.size();
   Eigen::MatrixXd rhs(n, 9);
   for (std::size_t j = 0; j < n; ++j) {
      rhs(j, 0) = data[j]->model.get_TYe(0, 0);
      rhs(j, 1) = data[j]->model.get_TYe(0, 1);
      rhs(j, 2) = data[j]->model.get_TYe(0, 2);
      rhs(j, 3) = data[j]->model.get_TYe(1, 0);
      rhs(j, 4) = data[j]->model.get_TYe(1, 1);
      rhs(j, 5) = data[j]->model.get_TYe(1, 2);
      rhs(j, 6) = data[j]->model.get_TYe(2, 0);
      rhs(j, 7) = data[j]->model.get_TYe(2, 1);
      rhs(j, 8) = data[j]->model.get_TYe(2, 2);
   }
   Eigen::MatrixXd solution = solver.solve(rhs);
   TYeCoeff1(0, 0) = solution(0, 0);
   TYeCoeff1(0, 1) = solution(0, 1);
   TYeCoeff1(0, 2) = solution(0, 2);
   TYeCoeff1(1, 0) = solution(0, 3);
   TYeCoeff1(1, 1) = solution(0, 4);
   TYeCoeff1(1, 2) = solution(0, 5);
   TYeCoeff1(2, 0) = solution(0, 6);
   TYeCoeff1(2, 1) = solution(0, 7);
   TYeCoeff1(2, 2) = solution(0, 8);
   TYeCoeff2(0, 0) = solution(1, 0);
   TYeCoeff2(0, 1) = solution(1, 1);
   TYeCoeff2(0, 2) = solution(1, 2);
   TYeCoeff2(1, 0) = solution(1, 3);
   TYeCoeff2(1, 1) = solution(1, 4);
   TYeCoeff2(1, 2) = solution(1, 5);
   TYeCoeff2(2, 0) = solution(1, 6);
   TYeCoeff2(2, 1) = solution(1, 7);
   TYeCoeff2(2, 2) = solution(1, 8);
}

void CNMSSMEFTHiggs_semi_analytic_solutions::calculate_TYu_coefficients(
   const Eigen::JacobiSVD<Eigen::MatrixXd>& solver, const CNMSSMEFTHiggs_semi_analytic_solutions::Data_vector_t& data)
{
   const std::size_t n = data.size();
   Eigen::MatrixXd rhs(n, 9);
   for (std::size_t j = 0; j < n; ++j) {
      rhs(j, 0) = data[j]->model.get_TYu(0, 0);
      rhs(j, 1) = data[j]->model.get_TYu(0, 1);
      rhs(j, 2) = data[j]->model.get_TYu(0, 2);
      rhs(j, 3) = data[j]->model.get_TYu(1, 0);
      rhs(j, 4) = data[j]->model.get_TYu(1, 1);
      rhs(j, 5) = data[j]->model.get_TYu(1, 2);
      rhs(j, 6) = data[j]->model.get_TYu(2, 0);
      rhs(j, 7) = data[j]->model.get_TYu(2, 1);
      rhs(j, 8) = data[j]->model.get_TYu(2, 2);
   }
   Eigen::MatrixXd solution = solver.solve(rhs);
   TYuCoeff1(0, 0) = solution(0, 0);
   TYuCoeff1(0, 1) = solution(0, 1);
   TYuCoeff1(0, 2) = solution(0, 2);
   TYuCoeff1(1, 0) = solution(0, 3);
   TYuCoeff1(1, 1) = solution(0, 4);
   TYuCoeff1(1, 2) = solution(0, 5);
   TYuCoeff1(2, 0) = solution(0, 6);
   TYuCoeff1(2, 1) = solution(0, 7);
   TYuCoeff1(2, 2) = solution(0, 8);
   TYuCoeff2(0, 0) = solution(1, 0);
   TYuCoeff2(0, 1) = solution(1, 1);
   TYuCoeff2(0, 2) = solution(1, 2);
   TYuCoeff2(1, 0) = solution(1, 3);
   TYuCoeff2(1, 1) = solution(1, 4);
   TYuCoeff2(1, 2) = solution(1, 5);
   TYuCoeff2(2, 0) = solution(1, 6);
   TYuCoeff2(2, 1) = solution(1, 7);
   TYuCoeff2(2, 2) = solution(1, 8);
}

void CNMSSMEFTHiggs_semi_analytic_solutions::calculate_TKappa_coefficients(
   const Eigen::JacobiSVD<Eigen::MatrixXd>& solver, const CNMSSMEFTHiggs_semi_analytic_solutions::Data_vector_t& data)
{
   const std::size_t n = data.size();
   Eigen::VectorXd rhs(n);
   for (std::size_t j = 0; j < n; ++j) {
      rhs(j) = data[j]->model.get_TKappa();
   }
   Eigen::VectorXd solution = solver.solve(rhs);
   TKappaCoeff1 = solution(0);
   TKappaCoeff2 = solution(1);
}

void CNMSSMEFTHiggs_semi_analytic_solutions::calculate_TLambdax_coefficients(
   const Eigen::JacobiSVD<Eigen::MatrixXd>& solver, const CNMSSMEFTHiggs_semi_analytic_solutions::Data_vector_t& data)
{
   const std::size_t n = data.size();
   Eigen::VectorXd rhs(n);
   for (std::size_t j = 0; j < n; ++j) {
      rhs(j) = data[j]->model.get_TLambdax();
   }
   Eigen::VectorXd solution = solver.solve(rhs);
   TLambdaxCoeff1 = solution(0);
   TLambdaxCoeff2 = solution(1);
}

void CNMSSMEFTHiggs_semi_analytic_solutions::calculate_md2_coefficients(
   const Eigen::JacobiSVD<Eigen::MatrixXd>& solver, const CNMSSMEFTHiggs_semi_analytic_solutions::Data_vector_t& data)
{
   const std::size_t n = data.size();
   Eigen::MatrixXd rhs(n, 9);
   for (std::size_t j = 0; j < n; ++j) {
      rhs(j, 0) = data[j]->model.get_md2(0, 0);
      rhs(j, 1) = data[j]->model.get_md2(0, 1);
      rhs(j, 2) = data[j]->model.get_md2(0, 2);
      rhs(j, 3) = data[j]->model.get_md2(1, 0);
      rhs(j, 4) = data[j]->model.get_md2(1, 1);
      rhs(j, 5) = data[j]->model.get_md2(1, 2);
      rhs(j, 6) = data[j]->model.get_md2(2, 0);
      rhs(j, 7) = data[j]->model.get_md2(2, 1);
      rhs(j, 8) = data[j]->model.get_md2(2, 2);
   }
   Eigen::MatrixXd solution = solver.solve(rhs);
   md2Coeff1(0, 0) = solution(0, 0);
   md2Coeff1(0, 1) = solution(0, 1);
   md2Coeff1(0, 2) = solution(0, 2);
   md2Coeff1(1, 0) = solution(0, 3);
   md2Coeff1(1, 1) = solution(0, 4);
   md2Coeff1(1, 2) = solution(0, 5);
   md2Coeff1(2, 0) = solution(0, 6);
   md2Coeff1(2, 1) = solution(0, 7);
   md2Coeff1(2, 2) = solution(0, 8);
   md2Coeff2(0, 0) = solution(1, 0);
   md2Coeff2(0, 1) = solution(1, 1);
   md2Coeff2(0, 2) = solution(1, 2);
   md2Coeff2(1, 0) = solution(1, 3);
   md2Coeff2(1, 1) = solution(1, 4);
   md2Coeff2(1, 2) = solution(1, 5);
   md2Coeff2(2, 0) = solution(1, 6);
   md2Coeff2(2, 1) = solution(1, 7);
   md2Coeff2(2, 2) = solution(1, 8);
   md2Coeff3(0, 0) = solution(2, 0);
   md2Coeff3(0, 1) = solution(2, 1);
   md2Coeff3(0, 2) = solution(2, 2);
   md2Coeff3(1, 0) = solution(2, 3);
   md2Coeff3(1, 1) = solution(2, 4);
   md2Coeff3(1, 2) = solution(2, 5);
   md2Coeff3(2, 0) = solution(2, 6);
   md2Coeff3(2, 1) = solution(2, 7);
   md2Coeff3(2, 2) = solution(2, 8);
   md2Coeff4(0, 0) = solution(3, 0);
   md2Coeff4(0, 1) = solution(3, 1);
   md2Coeff4(0, 2) = solution(3, 2);
   md2Coeff4(1, 0) = solution(3, 3);
   md2Coeff4(1, 1) = solution(3, 4);
   md2Coeff4(1, 2) = solution(3, 5);
   md2Coeff4(2, 0) = solution(3, 6);
   md2Coeff4(2, 1) = solution(3, 7);
   md2Coeff4(2, 2) = solution(3, 8);
}

void CNMSSMEFTHiggs_semi_analytic_solutions::calculate_me2_coefficients(
   const Eigen::JacobiSVD<Eigen::MatrixXd>& solver, const CNMSSMEFTHiggs_semi_analytic_solutions::Data_vector_t& data)
{
   const std::size_t n = data.size();
   Eigen::MatrixXd rhs(n, 9);
   for (std::size_t j = 0; j < n; ++j) {
      rhs(j, 0) = data[j]->model.get_me2(0, 0);
      rhs(j, 1) = data[j]->model.get_me2(0, 1);
      rhs(j, 2) = data[j]->model.get_me2(0, 2);
      rhs(j, 3) = data[j]->model.get_me2(1, 0);
      rhs(j, 4) = data[j]->model.get_me2(1, 1);
      rhs(j, 5) = data[j]->model.get_me2(1, 2);
      rhs(j, 6) = data[j]->model.get_me2(2, 0);
      rhs(j, 7) = data[j]->model.get_me2(2, 1);
      rhs(j, 8) = data[j]->model.get_me2(2, 2);
   }
   Eigen::MatrixXd solution = solver.solve(rhs);
   me2Coeff1(0, 0) = solution(0, 0);
   me2Coeff1(0, 1) = solution(0, 1);
   me2Coeff1(0, 2) = solution(0, 2);
   me2Coeff1(1, 0) = solution(0, 3);
   me2Coeff1(1, 1) = solution(0, 4);
   me2Coeff1(1, 2) = solution(0, 5);
   me2Coeff1(2, 0) = solution(0, 6);
   me2Coeff1(2, 1) = solution(0, 7);
   me2Coeff1(2, 2) = solution(0, 8);
   me2Coeff2(0, 0) = solution(1, 0);
   me2Coeff2(0, 1) = solution(1, 1);
   me2Coeff2(0, 2) = solution(1, 2);
   me2Coeff2(1, 0) = solution(1, 3);
   me2Coeff2(1, 1) = solution(1, 4);
   me2Coeff2(1, 2) = solution(1, 5);
   me2Coeff2(2, 0) = solution(1, 6);
   me2Coeff2(2, 1) = solution(1, 7);
   me2Coeff2(2, 2) = solution(1, 8);
   me2Coeff3(0, 0) = solution(2, 0);
   me2Coeff3(0, 1) = solution(2, 1);
   me2Coeff3(0, 2) = solution(2, 2);
   me2Coeff3(1, 0) = solution(2, 3);
   me2Coeff3(1, 1) = solution(2, 4);
   me2Coeff3(1, 2) = solution(2, 5);
   me2Coeff3(2, 0) = solution(2, 6);
   me2Coeff3(2, 1) = solution(2, 7);
   me2Coeff3(2, 2) = solution(2, 8);
   me2Coeff4(0, 0) = solution(3, 0);
   me2Coeff4(0, 1) = solution(3, 1);
   me2Coeff4(0, 2) = solution(3, 2);
   me2Coeff4(1, 0) = solution(3, 3);
   me2Coeff4(1, 1) = solution(3, 4);
   me2Coeff4(1, 2) = solution(3, 5);
   me2Coeff4(2, 0) = solution(3, 6);
   me2Coeff4(2, 1) = solution(3, 7);
   me2Coeff4(2, 2) = solution(3, 8);
}

void CNMSSMEFTHiggs_semi_analytic_solutions::calculate_mHd2_coefficients(
   const Eigen::JacobiSVD<Eigen::MatrixXd>& solver, const CNMSSMEFTHiggs_semi_analytic_solutions::Data_vector_t& data)
{
   const std::size_t n = data.size();
   Eigen::VectorXd rhs(n);
   for (std::size_t j = 0; j < n; ++j) {
      rhs(j) = data[j]->model.get_mHd2();
   }
   Eigen::VectorXd solution = solver.solve(rhs);
   mHd2Coeff1 = solution(0);
   mHd2Coeff2 = solution(1);
   mHd2Coeff3 = solution(2);
   mHd2Coeff4 = solution(3);
}

void CNMSSMEFTHiggs_semi_analytic_solutions::calculate_mHu2_coefficients(
   const Eigen::JacobiSVD<Eigen::MatrixXd>& solver, const CNMSSMEFTHiggs_semi_analytic_solutions::Data_vector_t& data)
{
   const std::size_t n = data.size();
   Eigen::VectorXd rhs(n);
   for (std::size_t j = 0; j < n; ++j) {
      rhs(j) = data[j]->model.get_mHu2();
   }
   Eigen::VectorXd solution = solver.solve(rhs);
   mHu2Coeff1 = solution(0);
   mHu2Coeff2 = solution(1);
   mHu2Coeff3 = solution(2);
   mHu2Coeff4 = solution(3);
}

void CNMSSMEFTHiggs_semi_analytic_solutions::calculate_ml2_coefficients(
   const Eigen::JacobiSVD<Eigen::MatrixXd>& solver, const CNMSSMEFTHiggs_semi_analytic_solutions::Data_vector_t& data)
{
   const std::size_t n = data.size();
   Eigen::MatrixXd rhs(n, 9);
   for (std::size_t j = 0; j < n; ++j) {
      rhs(j, 0) = data[j]->model.get_ml2(0, 0);
      rhs(j, 1) = data[j]->model.get_ml2(0, 1);
      rhs(j, 2) = data[j]->model.get_ml2(0, 2);
      rhs(j, 3) = data[j]->model.get_ml2(1, 0);
      rhs(j, 4) = data[j]->model.get_ml2(1, 1);
      rhs(j, 5) = data[j]->model.get_ml2(1, 2);
      rhs(j, 6) = data[j]->model.get_ml2(2, 0);
      rhs(j, 7) = data[j]->model.get_ml2(2, 1);
      rhs(j, 8) = data[j]->model.get_ml2(2, 2);
   }
   Eigen::MatrixXd solution = solver.solve(rhs);
   ml2Coeff1(0, 0) = solution(0, 0);
   ml2Coeff1(0, 1) = solution(0, 1);
   ml2Coeff1(0, 2) = solution(0, 2);
   ml2Coeff1(1, 0) = solution(0, 3);
   ml2Coeff1(1, 1) = solution(0, 4);
   ml2Coeff1(1, 2) = solution(0, 5);
   ml2Coeff1(2, 0) = solution(0, 6);
   ml2Coeff1(2, 1) = solution(0, 7);
   ml2Coeff1(2, 2) = solution(0, 8);
   ml2Coeff2(0, 0) = solution(1, 0);
   ml2Coeff2(0, 1) = solution(1, 1);
   ml2Coeff2(0, 2) = solution(1, 2);
   ml2Coeff2(1, 0) = solution(1, 3);
   ml2Coeff2(1, 1) = solution(1, 4);
   ml2Coeff2(1, 2) = solution(1, 5);
   ml2Coeff2(2, 0) = solution(1, 6);
   ml2Coeff2(2, 1) = solution(1, 7);
   ml2Coeff2(2, 2) = solution(1, 8);
   ml2Coeff3(0, 0) = solution(2, 0);
   ml2Coeff3(0, 1) = solution(2, 1);
   ml2Coeff3(0, 2) = solution(2, 2);
   ml2Coeff3(1, 0) = solution(2, 3);
   ml2Coeff3(1, 1) = solution(2, 4);
   ml2Coeff3(1, 2) = solution(2, 5);
   ml2Coeff3(2, 0) = solution(2, 6);
   ml2Coeff3(2, 1) = solution(2, 7);
   ml2Coeff3(2, 2) = solution(2, 8);
   ml2Coeff4(0, 0) = solution(3, 0);
   ml2Coeff4(0, 1) = solution(3, 1);
   ml2Coeff4(0, 2) = solution(3, 2);
   ml2Coeff4(1, 0) = solution(3, 3);
   ml2Coeff4(1, 1) = solution(3, 4);
   ml2Coeff4(1, 2) = solution(3, 5);
   ml2Coeff4(2, 0) = solution(3, 6);
   ml2Coeff4(2, 1) = solution(3, 7);
   ml2Coeff4(2, 2) = solution(3, 8);
}

void CNMSSMEFTHiggs_semi_analytic_solutions::calculate_mq2_coefficients(
   const Eigen::JacobiSVD<Eigen::MatrixXd>& solver, const CNMSSMEFTHiggs_semi_analytic_solutions::Data_vector_t& data)
{
   const std::size_t n = data.size();
   Eigen::MatrixXd rhs(n, 9);
   for (std::size_t j = 0; j < n; ++j) {
      rhs(j, 0) = data[j]->model.get_mq2(0, 0);
      rhs(j, 1) = data[j]->model.get_mq2(0, 1);
      rhs(j, 2) = data[j]->model.get_mq2(0, 2);
      rhs(j, 3) = data[j]->model.get_mq2(1, 0);
      rhs(j, 4) = data[j]->model.get_mq2(1, 1);
      rhs(j, 5) = data[j]->model.get_mq2(1, 2);
      rhs(j, 6) = data[j]->model.get_mq2(2, 0);
      rhs(j, 7) = data[j]->model.get_mq2(2, 1);
      rhs(j, 8) = data[j]->model.get_mq2(2, 2);
   }
   Eigen::MatrixXd solution = solver.solve(rhs);
   mq2Coeff1(0, 0) = solution(0, 0);
   mq2Coeff1(0, 1) = solution(0, 1);
   mq2Coeff1(0, 2) = solution(0, 2);
   mq2Coeff1(1, 0) = solution(0, 3);
   mq2Coeff1(1, 1) = solution(0, 4);
   mq2Coeff1(1, 2) = solution(0, 5);
   mq2Coeff1(2, 0) = solution(0, 6);
   mq2Coeff1(2, 1) = solution(0, 7);
   mq2Coeff1(2, 2) = solution(0, 8);
   mq2Coeff2(0, 0) = solution(1, 0);
   mq2Coeff2(0, 1) = solution(1, 1);
   mq2Coeff2(0, 2) = solution(1, 2);
   mq2Coeff2(1, 0) = solution(1, 3);
   mq2Coeff2(1, 1) = solution(1, 4);
   mq2Coeff2(1, 2) = solution(1, 5);
   mq2Coeff2(2, 0) = solution(1, 6);
   mq2Coeff2(2, 1) = solution(1, 7);
   mq2Coeff2(2, 2) = solution(1, 8);
   mq2Coeff3(0, 0) = solution(2, 0);
   mq2Coeff3(0, 1) = solution(2, 1);
   mq2Coeff3(0, 2) = solution(2, 2);
   mq2Coeff3(1, 0) = solution(2, 3);
   mq2Coeff3(1, 1) = solution(2, 4);
   mq2Coeff3(1, 2) = solution(2, 5);
   mq2Coeff3(2, 0) = solution(2, 6);
   mq2Coeff3(2, 1) = solution(2, 7);
   mq2Coeff3(2, 2) = solution(2, 8);
   mq2Coeff4(0, 0) = solution(3, 0);
   mq2Coeff4(0, 1) = solution(3, 1);
   mq2Coeff4(0, 2) = solution(3, 2);
   mq2Coeff4(1, 0) = solution(3, 3);
   mq2Coeff4(1, 1) = solution(3, 4);
   mq2Coeff4(1, 2) = solution(3, 5);
   mq2Coeff4(2, 0) = solution(3, 6);
   mq2Coeff4(2, 1) = solution(3, 7);
   mq2Coeff4(2, 2) = solution(3, 8);
}

void CNMSSMEFTHiggs_semi_analytic_solutions::calculate_ms2_coefficients(
   const Eigen::JacobiSVD<Eigen::MatrixXd>& solver, const CNMSSMEFTHiggs_semi_analytic_solutions::Data_vector_t& data)
{
   const std::size_t n = data.size();
   Eigen::VectorXd rhs(n);
   for (std::size_t j = 0; j < n; ++j) {
      rhs(j) = data[j]->model.get_ms2();
   }
   Eigen::VectorXd solution = solver.solve(rhs);
   ms2Coeff1 = solution(0);
   ms2Coeff2 = solution(1);
   ms2Coeff3 = solution(2);
   ms2Coeff4 = solution(3);
}

void CNMSSMEFTHiggs_semi_analytic_solutions::calculate_mu2_coefficients(
   const Eigen::JacobiSVD<Eigen::MatrixXd>& solver, const CNMSSMEFTHiggs_semi_analytic_solutions::Data_vector_t& data)
{
   const std::size_t n = data.size();
   Eigen::MatrixXd rhs(n, 9);
   for (std::size_t j = 0; j < n; ++j) {
      rhs(j, 0) = data[j]->model.get_mu2(0, 0);
      rhs(j, 1) = data[j]->model.get_mu2(0, 1);
      rhs(j, 2) = data[j]->model.get_mu2(0, 2);
      rhs(j, 3) = data[j]->model.get_mu2(1, 0);
      rhs(j, 4) = data[j]->model.get_mu2(1, 1);
      rhs(j, 5) = data[j]->model.get_mu2(1, 2);
      rhs(j, 6) = data[j]->model.get_mu2(2, 0);
      rhs(j, 7) = data[j]->model.get_mu2(2, 1);
      rhs(j, 8) = data[j]->model.get_mu2(2, 2);
   }
   Eigen::MatrixXd solution = solver.solve(rhs);
   mu2Coeff1(0, 0) = solution(0, 0);
   mu2Coeff1(0, 1) = solution(0, 1);
   mu2Coeff1(0, 2) = solution(0, 2);
   mu2Coeff1(1, 0) = solution(0, 3);
   mu2Coeff1(1, 1) = solution(0, 4);
   mu2Coeff1(1, 2) = solution(0, 5);
   mu2Coeff1(2, 0) = solution(0, 6);
   mu2Coeff1(2, 1) = solution(0, 7);
   mu2Coeff1(2, 2) = solution(0, 8);
   mu2Coeff2(0, 0) = solution(1, 0);
   mu2Coeff2(0, 1) = solution(1, 1);
   mu2Coeff2(0, 2) = solution(1, 2);
   mu2Coeff2(1, 0) = solution(1, 3);
   mu2Coeff2(1, 1) = solution(1, 4);
   mu2Coeff2(1, 2) = solution(1, 5);
   mu2Coeff2(2, 0) = solution(1, 6);
   mu2Coeff2(2, 1) = solution(1, 7);
   mu2Coeff2(2, 2) = solution(1, 8);
   mu2Coeff3(0, 0) = solution(2, 0);
   mu2Coeff3(0, 1) = solution(2, 1);
   mu2Coeff3(0, 2) = solution(2, 2);
   mu2Coeff3(1, 0) = solution(2, 3);
   mu2Coeff3(1, 1) = solution(2, 4);
   mu2Coeff3(1, 2) = solution(2, 5);
   mu2Coeff3(2, 0) = solution(2, 6);
   mu2Coeff3(2, 1) = solution(2, 7);
   mu2Coeff3(2, 2) = solution(2, 8);
   mu2Coeff4(0, 0) = solution(3, 0);
   mu2Coeff4(0, 1) = solution(3, 1);
   mu2Coeff4(0, 2) = solution(3, 2);
   mu2Coeff4(1, 0) = solution(3, 3);
   mu2Coeff4(1, 1) = solution(3, 4);
   mu2Coeff4(1, 2) = solution(3, 5);
   mu2Coeff4(2, 0) = solution(3, 6);
   mu2Coeff4(2, 1) = solution(3, 7);
   mu2Coeff4(2, 2) = solution(3, 8);
}

} // namespace flexiblesusy
