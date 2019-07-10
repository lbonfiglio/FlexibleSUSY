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

// File generated at Wed 10 Jul 2019 12:34:54

#ifndef NMSSM_EFFECTIVE_COUPLINGS_H
#define NMSSM_EFFECTIVE_COUPLINGS_H

#include "NMSSM_mass_eigenstates.hpp"
#include "lowe.h"
#include "physical_input.hpp"
#include "standard_model.hpp"

#include <complex>
#include <Eigen/Core>

namespace flexiblesusy {

namespace standard_model {
class Standard_model;
}

class NMSSM_effective_couplings {
public:
   NMSSM_effective_couplings(const NMSSM_mass_eigenstates&,
                                   const softsusy::QedQcd&,
                                   const Physical_input&);

   void do_run_couplings(bool flag) { rg_improve = flag; }
   bool do_run_couplings() const { return rg_improve; }
   void do_include_qcd_corrections(bool flag) { include_qcd_corrections = flag; }
   bool do_include_qcd_corrections() const { return include_qcd_corrections; }
   void set_physical_inputs(const Physical_input& inputs_) { physical_input = inputs_; }
   void set_low_energy_data(const softsusy::QedQcd& qedqcd_) { qedqcd = qedqcd_; }
   void set_model(const NMSSM_mass_eigenstates& model_);

   double get_hhVPVP_partial_width(int gO1) const;
   double get_hhVGVG_partial_width(int gO1) const;
   double get_AhVPVP_partial_width(int gO1) const;
   double get_AhVGVG_partial_width(int gO1) const;
   std::complex<double> get_eff_CphhVPVP(int gO1) const { return eff_CphhVPVP(gO1); }
   std::complex<double> get_eff_CphhVGVG(int gO1) const { return eff_CphhVGVG(gO1); }
   std::complex<double> get_eff_CpAhVPVP(int gO1) const { return eff_CpAhVPVP(gO1); }
   std::complex<double> get_eff_CpAhVGVG(int gO1) const { return eff_CpAhVGVG(gO1); }

   void calculate_effective_couplings();

   std::complex<double> CphhconjVWmVWm(int gI2) const;
   std::complex<double> CpbarFeFeAhPL(int gO2, int gI1, int gI2) const;
   std::complex<double> CpbarFeFehhPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarFdFdAhPL(int gO2, int gI1, int gI2) const;
   std::complex<double> CpbarFdFdhhPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarFuFuAhPL(int gO2, int gI1, int gI2) const;
   std::complex<double> CpbarFuFuhhPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CphhSdconjSd(int gt1, int gt2, int gt3) const;
   std::complex<double> CphhSuconjSu(int gt1, int gt2, int gt3) const;
   std::complex<double> CphhSeconjSe(int gt1, int gt2, int gt3) const;
   std::complex<double> CphhHpmconjHpm(int gt1, int gt2, int gt3) const;
   std::complex<double> CpbarChaChahhPL(int gt3, int gt1, int gt2) const;
   std::complex<double> CpAhSdconjSd(int gt1, int gt2, int gt3) const;
   std::complex<double> CpAhSuconjSu(int gt1, int gt2, int gt3) const;
   std::complex<double> CpAhSeconjSe(int gt1, int gt2, int gt3) const;
   std::complex<double> CpAhHpmconjHpm(int gt1, int gt2, int gt3) const;
   std::complex<double> CpbarChaChaAhPL(int gt3, int gt2, int gt1) const;
   void calculate_eff_CphhVPVP(int gO1);
   void calculate_eff_CphhVGVG(int gO1);
   void calculate_eff_CpAhVPVP(int gO1);
   void calculate_eff_CpAhVGVG(int gO1);

private:
   NMSSM_mass_eigenstates model;
   softsusy::QedQcd qedqcd;
   Physical_input physical_input;
   bool rg_improve;
   bool include_qcd_corrections;

   void copy_mixing_matrices_from_model();

   standard_model::Standard_model initialise_SM() const;
   void run_SM_strong_coupling_to(standard_model::Standard_model, double m);

   // higher order corrections to the amplitudes for
   // effective coupling to photons
   std::complex<double> scalar_scalar_qcd_factor(double, double) const;
   std::complex<double> scalar_fermion_qcd_factor(double, double) const;
   std::complex<double> pseudoscalar_fermion_qcd_factor(double, double) const;

   // higher order corrections to the leading order
   // effective couplings to gluons
   double number_of_active_flavours(double) const;
   double scalar_scaling_factor(double) const;
   double pseudoscalar_scaling_factor(double) const;

   Eigen::Matrix<double,6,6> ZD{Eigen::Matrix<double,6,6>::Zero()};
   Eigen::Matrix<double,3,3> ZV{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,6,6> ZU{Eigen::Matrix<double,6,6>::Zero()};
   Eigen::Matrix<double,6,6> ZE{Eigen::Matrix<double,6,6>::Zero()};
   Eigen::Matrix<double,3,3> ZH{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> ZA{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,2,2> ZP{Eigen::Matrix<double,2,2>::Zero()};
   Eigen::Matrix<std::complex<double>,5,5> ZN{Eigen::Matrix<std::complex<double>,5,5>::Zero()};
   Eigen::Matrix<std::complex<double>,2,2> UM{Eigen::Matrix<std::complex<double>,2,2>::Zero()};
   Eigen::Matrix<std::complex<double>,2,2> UP{Eigen::Matrix<std::complex<double>,2,2>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> ZEL{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> ZER{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> ZDL{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> ZDR{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> ZUL{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> ZUR{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<double,2,2> ZZ{Eigen::Matrix<double,2,2>::Zero()};

   Eigen::Array<std::complex<double>,3,1> eff_CphhVPVP;
   Eigen::Array<std::complex<double>,3,1> eff_CphhVGVG;
   Eigen::Array<std::complex<double>,3,1> eff_CpAhVPVP;
   Eigen::Array<std::complex<double>,3,1> eff_CpAhVGVG;

};

} // namespace flexiblesusy

#endif
