#ifndef TEST_MSSM_H
#define TEST_MSSM_H

#include "MSSM_mass_eigenstates.hpp"

#include "ew_input.hpp"
#include "wrappers.hpp"

void setup_MSSM_const(flexiblesusy::MSSM_mass_eigenstates& m,
                      const flexiblesusy::MSSM_input_parameters& input)
{
   using namespace flexiblesusy;

   const double ALPHASMZ = 0.1176;
   const double ALPHAMZ = 1.0 / 127.918;
   const double sinthWsq = 0.23122;
   const double alpha1 = 5.0 * ALPHAMZ / (3.0 * (1.0 - sinthWsq));
   const double alpha2 = ALPHAMZ / sinthWsq;
   const double g1 = Sqrt(4.0 * Pi * alpha1);
   const double g2 = Sqrt(4.0 * Pi * alpha2);
   const double g3 = Sqrt(4.0 * Pi * ALPHASMZ);

   const double root2 = sqrt(2.0);
   const double vev = 246.0;
   const double tanB = 20.;
   const double vu = vev*std::sin(std::atan(tanB));
   const double vd = vev*std::cos(std::atan(tanB));
   const double scale = Electroweak_constants::MZ;

   Eigen::Matrix<double,3,3> Yu(Eigen::Matrix<double,3,3>::Zero()),
      Yd(Eigen::Matrix<double,3,3>::Zero()),
      Ye(Eigen::Matrix<double,3,3>::Zero());
   Yu(0,0) = 0.001;
   Yu(1,1) = 0.001;
   Yu(2,2) = 165.0   * root2 / vu;
   Yd(0,0) = 0.001;
   Yd(1,1) = 0.001;
   Yd(2,2) = 2.9     * root2 / vd;
   Ye(0,0) = 5.11e-4 * root2 / vd;
   Ye(1,1) = 0.105  * root2 / vd;
   Ye(2,2) = 1.77699 * root2 / vd;

   m.set_scale(scale);
   m.set_loops(1);
   m.set_g1(g1);
   m.set_g2(g2);
   m.set_g3(g3);
   m.set_Yu(Yu);
   m.set_Yd(Yd);
   m.set_Ye(Ye);
   m.set_vu(vu);
   m.set_vd(vd);
   m.set_mHd2(250000.);
   m.set_mHu2(-200000.);
   m.set_MassB(200.);
   m.set_MassWB(400.);
   m.set_MassG(800.);

   Eigen::Matrix<double,3,3> md2(Eigen::Matrix<double,3,3>::Zero()),
mu2(Eigen::Matrix<double,3,3>::Zero()),
mq2(Eigen::Matrix<double,3,3>::Zero()),
me2(Eigen::Matrix<double,3,3>::Zero()),
ml2(Eigen::Matrix<double,3,3>::Zero());

   mu2(0,0) = 1e+6;
   mu2(1,1) = 1e+6;
   mu2(2,2) = 1e+6;
   mq2(0,0) = 1e+6;
   mq2(1,1) = 1e+6;
   mq2(2,2) = 1e+6;
   md2(0,0) = 1e+6;
   md2(1,1) = 1e+6;
   md2(2,2) = 1e+6;
   
   m.set_md2(md2);
   m.set_mq2(mq2);
   m.set_mu2(mu2);

   ml2(0,0) = 1e+6;
   ml2(1,1) = 1e+6;
   ml2(2,2) = 1e+6;
   me2(0,0) = 1e+6;
   me2(1,1) = 1e+6;
   me2(2,2) = 1e+6;
   m.set_ml2(ml2);
   m.set_me2(me2);

   m.solve_ewsb_tree_level();
}

#endif