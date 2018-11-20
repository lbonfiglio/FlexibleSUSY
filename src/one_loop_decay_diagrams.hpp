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

#ifndef ONE_LOOP_DECAY_DIAGRAMS_H
#define ONE_LOOP_DECAY_DIAGRAMS_H

#include "decay_amplitudes.hpp"

#include <complex>

namespace flexiblesusy {

Decay_amplitude_SSS calculate_diagram_SSS_t1g1n1_FFF(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, const std::complex<double>&, const
   std::complex<double>&, const std::complex<double>&, double);

Decay_amplitude_SSS calculate_diagram_SSS_t1g2n2_SSS(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSS calculate_diagram_SSS_t1g3n3_UUU(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSS calculate_diagram_SSS_t1g4n4_SSV(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSS calculate_diagram_SSS_t1g5n5_SVS(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSS calculate_diagram_SSS_t1g6n6_VSS(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSS calculate_diagram_SSS_t1g7n7_SVV(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSS calculate_diagram_SSS_t1g8n8_VSV(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSS calculate_diagram_SSS_t1g9n9_VVS(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSS calculate_diagram_SSS_t1g10n10_VVV(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSS calculate_diagram_SSS_t2g1n11_SS(double, double, double,
   double, double, const std::complex<double>&, const std::complex<double>&,
   double);

Decay_amplitude_SSS calculate_diagram_SSS_t2g2n12_VV(double, double, double,
   double, double, const std::complex<double>&, const std::complex<double>&,
   double, double);

Decay_amplitude_SSS calculate_diagram_SSS_t3g1n13_SS(double, double, double,
   double, double, const std::complex<double>&, const std::complex<double>&,
   double);

Decay_amplitude_SSS calculate_diagram_SSS_t3g2n14_VV(double, double, double,
   double, double, const std::complex<double>&, const std::complex<double>&,
   double, double);

Decay_amplitude_SSS calculate_diagram_SSS_t4g1n15_SS(double, double, double,
   double, double, const std::complex<double>&, const std::complex<double>&,
   double);

Decay_amplitude_SSS calculate_diagram_SSS_t4g2n16_VV(double, double, double,
   double, double, const std::complex<double>&, const std::complex<double>&,
   double, double);

Decay_amplitude_SVV calculate_diagram_SVV_t1g1n1_FFF(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, const std::complex<double>&, const
   std::complex<double>&, const std::complex<double>&, double);

Decay_amplitude_SVV calculate_diagram_SVV_t1g2n2_SSS(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SVV calculate_diagram_SVV_t1g3n3_UUU(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SVV calculate_diagram_SVV_t1g4n4_SSV(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SVV calculate_diagram_SVV_t1g5n5_SVS(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SVV calculate_diagram_SVV_t1g6n6_VSS(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SVV calculate_diagram_SVV_t1g7n7_SVV(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SVV calculate_diagram_SVV_t1g8n8_VSV(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SVV calculate_diagram_SVV_t1g9n9_VVS(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SVV calculate_diagram_SVV_t1g10n10_VVV(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double, double);

Decay_amplitude_SVV calculate_diagram_SVV_t2g1n11_SS(double, double, double,
   double, double, const std::complex<double>&, const std::complex<double>&,
   double);

Decay_amplitude_SVV calculate_diagram_SVV_t2g2n12_VV(double, double, double,
   double, double, const std::complex<double>&, const std::complex<double>&,
   const std::complex<double>&, const std::complex<double>&, double, double);

Decay_amplitude_SVV calculate_diagram_SVV_t3g1n13_SV(double, double, double,
   double, double, const std::complex<double>&, const std::complex<double>&,
   double);

Decay_amplitude_SVV calculate_diagram_SVV_t4g1n14_SV(double, double, double,
   double, double, const std::complex<double>&, const std::complex<double>&,
   double);

Decay_amplitude_SSV calculate_diagram_SSV_t1g1n1_FFF(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, const std::complex<double>&, const
   std::complex<double>&, const std::complex<double>&, double);

Decay_amplitude_SSV calculate_diagram_SSV_t1g2n2_SSS(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSV calculate_diagram_SSV_t1g3n3_UUU(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSV calculate_diagram_SSV_t1g4n4_SSV(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSV calculate_diagram_SSV_t1g5n5_SVS(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSV calculate_diagram_SSV_t1g6n6_VSS(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSV calculate_diagram_SSV_t1g7n7_SVV(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSV calculate_diagram_SSV_t1g8n8_VSV(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSV calculate_diagram_SSV_t1g9n9_VVS(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSV calculate_diagram_SSV_t1g10n10_VVV(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSV calculate_diagram_SSV_t2g1n11_SV(double, double, double,
   double, double, const std::complex<double>&, const std::complex<double>&,
   double);

Decay_amplitude_SSV calculate_diagram_SSV_t3g1n12_SV(double, double, double,
   double, double, const std::complex<double>&, const std::complex<double>&,
   double);

Decay_amplitude_SFF calculate_diagram_SFF_t1g1n1_FFS(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, const std::complex<double>&, const
   std::complex<double>&, const std::complex<double>&, double);

Decay_amplitude_SFF calculate_diagram_SFF_t1g2n2_SSF(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, const std::complex<double>&, const
   std::complex<double>&, double);

Decay_amplitude_SFF calculate_diagram_SFF_t1g3n3_FFV(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, const std::complex<double>&, const
   std::complex<double>&, const std::complex<double>&, double, double);

Decay_amplitude_SFF calculate_diagram_SFF_t1g4n4_SVF(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, const std::complex<double>&, const
   std::complex<double>&, double);

Decay_amplitude_SFF calculate_diagram_SFF_t1g5n5_VSF(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, const std::complex<double>&, const
   std::complex<double>&, double);

Decay_amplitude_SFF calculate_diagram_SFF_t1g6n6_VVF(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, const std::complex<double>&, const
   std::complex<double>&, double);

} // namespace flexiblesusy

#endif
