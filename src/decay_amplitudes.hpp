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

#ifndef DECAY_AMPLITUDES_H
#define DECAY_AMPLITUDES_H

#include <complex>

namespace flexiblesusy {

/**
 * @class Decay_amplitude_SSS
 * @brief generic amplitude for the decay of a scalar into two scalars
 */
struct Decay_amplitude_SSS {
   double m_decay{0.};
   double m_out_1{0.};
   double m_out_2{0.};
   std::complex<double> matrix_element{};

   double square() const;
};

/**
 * @class Decay_amplitude_SSV
 * @brief generic amplitude for the decay of a scalar into a scalar and vector
 */
struct Decay_amplitude_SSV {
   double m_decay{0.};
   double m_scalar{0.};
   double m_vector{0.};
   std::complex<double> matrix_element{};

   double square() const;
};

/**
 * @class Decay_amplitude_SVV
 * @brief generic amplitude for the decay of a scalar into two vectors
 */
struct Decay_amplitude_SVV {
   double m_decay{0.};
   double m_out_1{0.};
   double m_out_2{0.};
   std::complex<double> M1{};
   std::complex<double> M2{};

   double square() const;
};

/**
 * @class Decay_amplitude_SFF
 * @brief generic amplitude for the decay of a scalar into two fermions
 */
struct Decay_amplitude_SFF {
   double m_decay{0.};
   double m_out_1{0.};
   double m_out_2{0.};
   std::complex<double> matrix_element_left{};
   std::complex<double> matrix_element_right{};

   double square() const;
};

/**
 * @class Decay_amplitude_FFS
 * @brief generic amplitude for the decay of a fermion into a fermion and scalar
 */
struct Decay_amplitude_FFS {
   double m_decay{0.};
   double m_fermion{0.};
   double m_scalar{0.};
   std::complex<double> matrix_element_left{};
   std::complex<double> matrix_element_right{};

   double square() const;
};

/**
 * @class Decay_amplitude_FFV
 * @brief generic amplitude for the decay of a fermion into a fermion and vector
 */
struct Decay_amplitude_FFV {
   double m_decay{0.};
   double m_fermion{0.};
   double m_vector{0.};
   std::complex<double> matrix_element_gam_1{};
   std::complex<double> matrix_element_gam_2{};
   std::complex<double> matrix_element_p_1{};
   std::complex<double> matrix_element_p_2{};

   double square() const;
};

template <typename Amplitude>
double square_amplitude(const Amplitude& a)
{
   return a.square();
}

} // namespace flexiblesusy

#endif
