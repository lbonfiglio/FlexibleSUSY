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

/**
 * @todo implement formulas for most general forms of amplitudes
 */

#include "decay_amplitudes.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

double Decay_amplitude_SSS::square() const
{
   return AbsSqr(matrix_element);
}

// @todo handle massless vectors safely
double Decay_amplitude_SSV::square() const
{
   const double m_in_sq = Sqr(m_decay);
   const double m_s_sq = Sqr(m_scalar);
   const double m_v_sq = Sqr(m_vector);

   return 0.25 * (Sqr(m_in_sq) + Sqr(m_s_sq - m_v_sq)
                  - 2. * m_in_sq * (m_s_sq + m_v_sq))
      * AbsSqr(matrix_element) / m_v_sq;
}

// @todo implement fully general decomposition
// @todo handle massless vectors safely
double Decay_amplitude_SVV::square() const
{
   const double m_in_sq = Sqr(m_decay);
   const double m_in_4 = Power4(m_decay);
   const double m_1_sq = Sqr(m_out_1);
   const double m_1_4 = Power4(m_out_1);
   const double m_2_sq = Sqr(m_out_2);
   const double m_2_4 = Power4(m_out_2);

   const double m11 = 0.5 * (m_in_4 + m_1_4 + m_2_4 + 10. * m_1_sq * m_2_sq
                             - 2. * m_in_sq * (m_1_sq + m_2_sq)) * AbsSqr(M1)
      / (m_1_sq * m_2_sq);

   const double m22 = 0.125 * Sqr(m_in_4 + Sqr(m_1_sq - m_2_sq)
                                  - 2. * m_in_sq * (m_1_sq + m_2_sq))
      * AbsSqr(M2) / (m_1_sq * m_2_sq);

   const double m12 = 0.5 * (m_in_sq * m_in_4 - 3. * m_in_4 * (m_1_sq + m_2_sq)
                              - Sqr(m_1_sq - m_2_sq) * (m_1_sq + m_2_sq)
                              + m_1_sq * (3. * m_1_4 + 2. * m_1_sq * m_2_sq
                                          + 3. * m_2_4))
      * Re(M1 * Conj(M2)) / (m_1_sq * m_2_sq);

   return m11 + m22 + m12;
}

double Decay_amplitude_SFF::square() const
{
   const double m_in_sq = Sqr(m_decay);
   const double m_1_sq = Sqr(m_out_1);
   const double m_2_sq = Sqr(m_out_2);

   return (m_in_sq - m_1_sq - m_2_sq) *
      (AbsSqr(matrix_element_left) + AbsSqr(matrix_element_right))
      - 4. * m_out_1 * m_out_2 * Re(
         matrix_element_left * Conj(matrix_element_right));
}

double Decay_amplitude_FFS::square() const
{
   const double m_in_sq = Sqr(m_decay);
   const double m_f_sq = Sqr(m_fermion);
   const double m_s_sq = Sqr(m_scalar);

   return 0.5 * (m_in_sq + m_f_sq - m_s_sq) *
      (AbsSqr(matrix_element_left) + AbsSqr(matrix_element_right))
      + 2. * m_fermion * m_scalar * Re(
         matrix_element_left * Conj(matrix_element_right));
}

// @todo handle massless vectors safely
double Decay_amplitude_FFV::square() const
{
   const double m_in_sq = Sqr(m_decay);
   const double m_in_4 = Power4(m_decay);
   const double m_f_sq = Sqr(m_fermion);
   const double m_f_4 = Power4(m_fermion);
   const double m_v_sq = Sqr(m_vector);
   const double m_v_4 = Power4(m_vector);

   const double c1 = 0.5 * (m_in_4 + m_f_4 + m_f_sq * m_v_sq - 2. * m_v_4
                            + m_in_sq * (m_v_sq - 2. * m_f_sq)) / m_v_sq;

   const double c2 = 0.125 * ((m_in_sq + m_f_sq - m_v_sq) * (
                                 m_in_4  + Sqr(m_f_sq - m_v_sq)
                                 - 2. * m_in_sq * (m_f_sq + m_v_sq))) / m_v_sq;

   const double c3 = -3. * m_decay * m_fermion;

   const double c4 = -0.25 * (m_fermion * (m_in_4 + Sqr(m_f_sq - m_v_sq)
                                           - 2. * m_in_sq * (m_f_sq + m_v_sq)))
      / m_v_sq;

   const double c5 = -0.25 * (m_decay * (m_in_4 + Sqr(m_f_sq - m_v_sq)
                                         - 2. * m_in_sq * (m_f_sq + m_v_sq)))
      / m_v_sq;

   const double c6 = 0.25 * (m_decay * m_fermion * (
                                m_in_4 + Sqr(m_f_sq - m_v_sq) - 2. * m_in_sq * (
                                   m_f_sq + m_v_sq))) / m_v_sq;

   return c1 * (AbsSqr(matrix_element_gam_left) + AbsSqr(matrix_element_gam_right))
      + c2 * (AbsSqr(matrix_element_p_1) + AbsSqr(matrix_element_p_2))
      + 2. * c3 * Re(matrix_element_gam_left * Conj(matrix_element_gam_right))
      + 2. * c4 * Re(matrix_element_gam_left * Conj(matrix_element_p_1)
                     + matrix_element_gam_right * Conj(matrix_element_p_2))
      + 2. * c5 * Re(matrix_element_gam_left * Conj(matrix_element_p_2)
                     + matrix_element_gam_right * Conj(matrix_element_p_1))
      + 2. * c6 * Re(matrix_element_p_1 * Conj(matrix_element_p_2));
}

} // namespace flexiblesusy
