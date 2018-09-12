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
   return AbsSqr(form_factor);
}

double Decay_amplitude_SSV::square() const
{
   if (m_vector <= massless_vector_threshold) {
      return -2. * AbsSqr(form_factor) * (Sqr(m_decay) + Sqr(m_scalar));
   }

   const double m_in_sq = Sqr(m_decay);
   const double m_s_sq = Sqr(m_scalar);
   const double m_v_sq = Sqr(m_vector);

   return (Sqr(m_in_sq) + Sqr(m_v_sq - m_s_sq)
           - 2. * m_in_sq * (m_s_sq + m_v_sq))
      * AbsSqr(form_factor) / m_v_sq;
}

// @todo handle massless vectors safely
// @todo check expressions correct
double Decay_amplitude_SVV::square() const
{
   if (m_out_1 <= massless_vector_threshold &&
       m_out_2 <= massless_vector_threshold) {
      const double m0 = 4. * AbsSqr(form_factor_g);

      const double m2 = Sqr(m_decay) * (
         Re(form_factor_g * Conj(form_factor_12))
         + Re(form_factor_g * Conj(form_factor_21)));

      const double m4 = 0.5 * Power4(m_decay) * (
            Re(form_factor_11 * Conj(form_factor_22))
            + Re(form_factor_12 * Conj(form_factor_21))
            + AbsSqr(form_factor_eps));

      return m0 + m2 + m4;
   } else if (m_out_1 <= massless_vector_threshold) {
      const double m0 = 3. * AbsSqr(form_factor_g);

      const double m_s_sq = Sqr(m_decay);
      const double m_vec_sq = Sqr(m_out_2);

      const double m41 = 0.25 * Sqr(m_s_sq - m_vec_sq) * (
         2. * AbsSqr(form_factor_eps) - AbsSqr(form_factor_21));

      const double m42 = -0.5 * Sqr(m_s_sq - m_vec_sq) *
         Re(form_factor_g * Conj(form_factor_11)) / m_vec_sq;

      const double m43 = -0.25 * Power3(m_s_sq - m_vec_sq) *
         Re(form_factor_11 * Conj(form_factor_21)) / m_vec_sq;

      return m0 + m41 + m42 + m43;
   } else if (m_out_2 <= massless_vector_threshold) {
      const double m0 = 3. * AbsSqr(form_factor_g);

      const double m_s_sq = Sqr(m_decay);
      const double m_vec_sq = Sqr(m_out_1);

      const double m41 = 0.25 * Sqr(m_s_sq - m_vec_sq) * (
         2. * AbsSqr(form_factor_eps) - AbsSqr(form_factor_21));

      const double m42 = -0.5 * Sqr(m_s_sq - m_vec_sq) *
         Re(form_factor_g * Conj(form_factor_22)) / m_vec_sq;

      const double m43 = -0.25 * Power3(m_s_sq - m_vec_sq) *
         Re(form_factor_22 * Conj(form_factor_21)) / m_vec_sq;

      return m0 + m41 + m42 + m43;
   }

   const double m_s_sq = Sqr(m_decay);
   const double m_s_4 = Power4(m_decay);
   const double m_1_sq = Sqr(m_out_1);
   const double m_1_4 = Power4(m_out_1);
   const double m_2_sq = Sqr(m_out_2);
   const double m_2_4 = Power4(m_out_2);

   const double mgg = 0.25 * AbsSqr(form_factor_g) * (
      m_s_4 + m_1_4 + m_2_4 + 10. * m_1_sq * m_2_sq
      - 2. * m_s_sq * (m_1_sq + m_2_sq)) / (m_1_sq * m_2_sq);

   const double m33 = 0.0625 * Sqr(m_s_4 + Sqr(m_1_sq - m_2_sq)
                                   - 2. * m_s_sq * (m_1_sq + m_2_sq))
      * AbsSqr(form_factor_21) / (m_1_sq * m_2_sq);

   const double mepseps = AbsSqr(form_factor_eps) * (
      0.5 * Sqr(m_s_sq - m_1_sq - m_2_sq) - 2. * m_1_sq * m_2_sq);

   const double mg3 = 0.25 * (m_s_4 * m_s_sq - 3. * m_s_4 * (m_1_sq + m_2_sq)
                              - Sqr(m_1_sq - m_2_sq) * (m_1_sq + m_2_sq)
                              + m_s_sq * (3. * m_1_4 + 2. * m_1_sq * m_2_sq
                                          + 3. * m_2_4))
      * Re(form_factor_g * Conj(form_factor_21)) / (m_1_sq * m_2_sq);

   return mgg + m33 + mepseps + mg3;
}

double Decay_amplitude_SFF::square() const
{
   const double m_in_sq = Sqr(m_decay);
   const double m_1_sq = Sqr(m_out_1);
   const double m_2_sq = Sqr(m_out_2);

   return (m_in_sq - m_1_sq - m_2_sq) *
      (AbsSqr(form_factor_left) + AbsSqr(form_factor_right))
      - 4. * m_out_1 * m_out_2 * Re(
         form_factor_left * Conj(form_factor_right));
}

double Decay_amplitude_FFS::square() const
{
   const double m_in_sq = Sqr(m_decay);
   const double m_f_sq = Sqr(m_fermion);
   const double m_s_sq = Sqr(m_scalar);

   return 0.5 * (m_in_sq + m_f_sq - m_s_sq) *
      (AbsSqr(form_factor_left) + AbsSqr(form_factor_right))
      + 2. * m_fermion * m_scalar * Re(
         form_factor_left * Conj(form_factor_right));
}

// @todo handle massless vectors safely
double Decay_amplitude_FFV::square() const
{
   if (m_vector <= massless_vector_threshold) {
      return 0.;
   }

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

   return c1 * (AbsSqr(form_factor_gam_left) + AbsSqr(form_factor_gam_right))
      + c2 * (AbsSqr(form_factor_p_1) + AbsSqr(form_factor_p_2))
      + 2. * c3 * Re(form_factor_gam_left * Conj(form_factor_gam_right))
      + 2. * c4 * Re(form_factor_gam_left * Conj(form_factor_p_1)
                     + form_factor_gam_right * Conj(form_factor_p_2))
      + 2. * c5 * Re(form_factor_gam_left * Conj(form_factor_p_2)
                     + form_factor_gam_right * Conj(form_factor_p_1))
      + 2. * c6 * Re(form_factor_p_1 * Conj(form_factor_p_2));
}

} // namespace flexiblesusy
