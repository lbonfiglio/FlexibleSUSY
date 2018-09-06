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

#include "field_traits.hpp"

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
   std::complex<double> matrix_element_gam_left{};
   std::complex<double> matrix_element_gam_right{};
   std::complex<double> matrix_element_p_1{};
   std::complex<double> matrix_element_p_2{};

   double square() const;
};

namespace detail {

template <class Field_in, class Field_out_1, class Field_out_2,
          class Amplitude_type = void>
struct Two_body_decay_amplitude_type { };

template <class Field_in, class Field_out_1, class Field_out_2>
struct Two_body_decay_amplitude_type<
   Field_in, Field_out_1, Field_out_2,
   typename std::enable_if<cxx_qft::is_scalar<Field_in>::value &&
                           cxx_qft::is_scalar<Field_out_1>::value &&
                           cxx_qft::is_scalar<Field_out_2>::value>::type > {
   using type = Decay_amplitude_SSS;
};

template <class Field_in, class Field_out_1, class Field_out_2>
struct Two_body_decay_amplitude_type<
   Field_in, Field_out_1, Field_out_2,
   typename std::enable_if<cxx_qft::is_scalar<Field_in>::value &&
                           cxx_qft::is_scalar<Field_out_1>::value &&
                           cxx_qft::is_vector<Field_out_2>::value>::type> {
   using type = Decay_amplitude_SSV;
};

template <class Field_in, class Field_out_1, class Field_out_2>
struct Two_body_decay_amplitude_type<
   Field_in, Field_out_1, Field_out_2,
   typename std::enable_if<cxx_qft::is_scalar<Field_in>::value &&
                           cxx_qft::is_vector<Field_out_1>::value &&
                           cxx_qft::is_scalar<Field_out_2>::value>::type> {
   using type = Decay_amplitude_SSV;
};

template <class Field_in, class Field_out_1, class Field_out_2>
struct Two_body_decay_amplitude_type<
   Field_in, Field_out_1, Field_out_2,
   typename std::enable_if<cxx_qft::is_scalar<Field_in>::value &&
                           cxx_qft::is_vector<Field_out_1>::value &&
                           cxx_qft::is_vector<Field_out_2>::value>::type> {
   using type = Decay_amplitude_SVV;
};

template <class Field_in, class Field_out_1, class Field_out_2>
struct Two_body_decay_amplitude_type<
   Field_in, Field_out_1, Field_out_2,
   typename std::enable_if<cxx_qft::is_scalar<Field_in>::value &&
                           cxx_qft::is_fermion<Field_out_1>::value &&
                           cxx_qft::is_fermion<Field_out_2>::value>::type> {
   using type = Decay_amplitude_SFF;
};

template <class Field_in, class Field_out_1, class Field_out_2>
struct Two_body_decay_amplitude_type<
   Field_in, Field_out_1, Field_out_2,
   typename std::enable_if<cxx_qft::is_fermion<Field_in>::value &&
                           cxx_qft::is_fermion<Field_out_1>::value &&
                           cxx_qft::is_scalar<Field_out_2>::value>::type> {
   using type = Decay_amplitude_FFS;
};

template <class Field_in, class Field_out_1, class Field_out_2>
struct Two_body_decay_amplitude_type<
   Field_in, Field_out_1, Field_out_2,
   typename std::enable_if<cxx_qft::is_fermion<Field_in>::value &&
                           cxx_qft::is_scalar<Field_out_1>::value &&
                           cxx_qft::is_fermion<Field_out_2>::value>::type> {
   using type = Decay_amplitude_FFS;
};

template <class Field_in, class Field_out_1, class Field_out_2>
struct Two_body_decay_amplitude_type<
   Field_in, Field_out_1, Field_out_2,
   typename std::enable_if<cxx_qft::is_fermion<Field_in>::value &&
                           cxx_qft::is_fermion<Field_out_1>::value &&
                           cxx_qft::is_vector<Field_out_2>::value>::type> {
   using type = Decay_amplitude_FFV;
};

template <class Field_in, class Field_out_1, class Field_out_2>
struct Two_body_decay_amplitude_type<
   Field_in, Field_out_1, Field_out_2,
   typename std::enable_if<cxx_qft::is_fermion<Field_in>::value &&
                           cxx_qft::is_vector<Field_out_1>::value &&
                           cxx_qft::is_fermion<Field_out_2>::value>::type> {
   using type = Decay_amplitude_FFV;
};

} // namespace detail

/**
 * @class Decay_amplitude_type
 * @brief helper class to determine amplitude type for a given set of fields
 */
template <class... Fields>
struct Decay_amplitude_type {
   using type =
      typename std::enable_if<
      sizeof...(Fields) == 3,
      typename detail::Two_body_decay_amplitude_type<Fields...>::type >::type;
};

template <typename Amplitude>
double square_amplitude(const Amplitude& a)
{
   return a.square();
}

// @todo implement all of the following with the correct conventions
template <class Field_in, class Field_out_1, class Field_out_2, class Vertex>
typename std::enable_if<cxx_qft::is_scalar<Field_in>::value &&
                        cxx_qft::is_scalar<Field_out_1>::value &&
                        cxx_qft::is_scalar<Field_out_2>::value,
                        Decay_amplitude_SSS>::type
tree_level_decay_amplitude(double m_decay, double m_out_1, double m_out_2,
                           const Vertex& vertex)
{
   Decay_amplitude_SSS amplitude;

   amplitude.m_decay = m_decay;
   amplitude.m_out_1 = m_out_1;
   amplitude.m_out_2 = m_out_2;

   amplitude.matrix_element = vertex.value();

   return amplitude;
}

template <class Field_in, class Field_out_1, class Field_out_2, class Vertex>
typename std::enable_if<cxx_qft::is_scalar<Field_in>::value &&
                        cxx_qft::is_scalar<Field_out_1>::value &&
                        cxx_qft::is_vector<Field_out_2>::value,
                        Decay_amplitude_SSV>::type
tree_level_decay_amplitude(double m_decay, double m_scalar, double m_vector,
                           const Vertex& vertex)
{
   Decay_amplitude_SSV amplitude;

   amplitude.m_decay = m_decay;
   amplitude.m_scalar = m_scalar;
   amplitude.m_vector = m_vector;

   amplitude.matrix_element = vertex.value(1, 2);

   return amplitude;
}

template <class Field_in, class Field_out_1, class Field_out_2, class Vertex>
typename std::enable_if<cxx_qft::is_scalar<Field_in>::value &&
                        cxx_qft::is_vector<Field_out_1>::value &&
                        cxx_qft::is_scalar<Field_out_2>::value,
                        Decay_amplitude_SSV>::type
tree_level_decay_amplitude(double m_decay, double m_vector, double m_scalar,
                           const Vertex& vertex)
{
   return tree_level_decay_amplitude<Field_in, Field_out_2, Field_out_1>(
      m_decay, m_scalar, m_vector, vertex);
}

template <class Field_in, class Field_out_1, class Field_out_2, class Vertex>
typename std::enable_if<cxx_qft::is_scalar<Field_in>::value &&
                        cxx_qft::is_vector<Field_out_1>::value &&
                        cxx_qft::is_vector<Field_out_2>::value,
                        Decay_amplitude_SVV>::type
tree_level_decay_amplitude(double m_decay, double m_out_1, double m_out_2,
                           const Vertex& vertex)
{
   Decay_amplitude_SVV amplitude;

   amplitude.m_decay = m_decay;
   amplitude.m_out_1 = m_out_1;
   amplitude.m_out_2 = m_out_2;

   amplitude.M1 = vertex.value();

   return amplitude;
}

template <class Field_in, class Field_out_1, class Field_out_2, class Vertex>
typename std::enable_if<cxx_qft::is_scalar<Field_in>::value &&
                        cxx_qft::is_fermion<Field_out_1>::value &&
                        cxx_qft::is_fermion<Field_out_2>::value,
                        Decay_amplitude_SFF>::type
tree_level_decay_amplitude(double m_decay, double m_out_1, double m_out_2,
                           const Vertex& vertex)
{
   Decay_amplitude_SFF amplitude;

   amplitude.m_decay = m_decay;
   amplitude.m_out_1 = m_out_1;
   amplitude.m_out_2 = m_out_2;

   // @todo check sign convention in SARAH
   amplitude.matrix_element_left = vertex.left();
   amplitude.matrix_element_right = vertex.right();

   return amplitude;
}

template <class Field_in, class Field_out_1, class Field_out_2, class Vertex>
typename std::enable_if<cxx_qft::is_fermion<Field_in>::value &&
                        cxx_qft::is_fermion<Field_out_1>::value &&
                        cxx_qft::is_scalar<Field_out_2>::value,
                        Decay_amplitude_FFS>::type
tree_level_decay_amplitude(double m_decay, double m_fermion, double m_scalar,
                           const Vertex& vertex)
{
   Decay_amplitude_FFS amplitude;

   amplitude.m_decay = m_decay;
   amplitude.m_fermion = m_fermion;
   amplitude.m_scalar = m_scalar;

   amplitude.matrix_element_left = vertex.left();
   amplitude.matrix_element_right = vertex.right();

   return amplitude;
}

template <class Field_in, class Field_out_1, class Field_out_2, class Vertex>
typename std::enable_if<cxx_qft::is_fermion<Field_in>::value &&
                        cxx_qft::is_scalar<Field_out_1>::value &&
                        cxx_qft::is_fermion<Field_out_2>::value,
                        Decay_amplitude_FFS>::type
tree_level_decay_amplitude(double m_decay, double m_scalar, double m_fermion,
                           const Vertex& vertex)
{
   return tree_level_decay_amplitude<Field_in, Field_out_2, Field_out_1>(
      m_decay, m_fermion, m_scalar, vertex);
}

template <class Field_in, class Field_out_1, class Field_out_2, class Vertex>
typename std::enable_if<cxx_qft::is_fermion<Field_in>::value &&
                        cxx_qft::is_fermion<Field_out_1>::value &&
                        cxx_qft::is_vector<Field_out_2>::value,
                        Decay_amplitude_FFV>::type
tree_level_decay_amplitude(double m_decay, double m_fermion, double m_vector,
                           const Vertex& vertex)
{
   Decay_amplitude_FFV amplitude;

   amplitude.m_decay = m_decay;
   amplitude.m_fermion = m_fermion;
   amplitude.m_vector = m_vector;

   amplitude.matrix_element_gam_left = vertex.left();
   amplitude.matrix_element_gam_right = vertex.right();

   return amplitude;
}

template <class Field_in, class Field_out_1, class Field_out_2, class Vertex>
typename std::enable_if<cxx_qft::is_fermion<Field_in>::value &&
                        cxx_qft::is_vector<Field_out_1>::value &&
                        cxx_qft::is_fermion<Field_out_2>::value,
                        Decay_amplitude_FFV>::type
tree_level_decay_amplitude(double m_decay, double m_vector, double m_fermion,
                           const Vertex& vertex)
{
   return tree_level_decay_amplitude<Field_in, Field_out_2, Field_out_1>(
      m_decay, m_fermion, m_vector, vertex);
}

} // namespace flexiblesusy

#endif
