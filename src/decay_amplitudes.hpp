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
   std::complex<double> matrix_element_gam_1{};
   std::complex<double> matrix_element_gam_2{};
   std::complex<double> matrix_element_p_1{};
   std::complex<double> matrix_element_p_2{};

   double square() const;
};

namespace detail {

template <class Field_in, class Field_out_1, class Field_out_2, class Amplitude_type = void>
struct Two_body_decay_amplitude_type { };

template <class Field_in, class Field_out_1, class Field_out_2>
struct Two_body_decay_amplitude_type<Field_in, Field_out_1, Field_out_2,
                                     typename std::enable_if<cxx_qft::is_scalar<Field_in>::value &&
                                                             cxx_qft::is_scalar<Field_out_1>::value &&
                                                             cxx_qft::is_scalar<Field_out_2>::value>::type > {
   using type = Decay_amplitude_SSS;
};

template <class Field_in, class Field_out_1, class Field_out_2>
struct Two_body_decay_amplitude_type<Field_in, Field_out_1, Field_out_2,
                                     typename std::enable_if<cxx_qft::is_scalar<Field_in>::value &&
                                                             cxx_qft::is_scalar<Field_out_1>::value &&
                                                             cxx_qft::is_vector<Field_out_2>::value>::type> {
   using type = Decay_amplitude_SSV;
};

template <class Field_in, class Field_out_1, class Field_out_2>
struct Two_body_decay_amplitude_type<Field_in, Field_out_1, Field_out_2,
                                     typename std::enable_if<cxx_qft::is_scalar<Field_in>::value &&
                                                             cxx_qft::is_vector<Field_out_1>::value &&
                                                             cxx_qft::is_scalar<Field_out_2>::value>::type> {
   using type = Decay_amplitude_SSV;
};

template <class Field_in, class Field_out_1, class Field_out_2>
struct Two_body_decay_amplitude_type<Field_in, Field_out_1, Field_out_2,
                                     typename std::enable_if<cxx_qft::is_scalar<Field_in>::value &&
                                                             cxx_qft::is_vector<Field_out_1>::value &&
                                                             cxx_qft::is_vector<Field_out_2>::value>::type> {
   using type = Decay_amplitude_SVV;
};

template <class Field_in, class Field_out_1, class Field_out_2>
struct Two_body_decay_amplitude_type<Field_in, Field_out_1, Field_out_2,
                                     typename std::enable_if<cxx_qft::is_scalar<Field_in>::value &&
                                                             cxx_qft::is_fermion<Field_out_1>::value &&
                                                             cxx_qft::is_fermion<Field_out_2>::value>::type> {
   using type = Decay_amplitude_SFF;
};

template <class Field_in, class Field_out_1, class Field_out_2>
struct Two_body_decay_amplitude_type<Field_in, Field_out_1, Field_out_2,
                                     typename std::enable_if<cxx_qft::is_fermion<Field_in>::value &&
                                                             cxx_qft::is_fermion<Field_out_1>::value &&
                                                             cxx_qft::is_scalar<Field_out_2>::value>::type> {
   using type = Decay_amplitude_FFS;
};

template <class Field_in, class Field_out_1, class Field_out_2>
struct Two_body_decay_amplitude_type<Field_in, Field_out_1, Field_out_2,
                                     typename std::enable_if<cxx_qft::is_fermion<Field_in>::value &&
                                                             cxx_qft::is_scalar<Field_out_1>::value &&
                                                             cxx_qft::is_fermion<Field_out_2>::value>::type> {
   using type = Decay_amplitude_FFS;
};

template <class Field_in, class Field_out_1, class Field_out_2>
struct Two_body_decay_amplitude_type<Field_in, Field_out_1, Field_out_2,
                                     typename std::enable_if<cxx_qft::is_fermion<Field_in>::value &&
                                                             cxx_qft::is_fermion<Field_out_1>::value &&
                                                             cxx_qft::is_vector<Field_out_2>::value>::type> {
   using type = Decay_amplitude_FFV;
};

template <class Field_in, class Field_out_1, class Field_out_2>
struct Two_body_decay_amplitude_type<Field_in, Field_out_1, Field_out_2,
                                typename std::enable_if<cxx_qft::is_fermion<Field_in>::value &&
                                                        cxx_qft::is_vector<Field_out_1>::value &&
                                                        cxx_qft::is_fermion<Field_out_2>::value>::type> {
   using type = Decay_amplitude_FFV;
};

template <class... Fields>
struct Decay_amplitude_type {
   using type = typename std::enable_if<sizeof...(Fields) == 3,
                                        typename Two_body_decay_amplitude_type<Fields...>::type >::type;
};

} // namespace detail

template <class Field_in, class... Fields_out>
struct Decay_amplitude {
   using amplitude_type = typename detail::Decay_amplitude_type<Field_in, Fields_out...>::type;
   amplitude_type amplitude{};

   double square() const { return amplitude.square(); }
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
                        Decay_amplitude<Field_in, Field_out_1, Field_out_2> >::type
tree_level_decay_amplitude(const Vertex& vertex)
{
   Decay_amplitude_SSS amplitude;
   amplitude.matrix_element = vertex.value();

   Decay_amplitude<Field_in, Field_out_1, Field_out_2> result;
   result.amplitude = amplitude;

   return result;
}

template <class Field_in, class Field_out_1, class Field_out_2, class Vertex>
typename std::enable_if<cxx_qft::is_scalar<Field_in>::value &&
                        cxx_qft::is_scalar<Field_out_1>::value &&
                        cxx_qft::is_vector<Field_out_2>::value,
                        Decay_amplitude<Field_in, Field_out_1, Field_out_2> >::type
tree_level_decay_amplitude(const Vertex& vertex)
{
   Decay_amplitude_SSV amplitude;
   amplitude.matrix_element = vertex.value(1, 2);

   Decay_amplitude<Field_in, Field_out_1, Field_out_2> result;
   result.amplitude = amplitude;

   return result;
}

template <class Field_in, class Field_out_1, class Field_out_2, class Vertex>
typename std::enable_if<cxx_qft::is_scalar<Field_in>::value &&
                        cxx_qft::is_vector<Field_out_1>::value &&
                        cxx_qft::is_scalar<Field_out_2>::value,
                        Decay_amplitude<Field_in, Field_out_1, Field_out_2> >::type
tree_level_decay_amplitude(const Vertex& vertex)
{
   const auto reordered =
      tree_level_decay_amplitude<Field_in, Field_out_2, Field_out_1>(vertex);

   Decay_amplitude<Field_in, Field_out_1, Field_out_2> result;
   result.amplitude = reordered.amplitude;

   return result;
}

template <class Field_in, class Field_out_1, class Field_out_2, class Vertex>
typename std::enable_if<cxx_qft::is_scalar<Field_in>::value &&
                        cxx_qft::is_vector<Field_out_1>::value &&
                        cxx_qft::is_vector<Field_out_2>::value,
                        Decay_amplitude<Field_in, Field_out_1, Field_out_2> >::type
tree_level_decay_amplitude(const Vertex& vertex)
{
   Decay_amplitude_SVV amplitude;
   amplitude.M1 = vertex.value();

   Decay_amplitude<Field_in, Field_out_1, Field_out_2> result;
   result.amplitude = amplitude;

   return result;
}

template <class Field_in, class Field_out_1, class Field_out_2, class Vertex>
typename std::enable_if<cxx_qft::is_scalar<Field_in>::value &&
                        cxx_qft::is_fermion<Field_out_1>::value &&
                        cxx_qft::is_fermion<Field_out_2>::value,
                        Decay_amplitude<Field_in, Field_out_1, Field_out_2> >::type
tree_level_decay_amplitude(const Vertex& vertex)
{
   Decay_amplitude_SFF amplitude;

   amplitude.matrix_element_left = vertex.left();
   amplitude.matrix_element_right = vertex.right();

   Decay_amplitude<Field_in, Field_out_1, Field_out_2> result;
   result.amplitude = amplitude;

   return result;
}

template <class Field_in, class Field_out_1, class Field_out_2, class Vertex>
typename std::enable_if<cxx_qft::is_fermion<Field_in>::value &&
                        cxx_qft::is_fermion<Field_out_1>::value &&
                        cxx_qft::is_scalar<Field_out_2>::value,
                        Decay_amplitude<Field_in, Field_out_1, Field_out_2> >::type
tree_level_decay_amplitude(const Vertex& vertex)
{
   Decay_amplitude<Field_in, Field_out_1, Field_out_2> result;
   return result;
}

template <class Field_in, class Field_out_1, class Field_out_2, class Vertex>
typename std::enable_if<cxx_qft::is_fermion<Field_in>::value &&
                        cxx_qft::is_scalar<Field_out_1>::value &&
                        cxx_qft::is_fermion<Field_out_2>::value,
                        Decay_amplitude<Field_in, Field_out_1, Field_out_2> >::type
tree_level_decay_amplitude(const Vertex& vertex)
{
   const auto reordered =
      tree_level_decay_amplitude<Field_in, Field_out_2, Field_out_1>(vertex);

   Decay_amplitude<Field_in, Field_out_1, Field_out_2> result;
   result.amplitude = reordered.amplitude;

   return result;
}

template <class Field_in, class Field_out_1, class Field_out_2, class Vertex>
typename std::enable_if<cxx_qft::is_fermion<Field_in>::value &&
                        cxx_qft::is_fermion<Field_out_1>::value &&
                        cxx_qft::is_vector<Field_out_2>::value,
                        Decay_amplitude<Field_in, Field_out_1, Field_out_2> >::type
tree_level_decay_amplitude(const Vertex& vertex)
{
   Decay_amplitude<Field_in, Field_out_1, Field_out_2> result;
   return result;
}

template <class Field_in, class Field_out_1, class Field_out_2, class Vertex>
typename std::enable_if<cxx_qft::is_fermion<Field_in>::value &&
                        cxx_qft::is_vector<Field_out_1>::value &&
                        cxx_qft::is_fermion<Field_out_2>::value,
                        Decay_amplitude<Field_in, Field_out_1, Field_out_2> >::type
tree_level_decay_amplitude(const Vertex& vertex)
{
   const auto reordered =
      tree_level_decay_amplitude<Field_in, Field_out_2, Field_out_1>(vertex);

   Decay_amplitude<Field_in, Field_out_1, Field_out_2> result;
   result.amplitude = reordered.amplitude;

   return result;
}

} // namespace flexiblesusy

#endif
