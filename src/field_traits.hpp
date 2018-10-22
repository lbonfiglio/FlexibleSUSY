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

#ifndef FIELD_TRAITS_H
#define FIELD_TRAITS_H

#include <type_traits>

namespace flexiblesusy {

namespace field_traits {

enum class ParticleType {
   scalar,
   fermion,
   vector,
   ghost
};

template<typename Field>
struct is_massless {
   static constexpr bool value = Field::massless;
};

enum class ParticleColorRep {
   singlet,
   triplet,
   anti_triplet,
   sextet,
   octet
};
template<typename Field> struct is_triplet {
   static constexpr bool value = Field::color_rep == ParticleColorRep::triplet;
};
template<typename Field> struct is_anti_triplet {
   static constexpr bool value =
      Field::color_rep == ParticleColorRep::anti_triplet;
};
template<typename Field>
constexpr typename std::enable_if<
   is_triplet<Field>::value, ParticleColorRep
   >::type
color_conj() {
   return ParticleColorRep::anti_triplet;
}
template<typename Field>
constexpr typename std::enable_if<
   is_anti_triplet<Field>::value, ParticleColorRep
   >::type
color_conj() {
   return ParticleColorRep::triplet;
}
template<typename Field>
constexpr typename std::enable_if<
   !is_triplet<Field>::value && !is_anti_triplet<Field>::value, ParticleColorRep
   >::type
color_conj() {
   return Field::color_rep;
}

template<class Field>
struct bar {
   using index_bounds = typename Field::index_bounds;
   using sm_flags = typename Field::sm_flags;
   using lorentz_conjugate = Field;
   using type = bar<Field>;

   static constexpr int numberOfGenerations = Field::numberOfGenerations;
   static constexpr int numberOfFieldIndices = Field::numberOfFieldIndices;
   static constexpr double electric_charge = Field::electric_charge;
   static constexpr auto particle_type = Field::particle_type;
   static constexpr auto color_rep = color_conj<Field>();
   static constexpr auto massless = Field::massless;
};

template<class Field>
struct conj {
   using index_bounds = typename Field::index_bounds;
   using sm_flags = typename Field::sm_flags;
   using lorentz_conjugate = Field;
   using type = conj<Field>;

   static constexpr int numberOfGenerations = Field::numberOfGenerations;
   static constexpr int numberOfFieldIndices = Field::numberOfFieldIndices;
   static constexpr double electric_charge = Field::electric_charge;
   static constexpr auto particle_type = Field::particle_type;
   static constexpr auto color_rep = color_conj<Field>();
   static constexpr auto massless = Field::massless;
};

template<class Field>
struct bar<bar<Field> > {
   using type = Field;
};

template<class Field>
struct conj<conj<Field> > {
   using type = Field;
};

template<class Field>
struct remove_lorentz_conjugation {
   using type = Field;
};

template<class Field>
struct remove_lorentz_conjugation<bar<Field>> {
   using type = Field;
};

template<class Field>
struct remove_lorentz_conjugation<conj<Field>> {
   using type = Field;
};

template <class Field>
struct is_scalar : public std::false_type {};

template <class Field>
struct is_scalar<bar<Field> > : public is_scalar<Field> {};

template <class Field>
struct is_scalar<conj<Field> > : public is_scalar<Field> {};

template <class Field>
struct is_fermion : public std::false_type {};

template <class Field>
struct is_fermion<bar<Field> > : public is_fermion<Field> {};

template <class Field>
struct is_fermion<conj<Field> > : public is_fermion<Field> {};

template <class Field>
struct is_vector : public std::false_type {};

template <class Field>
struct is_vector<bar<Field> > : public is_vector<Field> {};

template <class Field>
struct is_vector<conj<Field> > : public is_vector<Field> {};

template <class Field>
struct is_ghost : public std::false_type {};

template <class Field>
struct is_ghost<bar<Field> > : public is_ghost<Field> {};

template <class Field>
struct is_ghost<conj<Field> > : public is_ghost<Field> {};

} // namespace field_traits

} // namespace flexiblesusy

#endif
