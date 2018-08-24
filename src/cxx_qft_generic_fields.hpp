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

#ifndef CXX_QFT_GENERIC_FIELDS_H
#define CXX_QFT_GENERIC_FIELDS_H

#include "multiindex.hpp"

#include <array>

#include <boost/array.hpp>
#include <boost/mpl/at.hpp>
#include <boost/version.hpp>
#include <boost/mpl/pair.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/vector_c.hpp>

#include <boost/range/iterator_range.hpp>

#include <boost/fusion/include/copy.hpp>
#include <boost/fusion/include/move.hpp>
#include <boost/fusion/adapted/boost_array.hpp>

namespace flexiblesusy {

namespace cxx_qft {

template<class Field>
struct field_indices {
   using type = std::array<int, Field::numberOfFieldIndices>;
};

namespace fields {

enum class ParticleType {
   scalar,
   fermion,
   vector,
   ghost
};

template<typename Field>
struct is_fermion {
   static constexpr bool value = Field::particle_type == ParticleType::fermion;
};

template<typename Field>
struct is_vector {
   static constexpr bool value = Field::particle_type == ParticleType::vector;
};

template<typename Field>
struct is_scalar {
   static constexpr bool value = Field::particle_type == ParticleType::scalar;
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

} // namespace fields

namespace detail {

template<class Begin, class End>
decltype(
   boost::make_iterator_range(
      multiindex<Begin, End>::begin(),
      multiindex<Begin, End>::end()
      )
   )
make_index_range( void )
{
   using index = multiindex<Begin, End>;

   return boost::make_iterator_range(
      index::begin(), index::end()
      );
}

} // namespace detail

template<class Field>
typename std::enable_if<
  Field::numberOfGenerations != 1,
  bool
>::type
isSMField(const typename field_indices<Field>::type& indices)
{
   boost::array<bool, Field::numberOfGenerations> sm_flags;

#if BOOST_VERSION >= 105800
   boost::fusion::move(typename Field::sm_flags(), sm_flags);
#else
   boost::fusion::copy(typename Field::sm_flags(), sm_flags);
#endif

   return sm_flags[indices[0]];
}

template<class Field>
typename std::enable_if<
  Field::numberOfGenerations == 1,
  bool
>::type
isSMField(const typename field_indices<Field>::type &)
{
   return boost::mpl::at_c<
      typename Field::sm_flags,
      0
      >::type::value;
}

template<class ObjectWithIndexBounds>
struct index_bounds {
   using type = typename ObjectWithIndexBounds::index_bounds;
};

template<class ObjectWithIndexBounds>
decltype(detail::make_index_range<
  typename ObjectWithIndexBounds::index_bounds::first,
  typename ObjectWithIndexBounds::index_bounds::second
>())
index_range( void )
{
  return detail::make_index_range<
    typename ObjectWithIndexBounds::index_bounds::first,
    typename ObjectWithIndexBounds::index_bounds::second
  >();
}

} // namespace cxx_qft

} // namespace flexiblesusy

#endif
