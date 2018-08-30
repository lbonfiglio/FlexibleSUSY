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

namespace cxx_qft {

namespace fields {

template <class T>
struct conj;

template <class T>
struct bar;

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

} // namespace fields

} // namespace cxx_qft

} // namespace flexiblesusy

#endif
