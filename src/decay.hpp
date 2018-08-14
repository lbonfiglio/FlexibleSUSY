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

#ifndef DECAY_H
#define DECAY_H

#include <algorithm>
#include <cstddef>
#include <initializer_list>
#include <map>
#include <vector>

namespace flexiblesusy {

namespace detail {

template <class Container>
std::size_t decay_hash_impl(Container pdgs)
{
   Container sorted(pdgs);
   std::sort(sorted.begin(), sorted.end());

   std::size_t hash_value = 0;
   return hash_value;
}

} // namespace detail

template <class T>
std::size_t decay_hash(T&& pdgs)
{
   return detail::decay_hash_impl(std::forward<T>(pdgs));
}

template <class T>
std::size_t decay_hash(std::initializer_list<T> pdgs)
{
   return detail::decay_hash_impl(std::vector<T>(pdgs));
}

class Decay {
public:
   Decay(int, std::initializer_list<int>, double);

   int get_initial_particle() const { return initial_pdg; }
   const std::vector<int>& get_final_state_particles() const {
      return product_pdgs;
   }
   std::size_t get_final_state_size() const { return product_pdgs.size(); }
   double get_width() const { return width; }

   void set_width(double w) { width = w; }

private:
   int initial_pdg{0};
   std::vector<int> product_pdgs{};
   double width{0.};
};

class Decays_list {
public:
   explicit Decays_list(int particle_id);
   ~Decays_list() = default;

   void clear();
   void set_decay(double width, std::initializer_list<int> products);
   const Decay& get_decay(std::initializer_list<int> products) const;
   double get_total_width() const { return total_width; }

private:
   int initial_pdg{0};
   std::map<std::size_t, Decay> decays{};
   double total_width{0.};
};

} // namespace flexiblesusy

#endif
