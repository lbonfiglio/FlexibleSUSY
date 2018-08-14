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

#include "decay.hpp"
#include "error.hpp"

#include <sstream>

namespace flexiblesusy {

Decay::Decay(
   int initial_pdg_, std::initializer_list<int> product_pdgs_, double width_)
   : initial_pdg(initial_pdg_)
   , product_pdgs(product_pdgs_)
   , width(width_)
{
}

Decays_list::Decays_list(int initial_pdg_)
   : initial_pdg(initial_pdg_)
{
}

void Decays_list::clear()
{
   decays.clear();
   total_width = 0.;
}

void Decays_list::set_decay(
   double width, std::initializer_list<int> product_pdgs)
{
   const auto decay_id = decay_hash(product_pdgs);

   const auto pos = decays.find(decay_id);
   // @todo more efficient implementation
   if (pos != std::end(decays)) {
      // @todo not currently checking initial states match
      // @todo if worried about failure of key creation,
      // could check that PDGs match
      total_width -= pos->second.get_width();
      pos->second.set_width(width);
   } else {
      decays.insert(std::make_pair(decay_id,
                                   Decay(initial_pdg, product_pdgs, width)));
   }
   total_width += width;
}

const Decay& Decays_list::get_decay(
   std::initializer_list<int> product_pdgs) const
{
   const auto decay_id = decay_hash(product_pdgs);

   const auto pos = decays.find(decay_id);

   if (pos == std::end(decays)) {
      std::ostringstream msg;
      msg << "Decay of particle " << initial_pdg
          << " into particles {";

      std::size_t count = 0;
      const std::size_t n_final = product_pdgs.size();
      for (auto it = std::begin(product_pdgs),
              end = std::end(product_pdgs); it != end; ++it, ++count) {
         msg << *it;
         if (count != n_final) {
            msg << ", ";
         }
      }

      msg << "} not found\n";

      throw OutOfBoundsError(msg.str());
   }

   return pos->second;
}

} // namespace flexiblesusy
