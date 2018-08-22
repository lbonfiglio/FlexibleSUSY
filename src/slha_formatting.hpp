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

#ifndef SLHA_FORMATTING_H
#define SLHA_FORMATTING_H

#include <string>
#include <boost/format.hpp>

namespace flexiblesusy {

inline boost::format format_mass(int pdg, double value, const std::string& name)
{
   return boost::format(" %9d   %16.8E   # %s\n") % pdg % value % name;
}

inline boost::format format_tensor(double value, const char* name)
{
   return boost::format("   %16.8E   # %s\n") % value % name;
}

inline boost::format format_tensor(double value, const std::string& name)
{
   return boost::format("   %16.8E   # %s\n") % value % name;
}

template <typename T, typename... Args>
boost::format format_tensor(int idx1, T idx2, Args... args)
{
   return boost::format(" %8d%s") % idx1 % format_tensor(idx2, args...);
}

inline boost::format format_matrix(
   int idx1, int idx2, double value, const std::string& name)
{
   return boost::format(" %2d %2d   %16.8E   # %s\n")
      % idx1 % idx2 % value % name;
}

inline boost::format format_vector(
   int idx, double value, const std::string& name)
{
   return boost::format(" %5d   %16.8E   # %s\n") % idx % value % name;
}

inline boost::format format_element(
   int idx, double value, const std::string& name)
{
   return boost::format(" %5d   %16.8E   # %s\n") % idx % value % name;
}

template <typename Container>
boost::format format_decay(double br, const Container& pids, const std::string& name)
{
   const std::size_t nda = pids.size();

   boost::format formatted_ids("");
   for (std::size_t i = nda - 1; i >= 0; --i) {
      formatted_ids = boost::format(" %9d%s") % pids[i] % formatted_ids;
   }

   return boost::format("   %16.8E  %2d  %s  # %s\n") % br % nda % formatted_ids % name;
}

inline boost::format format_total_width(int pdg, double width, const std::string& name)
{
   return boost::format("%9d   %16.8E   # %s\n") % pdg % width % name;
}

inline boost::format format_number(double value, const std::string& name)
{
   return boost::format("         %16.8E   # %s\n") % value % name;
}

inline boost::format format_scale(double value)
{
   return boost::format("%9.8E") % value;
}

inline boost::format format_spinfo(int idx, const std::string& entry)
{
   return boost::format(" %5d   %s\n") % idx % entry;
}

} // namespace flexiblesusy

#endif
