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

// File generated at Tue 23 Jul 2019 15:16:39

/**
 * @file NMSSM_two_scale_model.hpp
 * @brief contains class for model with routines needed to solve boundary
 *        value problem using the two_scale solver by solving EWSB
 *        and determine the pole masses and mixings
 *
 * This file was generated at Tue 23 Jul 2019 15:16:39 with FlexibleSUSY
 * 2.3.0 (git commit: b132b3936f5e72e9f34d3f552dab7c81fb2c369e) and SARAH 4.12.3 .
 */

#ifndef NMSSM_TWO_SCALE_H
#define NMSSM_TWO_SCALE_H

#include "NMSSM_model.hpp"
#include "NMSSM_mass_eigenstates.hpp"

#include "model.hpp"

namespace flexiblesusy {

class Two_scale;
/**
 * @class NMSSM<Two_scale>
 * @brief model class with routines for determining masses and mixings and EWSB
 */
template<>
class NMSSM<Two_scale> : public Model, public NMSSM_mass_eigenstates {
public:
   explicit NMSSM(const NMSSM_input_parameters& input_ = NMSSM_input_parameters());
   NMSSM(const NMSSM&) = default;
   NMSSM(NMSSM&&) = default;
   virtual ~NMSSM() = default;
   NMSSM& operator=(const NMSSM&) = default;
   NMSSM& operator=(NMSSM&&) = default;

   // interface functions
   virtual void calculate_spectrum() override;
   virtual void clear_problems() override;
   virtual std::string name() const override;
   virtual void run_to(double scale, double eps = -1.0) override;
   virtual void print(std::ostream& out = std::cerr) const override;
   virtual void set_precision(double) override;
};

std::ostream& operator<<(std::ostream&, const NMSSM<Two_scale>&);

} // namespace flexiblesusy

#endif
