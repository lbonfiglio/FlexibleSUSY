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

#ifndef TWO_BODY_DECAYS_H
#define TWO_BODY_DECAYS_H

#include <complex>

/**
 * @file two_body_decays.hpp
 * @brief contains functions for calculating two-body decay widths
 */

namespace flexiblesusy {

namespace two_body_decays {

double kallen(double x, double y, double z);

double scalar_to_fermion_fermion(double ms, double mf1, double mf2,
                                 std::complex<double> cpl_left,
                                 std::complex<double> cpl_right);

} // namespace two_body_decays

} // namespace flexiblesusy

#endif
