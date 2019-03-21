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

// File generated at Thu 21 Mar 2019 19:41:42

#ifndef CNMSSM_SEMI_ANALYTIC_CONVERGENCE_TESTER_H
#define CNMSSM_SEMI_ANALYTIC_CONVERGENCE_TESTER_H

#include "CNMSSM_convergence_tester.hpp"
#include "CNMSSM_semi_analytic_model.hpp"

#include "convergence_tester_drbar.hpp"

namespace flexiblesusy {

class Semi_analytic;

template<>
class CNMSSM_convergence_tester<Semi_analytic> : public Convergence_tester_DRbar<CNMSSM<Semi_analytic> > {
public:
   using Scale_getter = Convergence_tester_DRbar<CNMSSM<Semi_analytic>>::Scale_getter;

   CNMSSM_convergence_tester(CNMSSM<Semi_analytic>*, double, const Scale_getter& sg = Scale_getter());
   virtual ~CNMSSM_convergence_tester() = default;

protected:
   virtual double max_rel_diff() const;
};

} // namespace flexiblesusy

#endif
