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

#include "build_info.hpp"
#include "config.h"

#include <boost/version.hpp>
#include <Eigen/Core>
#include <gsl/gsl_version.h>

#include <iostream>

namespace flexiblesusy {

void print_all_info(std::ostream& ostr)
{
   ostr <<
      "Version information\n"
      "===================\n\n";
   print_version_info(ostr);

   ostr <<
      "\n"
      "System information\n"
      "==================\n\n";
   print_system_info(ostr);

   ostr <<
      "\n"
      "Build information\n"
      "=================\n\n";
   print_build_info(ostr);
}

void print_flexiblesusy_version(std::ostream& ostr)
{
   ostr << FLEXIBLESUSY_VERSION;
}

void print_version_info(std::ostream& ostr)
{
   const int boost_major = BOOST_VERSION / 100000;
   const int boost_minor = (BOOST_VERSION / 100) % 1000;
   const int boost_patch = BOOST_VERSION % 100;

   ostr <<
      "FlexibleSUSY version:                   " FLEXIBLESUSY_VERSION "\n"
      "FlexibleSUSY git commit:                " GIT_COMMIT "\n"
      "SARAH version:                          " SARAH_VERSION "\n"
      "Mathematica version:                    " << MATHEMATICA_VERSION << "\n"
      "Boost version:                          "
        << boost_major << '.' << boost_minor << '.' << boost_patch << "\n"
      "Eigen version:                          " << EIGEN_WORLD_VERSION
        << '.' << EIGEN_MAJOR_VERSION << '.' << EIGEN_MINOR_VERSION << "\n"
      "GSL version:                            " << GSL_VERSION "\n"
      "GM2Calc version:                        " << GM2CALC_VERSION "\n"
      "Himalaya version:                       " << HIMALAYA_VERSION "\n"
      ;
}

void print_build_info(std::ostream& ostr)
{
   ostr <<
      "Two-scale solver:                       "
#ifdef ENABLE_TWO_SCALE_SOLVER
      "enabled"
#else
      "disabled"
#endif
      "\n"
      "Semi-analytic solver:                   "
#ifdef ENABLE_SEMI_ANALYTIC_SOLVER
      "enabled"
#else
      "disabled"
#endif
      "\n"
      "Colored output:                         "
#ifdef ENABLE_COLORS
      "yes"
#else
      "no"
#endif
      "\n"
      "Debug output:                           "
#ifdef ENABLE_DEBUG
      "yes"
#else
      "no"
#endif
      "\n"
      "Silent output:                          "
#ifdef ENABLE_SILENT
      "yes"
#else
      "no"
#endif
      "\n"
      "Verbose output:                         "
#ifdef ENABLE_VERBOSE
      "yes"
#else
      "no"
#endif
      "\n"
      "Use Boost.Numeric.Odeint:               "
#ifdef ENABLE_ODEINT
      "yes"
#else
      "no"
#endif
      "\n"
      "Use FFlite:                             "
#ifdef ENABLE_FFLITE
      "yes"
#else
      "no"
#endif
      "\n"
      "Use GM2Calc:                            "
#ifdef ENABLE_GM2Calc
      "yes"
#else
      "no"
#endif
      "\n"
      "Use Himalaya:                           "
#ifdef ENABLE_HIMALAYA
      "yes"
#else
      "no"
#endif
      "\n"
      "Use LAPACK:                             "
#ifdef ENABLE_LAPACK
      "yes"
#else
      "no"
#endif
      "\n"
      "Use LibraryLink:                        "
#ifdef ENABLE_LIBRARYLINK
      "yes"
#else
      "no"
#endif
      "\n"
      "Use LoopTools:                          "
#ifdef ENABLE_LOOPTOOLS
      "yes"
#else
      "no"
#endif
      "\n"
      "Use MKL ILP64 workaround:               "
#ifdef ENABLE_ILP64MKL_WORKAROUND
      "yes"
#else
      "no"
#endif
      "\n"
      "Use multi-threading:                    "
#ifdef ENABLE_THREADS
      "yes"
#else
      "no"
#endif
      "\n"
      "Use SQLite3:                            "
#ifdef ENABLE_SQLITE
      "yes"
#else
      "no"
#endif
      "\n"
      "Mass eigenvalue error check:            "
#ifdef CHECK_EIGENVALUE_ERROR
      "yes"
#else
      "no"
#endif
      "\n"
      ;
}

void print_system_info(std::ostream& ostr)
{
   ostr <<
      "Operating system (uname -s):            " OPERATING_SYSTEM "\n"
      "Kernel version (uname -r):              " KERNEL_VERSION "\n"
      ;
}

} // namespace flexiblesusy
