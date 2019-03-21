#  ====================================================================
#  This file is part of FlexibleSUSY.
#
#  FlexibleSUSY is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published
#  by the Free Software Foundation, either version 3 of the License,
#  or (at your option) any later version.
#
#  FlexibleSUSY is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with FlexibleSUSY.  If not, see
#  <http://www.gnu.org/licenses/>.
#  ====================================================================

CNMSSMEFTHiggs_INCLUDE_MK += $(DIR)/semi_analytic.mk

LIBCNMSSMEFTHiggs_SRC += \
		$(DIR)/CNMSSMEFTHiggs_semi_analytic_convergence_tester.cpp \
		$(DIR)/CNMSSMEFTHiggs_semi_analytic_ewsb_solver.cpp \
		$(DIR)/CNMSSMEFTHiggs_semi_analytic_high_scale_constraint.cpp \
		$(DIR)/CNMSSMEFTHiggs_semi_analytic_initial_guesser.cpp \
		$(DIR)/CNMSSMEFTHiggs_semi_analytic_low_scale_constraint.cpp \
		$(DIR)/CNMSSMEFTHiggs_semi_analytic_model.cpp \
		$(DIR)/CNMSSMEFTHiggs_semi_analytic_soft_parameters_constraint.cpp \
		$(DIR)/CNMSSMEFTHiggs_semi_analytic_solutions.cpp \
		$(DIR)/CNMSSMEFTHiggs_semi_analytic_spectrum_generator.cpp \
		$(DIR)/CNMSSMEFTHiggs_semi_analytic_susy_convergence_tester.cpp \
		$(DIR)/CNMSSMEFTHiggs_semi_analytic_susy_scale_constraint.cpp
LIBCNMSSMEFTHiggs_HDR += \
		$(DIR)/CNMSSMEFTHiggs_semi_analytic_convergence_tester.hpp \
		$(DIR)/CNMSSMEFTHiggs_semi_analytic_ewsb_solver.hpp \
		$(DIR)/CNMSSMEFTHiggs_semi_analytic_high_scale_constraint.hpp \
		$(DIR)/CNMSSMEFTHiggs_semi_analytic_initial_guesser.hpp \
		$(DIR)/CNMSSMEFTHiggs_semi_analytic_low_scale_constraint.hpp \
		$(DIR)/CNMSSMEFTHiggs_semi_analytic_model.hpp \
		$(DIR)/CNMSSMEFTHiggs_semi_analytic_soft_parameters_constraint.hpp \
		$(DIR)/CNMSSMEFTHiggs_semi_analytic_solutions.hpp \
		$(DIR)/CNMSSMEFTHiggs_semi_analytic_spectrum_generator.hpp \
		$(DIR)/CNMSSMEFTHiggs_semi_analytic_susy_convergence_tester.hpp \
		$(DIR)/CNMSSMEFTHiggs_semi_analytic_susy_scale_constraint.hpp \
		$(DIR)/CNMSSMEFTHiggs_soft_parameters_constraint.hpp \
		$(DIR)/CNMSSMEFTHiggs_susy_convergence_tester.hpp

LIBCNMSSMEFTHiggs_SRC += \
      		models/CNMSSMEFTHiggs/CNMSSMEFTHiggs_semi_analytic_matching_constraint.cpp \
		models/CNMSSMEFTHiggs/CNMSSMEFTHiggs_standard_model_semi_analytic_matching.cpp

LIBCNMSSMEFTHiggs_HDR += \
      		models/CNMSSMEFTHiggs/CNMSSMEFTHiggs_matching_constraint.hpp \
		models/CNMSSMEFTHiggs/CNMSSMEFTHiggs_semi_analytic_matching_constraint.hpp \
		models/CNMSSMEFTHiggs/CNMSSMEFTHiggs_standard_model_semi_analytic_matching.hpp
