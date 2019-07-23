DIR          := models/NMSSM
MODNAME      := NMSSM
SARAH_MODEL  := NMSSM
WITH_$(MODNAME) := yes
MODNMSSM_MOD := SM MSSM_higgs NMSSM_higgs
MODNMSSM_DEP := $(patsubst %,model_specific/%,$(MODNMSSM_MOD))
MODNMSSM_INC := $(patsubst %,-Imodel_specific/%,$(MODNMSSM_MOD))
MODNMSSM_LIB := $(foreach M,$(MODNMSSM_MOD),model_specific/$M/libmodel_specific_$M$(MODULE_LIBEXT))

MODNMSSM_SUBMOD  := $(DIR)/cxx_qft
MODNMSSM_SUBMOD_INC := $(patsubst %,-I%,$(MODNMSSM_SUBMOD))

NMSSM_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)
NMSSM_INSTALL_CXXQFT_DIR := \
		$(NMSSM_INSTALL_DIR)/cxx_qft

NMSSM_MK     := \
		$(DIR)/module.mk

NMSSM_SUSY_BETAS_MK := \
		$(DIR)/susy_betas.mk

NMSSM_SOFT_BETAS_MK := \
		$(DIR)/soft_betas.mk

NMSSM_CXX_QFT_VERTICES_MK := \
		$(DIR)/cxx_qft/vertices.mk

NMSSM_FlexibleEFTHiggs_MK := \
		$(DIR)/FlexibleEFTHiggs.mk

NMSSM_INCLUDE_MK := \
		$(NMSSM_SUSY_BETAS_MK) \
		$(NMSSM_SOFT_BETAS_MK) \
		$(NMSSM_CXX_QFT_VERTICES_MK)

NMSSM_SLHA_INPUT := \
		$(DIR)/LesHouches.in.NMSSM_generated \


NMSSM_REFERENCES := \
		$(DIR)/NMSSM_references.tex

NMSSM_GNUPLOT := \
		$(DIR)/NMSSM_plot_rgflow.gnuplot \
		$(DIR)/NMSSM_plot_spectrum.gnuplot

NMSSM_TARBALL := \
		$(MODNAME).tar.gz

LIBNMSSM_SRC := \
		$(DIR)/NMSSM_a_muon.cpp \
		$(DIR)/NMSSM_edm.cpp \
		$(DIR)/NMSSM_effective_couplings.cpp \
		$(DIR)/NMSSM_info.cpp \
		$(DIR)/NMSSM_input_parameters.cpp \
		$(DIR)/NMSSM_mass_eigenstates.cpp \
		$(DIR)/NMSSM_observables.cpp \
		$(DIR)/NMSSM_physical.cpp \
		$(DIR)/NMSSM_slha_io.cpp \
		$(DIR)/NMSSM_soft_parameters.cpp \
		$(DIR)/NMSSM_susy_parameters.cpp \
		$(DIR)/NMSSM_utilities.cpp \
		$(DIR)/NMSSM_weinberg_angle.cpp

EXENMSSM_SRC := \
		$(DIR)/run_NMSSM.cpp \
		$(DIR)/run_cmd_line_NMSSM.cpp \
		$(DIR)/scan_NMSSM.cpp
LLNMSSM_LIB  :=
LLNMSSM_OBJ  :=
LLNMSSM_SRC  := \
		$(DIR)/NMSSM_librarylink.cpp

LLNMSSM_MMA  := \
		$(DIR)/NMSSM_librarylink.m \
		$(DIR)/run_NMSSM.m

LIBNMSSM_HDR := \
		$(DIR)/NMSSM_a_muon.hpp \
		$(DIR)/NMSSM_convergence_tester.hpp \
		$(DIR)/NMSSM_edm.hpp \
		$(DIR)/NMSSM_effective_couplings.hpp \
		$(DIR)/NMSSM_ewsb_solver.hpp \
		$(DIR)/NMSSM_ewsb_solver_interface.hpp \
		$(DIR)/NMSSM_high_scale_constraint.hpp \
		$(DIR)/NMSSM_info.hpp \
		$(DIR)/NMSSM_initial_guesser.hpp \
		$(DIR)/NMSSM_input_parameters.hpp \
		$(DIR)/NMSSM_low_scale_constraint.hpp \
		$(DIR)/NMSSM_mass_eigenstates.hpp \
		$(DIR)/NMSSM_model.hpp \
		$(DIR)/NMSSM_model_slha.hpp \
		$(DIR)/NMSSM_observables.hpp \
		$(DIR)/NMSSM_physical.hpp \
		$(DIR)/NMSSM_slha_io.hpp \
		$(DIR)/NMSSM_spectrum_generator.hpp \
		$(DIR)/NMSSM_spectrum_generator_interface.hpp \
		$(DIR)/NMSSM_soft_parameters.hpp \
		$(DIR)/NMSSM_susy_parameters.hpp \
		$(DIR)/NMSSM_susy_scale_constraint.hpp \
		$(DIR)/NMSSM_utilities.hpp \
		$(DIR)/NMSSM_weinberg_angle.hpp

LIBNMSSM_CXXQFT_HDR := \
		$(DIR)/cxx_qft/NMSSM_qft.hpp \
		$(DIR)/cxx_qft/NMSSM_fields.hpp \
		$(DIR)/cxx_qft/NMSSM_vertices.hpp \
		$(DIR)/cxx_qft/NMSSM_context_base.hpp \
		$(DIR)/cxx_qft/NMSSM_npointfunctions.hpp

ifneq ($(findstring two_scale,$(SOLVERS)),)
-include $(DIR)/two_scale.mk
endif
ifneq ($(findstring lattice,$(SOLVERS)),)
-include $(DIR)/lattice.mk
endif
ifneq ($(findstring semi_analytic,$(SOLVERS)),)
-include $(DIR)/semi_analytic.mk
endif

ifneq ($(MAKECMDGOALS),showbuild)
ifneq ($(MAKECMDGOALS),tag)
ifneq ($(MAKECMDGOALS),release)
ifneq ($(MAKECMDGOALS),doc)
-include $(NMSSM_SUSY_BETAS_MK)
-include $(NMSSM_SOFT_BETAS_MK)
-include $(NMSSM_CXX_QFT_VERTICES_MK)
-include $(NMSSM_FlexibleEFTHiggs_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(NMSSM_SUSY_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(NMSSM_SOFT_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(NMSSM_CXX_QFT_VERTICES_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(NMSSM_FlexibleEFTHiggs_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
endif
endif
endif
endif
endif
endif
endif
endif
endif
endif

# remove duplicates in case all solvers are used
LIBNMSSM_SRC := $(sort $(LIBNMSSM_SRC))
EXENMSSM_SRC := $(sort $(EXENMSSM_SRC))

LIBNMSSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBNMSSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBNMSSM_SRC)))

EXENMSSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXENMSSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXENMSSM_SRC)))

EXENMSSM_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXENMSSM_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXENMSSM_SRC)))

LIBNMSSM_DEP := \
		$(LIBNMSSM_OBJ:.o=.d)

EXENMSSM_DEP := \
		$(EXENMSSM_OBJ:.o=.d)

LLNMSSM_DEP  := \
		$(patsubst %.cpp, %.d, $(filter %.cpp, $(LLNMSSM_SRC)))

LLNMSSM_OBJ  := $(LLNMSSM_SRC:.cpp=.o)
LLNMSSM_LIB  := $(LLNMSSM_SRC:.cpp=$(LIBLNK_LIBEXT))

LIBNMSSM     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

METACODE_STAMP_NMSSM := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_NMSSM := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL_DIR) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib \
		clean-$(MODNAME)-obj distclean-$(MODNAME) \
		run-metacode-$(MODNAME) pack-$(MODNAME)-src

all-$(MODNAME): $(LIBNMSSM) $(EXENMSSM_EXE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(NMSSM_INSTALL_DIR)
		install -d $(NMSSM_INSTALL_CXXQFT_DIR)
		install -m u=rw,g=r,o=r $(LIBNMSSM_SRC) $(NMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBNMSSM_HDR) $(NMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBNMSSM_CXXQFT_HDR) $(NMSSM_INSTALL_CXXQFT_DIR)
		install -m u=rw,g=r,o=r $(EXENMSSM_SRC) $(NMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLNMSSM_SRC) $(NMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLNMSSM_MMA) $(NMSSM_INSTALL_DIR)
		$(INSTALL_STRIPPED) $(NMSSM_MK) $(NMSSM_INSTALL_DIR) -m u=rw,g=r,o=r
		install -m u=rw,g=r,o=r $(NMSSM_INCLUDE_MK) $(NMSSM_INSTALL_DIR)
ifneq ($(NMSSM_SLHA_INPUT),)
		install -m u=rw,g=r,o=r $(NMSSM_SLHA_INPUT) $(NMSSM_INSTALL_DIR)
endif
		install -m u=rw,g=r,o=r $(NMSSM_REFERENCES) $(NMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(NMSSM_GNUPLOT) $(NMSSM_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIBNMSSM_DEP)
		-rm -f $(EXENMSSM_DEP)
		-rm -f $(LLNMSSM_DEP)

clean-$(MODNAME)-lib:
		-rm -f $(LIBNMSSM)
		-rm -f $(LLNMSSM_LIB)

clean-$(MODNAME)-obj:
		-rm -f $(LIBNMSSM_OBJ)
		-rm -f $(EXENMSSM_OBJ)
		-rm -f $(LLNMSSM_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		-rm -f $(LIBNMSSM_SRC)
		-rm -f $(LIBNMSSM_HDR)
		-rm -f $(LIBNMSSM_CXXQFT_HDR)
		-rm -f $(EXENMSSM_SRC)
		-rm -f $(LLNMSSM_SRC)
		-rm -f $(LLNMSSM_MMA)
		-rm -f $(METACODE_STAMP_NMSSM)
		-rm -f $(NMSSM_INCLUDE_MK)
		-rm -f $(NMSSM_SLHA_INPUT)
		-rm -f $(NMSSM_REFERENCES)
		-rm -f $(NMSSM_GNUPLOT)

distclean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		-rm -f $(EXENMSSM_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-generated:: clean-$(MODNAME)-src

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		tar -czf $(NMSSM_TARBALL) \
		$(LIBNMSSM_SRC) $(LIBNMSSM_HDR) $(LIBNMSSM_CXXQFT_HDR) \
		$(EXENMSSM_SRC) \
		$(LLNMSSM_SRC) $(LLNMSSM_MMA) \
		$(NMSSM_MK) $(NMSSM_INCLUDE_MK) \
		$(NMSSM_SLHA_INPUT) $(NMSSM_REFERENCES) \
		$(NMSSM_GNUPLOT)

$(LIBNMSSM_SRC) $(LIBNMSSM_HDR) $(LIBNMSSM_CXXQFT_HDR) $(EXENMSSM_SRC) $(LLNMSSM_SRC) $(LLNMSSM_MMA) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_NMSSM)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_NMSSM): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_NMSSM)
		"$(MATH)" -run "Get[\"$<\"]; Quit[]" || (echo "Error: The code generation failed!"; exit 1)
		@touch "$(METACODE_STAMP_NMSSM)"
		@echo "Note: to regenerate NMSSM source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_NMSSM)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_NMSSM):
		@true
endif

$(LIBNMSSM_DEP) $(EXENMSSM_DEP) $(LLNMSSM_DEP) $(LIBNMSSM_OBJ) $(EXENMSSM_OBJ) $(LLNMSSM_OBJ) $(LLNMSSM_LIB): \
	CPPFLAGS += $(MODNMSSM_SUBMOD_INC) $(MODNMSSM_INC) $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(TSILFLAGS) $(HIMALAYAFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBNMSSM_DEP) $(EXENMSSM_DEP) $(LLNMSSM_DEP) $(LIBNMSSM_OBJ) $(EXENMSSM_OBJ) $(LLNMSSM_OBJ) $(LLNMSSM_LIB): \
	CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LLNMSSM_OBJ) $(LLNMSSM_LIB): \
	CPPFLAGS += $(LLFLAGS)

$(LIBNMSSM): $(LIBNMSSM_OBJ)
		$(MODULE_MAKE_LIB_CMD) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBNMSSM) $(MODNMSSM_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS)

$(LLNMSSM_LIB): $(LLNMSSM_OBJ) $(LIBNMSSM) $(MODNMSSM_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(LIBLNK_MAKE_LIB_CMD) $@ $(CPPFLAGS) $(CFLAGS) $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS) $(LLLIBS)

ALLDEP += $(LIBNMSSM_DEP) $(EXENMSSM_DEP)
ALLSRC += $(LIBNMSSM_SRC) $(EXENMSSM_SRC)
ALLLIB += $(LIBNMSSM)
ALLEXE += $(EXENMSSM_EXE)
ALLMODDEP += $(MODNMSSM_DEP)

ifeq ($(ENABLE_LIBRARYLINK),yes)
ALLDEP += $(LLNMSSM_DEP)
ALLSRC += $(LLNMSSM_SRC)
ALLLL  += $(LLNMSSM_LIB)
endif
