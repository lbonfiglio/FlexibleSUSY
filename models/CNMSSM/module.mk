DIR          := models/CNMSSM
MODNAME      := CNMSSM
SARAH_MODEL  := NMSSM
WITH_$(MODNAME) := yes
MODCNMSSM_MOD := SM MSSM_higgs NMSSM_higgs
MODCNMSSM_DEP := $(patsubst %,model_specific/%,$(MODCNMSSM_MOD))
MODCNMSSM_INC := $(patsubst %,-Imodel_specific/%,$(MODCNMSSM_MOD))
MODCNMSSM_LIB := $(foreach M,$(MODCNMSSM_MOD),model_specific/$M/libmodel_specific_$M$(MODULE_LIBEXT))

MODCNMSSM_SUBMOD  := $(DIR)/cxx_qft
MODCNMSSM_SUBMOD_INC := $(patsubst %,-I%,$(MODCNMSSM_SUBMOD))

CNMSSM_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)
CNMSSM_INSTALL_CXXQFT_DIR := \
		$(CNMSSM_INSTALL_DIR)/cxx_qft

CNMSSM_MK     := \
		$(DIR)/module.mk

CNMSSM_SUSY_BETAS_MK := \
		$(DIR)/susy_betas.mk

CNMSSM_SOFT_BETAS_MK := \
		$(DIR)/soft_betas.mk

CNMSSM_CXX_QFT_VERTICES_MK := \
		$(DIR)/cxx_qft/vertices.mk

CNMSSM_FlexibleEFTHiggs_MK := \
		$(DIR)/FlexibleEFTHiggs.mk

CNMSSM_INCLUDE_MK := \
		$(CNMSSM_SUSY_BETAS_MK) \
		$(CNMSSM_SOFT_BETAS_MK) \
		$(CNMSSM_CXX_QFT_VERTICES_MK)

CNMSSM_SLHA_INPUT := \
		$(DIR)/LesHouches.in.CNMSSM_generated \
		$(DIR)/LesHouches.in.CNMSSM

CNMSSM_REFERENCES := \
		$(DIR)/CNMSSM_references.tex

CNMSSM_GNUPLOT := \
		$(DIR)/CNMSSM_plot_rgflow.gnuplot \
		$(DIR)/CNMSSM_plot_spectrum.gnuplot

CNMSSM_TARBALL := \
		$(MODNAME).tar.gz

LIBCNMSSM_SRC := \
		$(DIR)/CNMSSM_a_muon.cpp \
		$(DIR)/CNMSSM_edm.cpp \
		$(DIR)/CNMSSM_effective_couplings.cpp \
		$(DIR)/CNMSSM_info.cpp \
		$(DIR)/CNMSSM_input_parameters.cpp \
		$(DIR)/CNMSSM_mass_eigenstates.cpp \
		$(DIR)/CNMSSM_observables.cpp \
		$(DIR)/CNMSSM_physical.cpp \
		$(DIR)/CNMSSM_slha_io.cpp \
		$(DIR)/CNMSSM_soft_parameters.cpp \
		$(DIR)/CNMSSM_susy_parameters.cpp \
		$(DIR)/CNMSSM_utilities.cpp \
		$(DIR)/CNMSSM_weinberg_angle.cpp

EXECNMSSM_SRC := \
		$(DIR)/run_CNMSSM.cpp \
		$(DIR)/run_cmd_line_CNMSSM.cpp \
		$(DIR)/scan_CNMSSM.cpp
LLCNMSSM_LIB  :=
LLCNMSSM_OBJ  :=
LLCNMSSM_SRC  := \
		$(DIR)/CNMSSM_librarylink.cpp

LLCNMSSM_MMA  := \
		$(DIR)/CNMSSM_librarylink.m \
		$(DIR)/run_CNMSSM.m

LIBCNMSSM_HDR := \
		$(DIR)/CNMSSM_a_muon.hpp \
		$(DIR)/CNMSSM_convergence_tester.hpp \
		$(DIR)/CNMSSM_edm.hpp \
		$(DIR)/CNMSSM_effective_couplings.hpp \
		$(DIR)/CNMSSM_ewsb_solver.hpp \
		$(DIR)/CNMSSM_ewsb_solver_interface.hpp \
		$(DIR)/CNMSSM_high_scale_constraint.hpp \
		$(DIR)/CNMSSM_info.hpp \
		$(DIR)/CNMSSM_initial_guesser.hpp \
		$(DIR)/CNMSSM_input_parameters.hpp \
		$(DIR)/CNMSSM_low_scale_constraint.hpp \
		$(DIR)/CNMSSM_mass_eigenstates.hpp \
		$(DIR)/CNMSSM_model.hpp \
		$(DIR)/CNMSSM_model_slha.hpp \
		$(DIR)/CNMSSM_observables.hpp \
		$(DIR)/CNMSSM_physical.hpp \
		$(DIR)/CNMSSM_slha_io.hpp \
		$(DIR)/CNMSSM_spectrum_generator.hpp \
		$(DIR)/CNMSSM_spectrum_generator_interface.hpp \
		$(DIR)/CNMSSM_soft_parameters.hpp \
		$(DIR)/CNMSSM_susy_parameters.hpp \
		$(DIR)/CNMSSM_susy_scale_constraint.hpp \
		$(DIR)/CNMSSM_utilities.hpp \
		$(DIR)/CNMSSM_weinberg_angle.hpp

LIBCNMSSM_CXXQFT_HDR := \
		$(DIR)/cxx_qft/CNMSSM_qft.hpp \
		$(DIR)/cxx_qft/CNMSSM_fields.hpp \
		$(DIR)/cxx_qft/CNMSSM_vertices.hpp \
		$(DIR)/cxx_qft/CNMSSM_context_base.hpp \
		$(DIR)/cxx_qft/CNMSSM_npointfunctions.hpp

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
-include $(CNMSSM_SUSY_BETAS_MK)
-include $(CNMSSM_SOFT_BETAS_MK)
-include $(CNMSSM_CXX_QFT_VERTICES_MK)
-include $(CNMSSM_FlexibleEFTHiggs_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(CNMSSM_SUSY_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(CNMSSM_SOFT_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(CNMSSM_CXX_QFT_VERTICES_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(CNMSSM_FlexibleEFTHiggs_MK): run-metacode-$(MODNAME)
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
LIBCNMSSM_SRC := $(sort $(LIBCNMSSM_SRC))
EXECNMSSM_SRC := $(sort $(EXECNMSSM_SRC))

LIBCNMSSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBCNMSSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBCNMSSM_SRC)))

EXECNMSSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXECNMSSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXECNMSSM_SRC)))

EXECNMSSM_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXECNMSSM_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXECNMSSM_SRC)))

LIBCNMSSM_DEP := \
		$(LIBCNMSSM_OBJ:.o=.d)

EXECNMSSM_DEP := \
		$(EXECNMSSM_OBJ:.o=.d)

LLCNMSSM_DEP  := \
		$(patsubst %.cpp, %.d, $(filter %.cpp, $(LLCNMSSM_SRC)))

LLCNMSSM_OBJ  := $(LLCNMSSM_SRC:.cpp=.o)
LLCNMSSM_LIB  := $(LLCNMSSM_SRC:.cpp=$(LIBLNK_LIBEXT))

LIBCNMSSM     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

METACODE_STAMP_CNMSSM := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_CNMSSM := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL_DIR) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib \
		clean-$(MODNAME)-obj distclean-$(MODNAME) \
		run-metacode-$(MODNAME) pack-$(MODNAME)-src

all-$(MODNAME): $(LIBCNMSSM) $(EXECNMSSM_EXE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(CNMSSM_INSTALL_DIR)
		install -d $(CNMSSM_INSTALL_CXXQFT_DIR)
		install -m u=rw,g=r,o=r $(LIBCNMSSM_SRC) $(CNMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBCNMSSM_HDR) $(CNMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBCNMSSM_CXXQFT_HDR) $(CNMSSM_INSTALL_CXXQFT_DIR)
		install -m u=rw,g=r,o=r $(EXECNMSSM_SRC) $(CNMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLCNMSSM_SRC) $(CNMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLCNMSSM_MMA) $(CNMSSM_INSTALL_DIR)
		$(INSTALL_STRIPPED) $(CNMSSM_MK) $(CNMSSM_INSTALL_DIR) -m u=rw,g=r,o=r
		install -m u=rw,g=r,o=r $(CNMSSM_INCLUDE_MK) $(CNMSSM_INSTALL_DIR)
ifneq ($(CNMSSM_SLHA_INPUT),)
		install -m u=rw,g=r,o=r $(CNMSSM_SLHA_INPUT) $(CNMSSM_INSTALL_DIR)
endif
		install -m u=rw,g=r,o=r $(CNMSSM_REFERENCES) $(CNMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(CNMSSM_GNUPLOT) $(CNMSSM_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIBCNMSSM_DEP)
		-rm -f $(EXECNMSSM_DEP)
		-rm -f $(LLCNMSSM_DEP)

clean-$(MODNAME)-lib:
		-rm -f $(LIBCNMSSM)
		-rm -f $(LLCNMSSM_LIB)

clean-$(MODNAME)-obj:
		-rm -f $(LIBCNMSSM_OBJ)
		-rm -f $(EXECNMSSM_OBJ)
		-rm -f $(LLCNMSSM_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		-rm -f $(LIBCNMSSM_SRC)
		-rm -f $(LIBCNMSSM_HDR)
		-rm -f $(LIBCNMSSM_CXXQFT_HDR)
		-rm -f $(EXECNMSSM_SRC)
		-rm -f $(LLCNMSSM_SRC)
		-rm -f $(LLCNMSSM_MMA)
		-rm -f $(METACODE_STAMP_CNMSSM)
		-rm -f $(CNMSSM_INCLUDE_MK)
		-rm -f $(CNMSSM_SLHA_INPUT)
		-rm -f $(CNMSSM_REFERENCES)
		-rm -f $(CNMSSM_GNUPLOT)

distclean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		-rm -f $(EXECNMSSM_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-generated:: clean-$(MODNAME)-src

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		tar -czf $(CNMSSM_TARBALL) \
		$(LIBCNMSSM_SRC) $(LIBCNMSSM_HDR) $(LIBCNMSSM_CXXQFT_HDR) \
		$(EXECNMSSM_SRC) \
		$(LLCNMSSM_SRC) $(LLCNMSSM_MMA) \
		$(CNMSSM_MK) $(CNMSSM_INCLUDE_MK) \
		$(CNMSSM_SLHA_INPUT) $(CNMSSM_REFERENCES) \
		$(CNMSSM_GNUPLOT)

$(LIBCNMSSM_SRC) $(LIBCNMSSM_HDR) $(LIBCNMSSM_CXXQFT_HDR) $(EXECNMSSM_SRC) $(LLCNMSSM_SRC) $(LLCNMSSM_MMA) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_CNMSSM)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_CNMSSM): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_CNMSSM)
		"$(MATH)" -run "Get[\"$<\"]; Quit[]" || (echo "Error: The code generation failed!"; exit 1)
		@touch "$(METACODE_STAMP_CNMSSM)"
		@echo "Note: to regenerate CNMSSM source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_CNMSSM)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_CNMSSM):
		@true
endif

$(LIBCNMSSM_DEP) $(EXECNMSSM_DEP) $(LLCNMSSM_DEP) $(LIBCNMSSM_OBJ) $(EXECNMSSM_OBJ) $(LLCNMSSM_OBJ) $(LLCNMSSM_LIB): \
	CPPFLAGS += $(MODCNMSSM_SUBMOD_INC) $(MODCNMSSM_INC) $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(TSILFLAGS) $(HIMALAYAFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBCNMSSM_DEP) $(EXECNMSSM_DEP) $(LLCNMSSM_DEP) $(LIBCNMSSM_OBJ) $(EXECNMSSM_OBJ) $(LLCNMSSM_OBJ) $(LLCNMSSM_LIB): \
	CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LLCNMSSM_OBJ) $(LLCNMSSM_LIB): \
	CPPFLAGS += $(LLFLAGS)

$(LIBCNMSSM): $(LIBCNMSSM_OBJ)
		$(MODULE_MAKE_LIB_CMD) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBCNMSSM) $(MODCNMSSM_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS)

$(LLCNMSSM_LIB): $(LLCNMSSM_OBJ) $(LIBCNMSSM) $(MODCNMSSM_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(LIBLNK_MAKE_LIB_CMD) $@ $(CPPFLAGS) $(CFLAGS) $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS) $(LLLIBS)

ALLDEP += $(LIBCNMSSM_DEP) $(EXECNMSSM_DEP)
ALLSRC += $(LIBCNMSSM_SRC) $(EXECNMSSM_SRC)
ALLLIB += $(LIBCNMSSM)
ALLEXE += $(EXECNMSSM_EXE)
ALLMODDEP += $(MODCNMSSM_DEP)

ifeq ($(ENABLE_LIBRARYLINK),yes)
ALLDEP += $(LLCNMSSM_DEP)
ALLSRC += $(LLCNMSSM_SRC)
ALLLL  += $(LLCNMSSM_LIB)
endif
