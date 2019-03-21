DIR          := models/CNMSSMEFTHiggs
MODNAME      := CNMSSMEFTHiggs
SARAH_MODEL  := NMSSM
WITH_$(MODNAME) := yes
MODCNMSSMEFTHiggs_MOD := SM MSSM_higgs NMSSM_higgs
MODCNMSSMEFTHiggs_DEP := $(patsubst %,model_specific/%,$(MODCNMSSMEFTHiggs_MOD))
MODCNMSSMEFTHiggs_INC := $(patsubst %,-Imodel_specific/%,$(MODCNMSSMEFTHiggs_MOD))
MODCNMSSMEFTHiggs_LIB := $(foreach M,$(MODCNMSSMEFTHiggs_MOD),model_specific/$M/libmodel_specific_$M$(MODULE_LIBEXT))

MODCNMSSMEFTHiggs_SUBMOD  := $(DIR)/cxx_qft
MODCNMSSMEFTHiggs_SUBMOD_INC := $(patsubst %,-I%,$(MODCNMSSMEFTHiggs_SUBMOD))

CNMSSMEFTHiggs_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)
CNMSSMEFTHiggs_INSTALL_CXXQFT_DIR := \
		$(CNMSSMEFTHiggs_INSTALL_DIR)/cxx_qft

CNMSSMEFTHiggs_MK     := \
		$(DIR)/module.mk

CNMSSMEFTHiggs_SUSY_BETAS_MK := \
		$(DIR)/susy_betas.mk

CNMSSMEFTHiggs_SOFT_BETAS_MK := \
		$(DIR)/soft_betas.mk

CNMSSMEFTHiggs_CXX_QFT_VERTICES_MK := \
		$(DIR)/cxx_qft/vertices.mk

CNMSSMEFTHiggs_FlexibleEFTHiggs_MK := \
		$(DIR)/FlexibleEFTHiggs.mk

CNMSSMEFTHiggs_INCLUDE_MK := \
		$(CNMSSMEFTHiggs_SUSY_BETAS_MK) \
		$(CNMSSMEFTHiggs_SOFT_BETAS_MK) \
		$(CNMSSMEFTHiggs_CXX_QFT_VERTICES_MK)

CNMSSMEFTHiggs_SLHA_INPUT := \
		$(DIR)/LesHouches.in.CNMSSMEFTHiggs_generated \


CNMSSMEFTHiggs_REFERENCES := \
		$(DIR)/CNMSSMEFTHiggs_references.tex

CNMSSMEFTHiggs_GNUPLOT := \
		$(DIR)/CNMSSMEFTHiggs_plot_rgflow.gnuplot \
		$(DIR)/CNMSSMEFTHiggs_plot_spectrum.gnuplot

CNMSSMEFTHiggs_TARBALL := \
		$(MODNAME).tar.gz

LIBCNMSSMEFTHiggs_SRC := \
		$(DIR)/CNMSSMEFTHiggs_a_muon.cpp \
		$(DIR)/CNMSSMEFTHiggs_edm.cpp \
		$(DIR)/CNMSSMEFTHiggs_effective_couplings.cpp \
		$(DIR)/CNMSSMEFTHiggs_info.cpp \
		$(DIR)/CNMSSMEFTHiggs_input_parameters.cpp \
		$(DIR)/CNMSSMEFTHiggs_mass_eigenstates.cpp \
		$(DIR)/CNMSSMEFTHiggs_observables.cpp \
		$(DIR)/CNMSSMEFTHiggs_physical.cpp \
		$(DIR)/CNMSSMEFTHiggs_slha_io.cpp \
		$(DIR)/CNMSSMEFTHiggs_soft_parameters.cpp \
		$(DIR)/CNMSSMEFTHiggs_susy_parameters.cpp \
		$(DIR)/CNMSSMEFTHiggs_utilities.cpp \
		$(DIR)/CNMSSMEFTHiggs_weinberg_angle.cpp

EXECNMSSMEFTHiggs_SRC := \
		$(DIR)/run_CNMSSMEFTHiggs.cpp \
		$(DIR)/run_cmd_line_CNMSSMEFTHiggs.cpp \
		$(DIR)/scan_CNMSSMEFTHiggs.cpp
LLCNMSSMEFTHiggs_LIB  :=
LLCNMSSMEFTHiggs_OBJ  :=
LLCNMSSMEFTHiggs_SRC  := \
		$(DIR)/CNMSSMEFTHiggs_librarylink.cpp

LLCNMSSMEFTHiggs_MMA  := \
		$(DIR)/CNMSSMEFTHiggs_librarylink.m \
		$(DIR)/run_CNMSSMEFTHiggs.m

LIBCNMSSMEFTHiggs_HDR := \
		$(DIR)/CNMSSMEFTHiggs_a_muon.hpp \
		$(DIR)/CNMSSMEFTHiggs_convergence_tester.hpp \
		$(DIR)/CNMSSMEFTHiggs_edm.hpp \
		$(DIR)/CNMSSMEFTHiggs_effective_couplings.hpp \
		$(DIR)/CNMSSMEFTHiggs_ewsb_solver.hpp \
		$(DIR)/CNMSSMEFTHiggs_ewsb_solver_interface.hpp \
		$(DIR)/CNMSSMEFTHiggs_high_scale_constraint.hpp \
		$(DIR)/CNMSSMEFTHiggs_info.hpp \
		$(DIR)/CNMSSMEFTHiggs_initial_guesser.hpp \
		$(DIR)/CNMSSMEFTHiggs_input_parameters.hpp \
		$(DIR)/CNMSSMEFTHiggs_low_scale_constraint.hpp \
		$(DIR)/CNMSSMEFTHiggs_mass_eigenstates.hpp \
		$(DIR)/CNMSSMEFTHiggs_model.hpp \
		$(DIR)/CNMSSMEFTHiggs_model_slha.hpp \
		$(DIR)/CNMSSMEFTHiggs_observables.hpp \
		$(DIR)/CNMSSMEFTHiggs_physical.hpp \
		$(DIR)/CNMSSMEFTHiggs_slha_io.hpp \
		$(DIR)/CNMSSMEFTHiggs_spectrum_generator.hpp \
		$(DIR)/CNMSSMEFTHiggs_spectrum_generator_interface.hpp \
		$(DIR)/CNMSSMEFTHiggs_soft_parameters.hpp \
		$(DIR)/CNMSSMEFTHiggs_susy_parameters.hpp \
		$(DIR)/CNMSSMEFTHiggs_susy_scale_constraint.hpp \
		$(DIR)/CNMSSMEFTHiggs_utilities.hpp \
		$(DIR)/CNMSSMEFTHiggs_weinberg_angle.hpp

LIBCNMSSMEFTHiggs_CXXQFT_HDR := \
		$(DIR)/cxx_qft/CNMSSMEFTHiggs_qft.hpp \
		$(DIR)/cxx_qft/CNMSSMEFTHiggs_fields.hpp \
		$(DIR)/cxx_qft/CNMSSMEFTHiggs_vertices.hpp \
		$(DIR)/cxx_qft/CNMSSMEFTHiggs_context_base.hpp \
		$(DIR)/cxx_qft/CNMSSMEFTHiggs_npointfunctions.hpp

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
-include $(CNMSSMEFTHiggs_SUSY_BETAS_MK)
-include $(CNMSSMEFTHiggs_SOFT_BETAS_MK)
-include $(CNMSSMEFTHiggs_CXX_QFT_VERTICES_MK)
-include $(CNMSSMEFTHiggs_FlexibleEFTHiggs_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(CNMSSMEFTHiggs_SUSY_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(CNMSSMEFTHiggs_SOFT_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(CNMSSMEFTHiggs_CXX_QFT_VERTICES_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(CNMSSMEFTHiggs_FlexibleEFTHiggs_MK): run-metacode-$(MODNAME)
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
LIBCNMSSMEFTHiggs_SRC := $(sort $(LIBCNMSSMEFTHiggs_SRC))
EXECNMSSMEFTHiggs_SRC := $(sort $(EXECNMSSMEFTHiggs_SRC))

LIBCNMSSMEFTHiggs_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBCNMSSMEFTHiggs_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBCNMSSMEFTHiggs_SRC)))

EXECNMSSMEFTHiggs_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXECNMSSMEFTHiggs_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXECNMSSMEFTHiggs_SRC)))

EXECNMSSMEFTHiggs_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXECNMSSMEFTHiggs_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXECNMSSMEFTHiggs_SRC)))

LIBCNMSSMEFTHiggs_DEP := \
		$(LIBCNMSSMEFTHiggs_OBJ:.o=.d)

EXECNMSSMEFTHiggs_DEP := \
		$(EXECNMSSMEFTHiggs_OBJ:.o=.d)

LLCNMSSMEFTHiggs_DEP  := \
		$(patsubst %.cpp, %.d, $(filter %.cpp, $(LLCNMSSMEFTHiggs_SRC)))

LLCNMSSMEFTHiggs_OBJ  := $(LLCNMSSMEFTHiggs_SRC:.cpp=.o)
LLCNMSSMEFTHiggs_LIB  := $(LLCNMSSMEFTHiggs_SRC:.cpp=$(LIBLNK_LIBEXT))

LIBCNMSSMEFTHiggs     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

METACODE_STAMP_CNMSSMEFTHiggs := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_CNMSSMEFTHiggs := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL_DIR) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib \
		clean-$(MODNAME)-obj distclean-$(MODNAME) \
		run-metacode-$(MODNAME) pack-$(MODNAME)-src

all-$(MODNAME): $(LIBCNMSSMEFTHiggs) $(EXECNMSSMEFTHiggs_EXE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(CNMSSMEFTHiggs_INSTALL_DIR)
		install -d $(CNMSSMEFTHiggs_INSTALL_CXXQFT_DIR)
		install -m u=rw,g=r,o=r $(LIBCNMSSMEFTHiggs_SRC) $(CNMSSMEFTHiggs_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBCNMSSMEFTHiggs_HDR) $(CNMSSMEFTHiggs_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBCNMSSMEFTHiggs_CXXQFT_HDR) $(CNMSSMEFTHiggs_INSTALL_CXXQFT_DIR)
		install -m u=rw,g=r,o=r $(EXECNMSSMEFTHiggs_SRC) $(CNMSSMEFTHiggs_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLCNMSSMEFTHiggs_SRC) $(CNMSSMEFTHiggs_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLCNMSSMEFTHiggs_MMA) $(CNMSSMEFTHiggs_INSTALL_DIR)
		$(INSTALL_STRIPPED) $(CNMSSMEFTHiggs_MK) $(CNMSSMEFTHiggs_INSTALL_DIR) -m u=rw,g=r,o=r
		install -m u=rw,g=r,o=r $(CNMSSMEFTHiggs_INCLUDE_MK) $(CNMSSMEFTHiggs_INSTALL_DIR)
ifneq ($(CNMSSMEFTHiggs_SLHA_INPUT),)
		install -m u=rw,g=r,o=r $(CNMSSMEFTHiggs_SLHA_INPUT) $(CNMSSMEFTHiggs_INSTALL_DIR)
endif
		install -m u=rw,g=r,o=r $(CNMSSMEFTHiggs_REFERENCES) $(CNMSSMEFTHiggs_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(CNMSSMEFTHiggs_GNUPLOT) $(CNMSSMEFTHiggs_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIBCNMSSMEFTHiggs_DEP)
		-rm -f $(EXECNMSSMEFTHiggs_DEP)
		-rm -f $(LLCNMSSMEFTHiggs_DEP)

clean-$(MODNAME)-lib:
		-rm -f $(LIBCNMSSMEFTHiggs)
		-rm -f $(LLCNMSSMEFTHiggs_LIB)

clean-$(MODNAME)-obj:
		-rm -f $(LIBCNMSSMEFTHiggs_OBJ)
		-rm -f $(EXECNMSSMEFTHiggs_OBJ)
		-rm -f $(LLCNMSSMEFTHiggs_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		-rm -f $(LIBCNMSSMEFTHiggs_SRC)
		-rm -f $(LIBCNMSSMEFTHiggs_HDR)
		-rm -f $(LIBCNMSSMEFTHiggs_CXXQFT_HDR)
		-rm -f $(EXECNMSSMEFTHiggs_SRC)
		-rm -f $(LLCNMSSMEFTHiggs_SRC)
		-rm -f $(LLCNMSSMEFTHiggs_MMA)
		-rm -f $(METACODE_STAMP_CNMSSMEFTHiggs)
		-rm -f $(CNMSSMEFTHiggs_INCLUDE_MK)
		-rm -f $(CNMSSMEFTHiggs_SLHA_INPUT)
		-rm -f $(CNMSSMEFTHiggs_REFERENCES)
		-rm -f $(CNMSSMEFTHiggs_GNUPLOT)

distclean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		-rm -f $(EXECNMSSMEFTHiggs_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-generated:: clean-$(MODNAME)-src

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		tar -czf $(CNMSSMEFTHiggs_TARBALL) \
		$(LIBCNMSSMEFTHiggs_SRC) $(LIBCNMSSMEFTHiggs_HDR) $(LIBCNMSSMEFTHiggs_CXXQFT_HDR) \
		$(EXECNMSSMEFTHiggs_SRC) \
		$(LLCNMSSMEFTHiggs_SRC) $(LLCNMSSMEFTHiggs_MMA) \
		$(CNMSSMEFTHiggs_MK) $(CNMSSMEFTHiggs_INCLUDE_MK) \
		$(CNMSSMEFTHiggs_SLHA_INPUT) $(CNMSSMEFTHiggs_REFERENCES) \
		$(CNMSSMEFTHiggs_GNUPLOT)

$(LIBCNMSSMEFTHiggs_SRC) $(LIBCNMSSMEFTHiggs_HDR) $(LIBCNMSSMEFTHiggs_CXXQFT_HDR) $(EXECNMSSMEFTHiggs_SRC) $(LLCNMSSMEFTHiggs_SRC) $(LLCNMSSMEFTHiggs_MMA) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_CNMSSMEFTHiggs)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_CNMSSMEFTHiggs): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_CNMSSMEFTHiggs)
		"$(MATH)" -run "Get[\"$<\"]; Quit[]" || (echo "Error: The code generation failed!"; exit 1)
		@touch "$(METACODE_STAMP_CNMSSMEFTHiggs)"
		@echo "Note: to regenerate CNMSSMEFTHiggs source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_CNMSSMEFTHiggs)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_CNMSSMEFTHiggs):
		@true
endif

$(LIBCNMSSMEFTHiggs_DEP) $(EXECNMSSMEFTHiggs_DEP) $(LLCNMSSMEFTHiggs_DEP) $(LIBCNMSSMEFTHiggs_OBJ) $(EXECNMSSMEFTHiggs_OBJ) $(LLCNMSSMEFTHiggs_OBJ) $(LLCNMSSMEFTHiggs_LIB): \
	CPPFLAGS += $(MODCNMSSMEFTHiggs_SUBMOD_INC) $(MODCNMSSMEFTHiggs_INC) $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(TSILFLAGS) $(HIMALAYAFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBCNMSSMEFTHiggs_DEP) $(EXECNMSSMEFTHiggs_DEP) $(LLCNMSSMEFTHiggs_DEP) $(LIBCNMSSMEFTHiggs_OBJ) $(EXECNMSSMEFTHiggs_OBJ) $(LLCNMSSMEFTHiggs_OBJ) $(LLCNMSSMEFTHiggs_LIB): \
	CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LLCNMSSMEFTHiggs_OBJ) $(LLCNMSSMEFTHiggs_LIB): \
	CPPFLAGS += $(LLFLAGS)

$(LIBCNMSSMEFTHiggs): $(LIBCNMSSMEFTHiggs_OBJ)
		$(MODULE_MAKE_LIB_CMD) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBCNMSSMEFTHiggs) $(MODCNMSSMEFTHiggs_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS)

$(LLCNMSSMEFTHiggs_LIB): $(LLCNMSSMEFTHiggs_OBJ) $(LIBCNMSSMEFTHiggs) $(MODCNMSSMEFTHiggs_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(LIBLNK_MAKE_LIB_CMD) $@ $(CPPFLAGS) $(CFLAGS) $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS) $(LLLIBS)

ALLDEP += $(LIBCNMSSMEFTHiggs_DEP) $(EXECNMSSMEFTHiggs_DEP)
ALLSRC += $(LIBCNMSSMEFTHiggs_SRC) $(EXECNMSSMEFTHiggs_SRC)
ALLLIB += $(LIBCNMSSMEFTHiggs)
ALLEXE += $(EXECNMSSMEFTHiggs_EXE)
ALLMODDEP += $(MODCNMSSMEFTHiggs_DEP)

ifeq ($(ENABLE_LIBRARYLINK),yes)
ALLDEP += $(LLCNMSSMEFTHiggs_DEP)
ALLSRC += $(LLCNMSSMEFTHiggs_SRC)
ALLLL  += $(LLCNMSSMEFTHiggs_LIB)
endif
