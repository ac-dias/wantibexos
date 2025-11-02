include makefile.inc

FC_FLAGS := $(OMP) $(COND) $(EXTRA)
L_FLAGS  := $(OMP) $(LIBS) $(COND) $(EXTRA)

SRC_DIR   := ./src
BUILD_DIR := ./build
BIN_DIR   := $(DIR)
UTILS_DIR := ./utils

SUBROUTINE_NAMES := \
    berry_curvature_subs \
    boltzmann_subs \
    bse_subs \
    bse_subs_kpath \
    bse_subs_temp \
    special_funct \
    ei_spec_funct \
    coulomb_pot \
    diel-pp-subs \
    dos_subs \
    efmass-subs \
    general_subs \
    hamiltonians \
    hamiltonian_tb \
    mhkpack_subs \
    module_input_read \
    optics \
    spin_txt_subs \
    emission_subs \
    coulomb_pot_gw \
    gw_subs

SUBPROGRAM_NAMES := \
    bands-kpath-tool \
    berry_curvature_bz-tool \
    berry_curvature_kpath-tool \
    bse_diel-tool \
    bse_diel-tool-pol \
    bse_kpath-tool \
    bse_kpath-tool-temp \
    bse_solver-tool-diel \
    bse_solver-tool-diel-temp \
    diel-pp-bse \
    diel-pp-bse-pol \
    tdos-tool-stxt \
    diel-pp \
    diel-pp-pol \
    efmass \
    exciton_lifetime \
    sp_diel-tool \
    sp_diel-tool-pol \
    sp_opt_bz-tool \
    sp_solver-tool-diel \
    boltzmann_transport \
    emission_PL

SUBROUTINE_OBJS := $(foreach name,$(SUBROUTINE_NAMES),$(BUILD_DIR)/subroutines/$(name).o)
SUBPROGRAM_OBJS := $(foreach name,$(SUBPROGRAM_NAMES),$(BUILD_DIR)/subprograms/$(name).o)

MODULE_OBJECTS := $(SUBROUTINE_OBJS) $(SUBPROGRAM_OBJS)

MAIN_SRC  := $(SRC_DIR)/wtb_main.F90
MAIN_EXEC := $(BIN_DIR)/wtb.x

all: $(MAIN_EXEC) pp
	@echo "--- Build Complete ---"

$(MAIN_EXEC): $(MODULE_OBJECTS) $(MAIN_SRC) makefile.inc
	@mkdir -p $(BIN_DIR)
	@echo "--- Linking Main Executable: $@ ---"
	$(FOR) $(MAIN_SRC) $(MODULE_OBJECTS) -o $@ $(L_FLAGS)
	@cp $@ ./build/wtb.x

$(BUILD_DIR)/subroutines/%.o: $(SRC_DIR)/subroutines/%.F90 makefile.inc
	@mkdir -p $(dir $@)
	@echo "Compiling Subroutine: $< -> $@"
	$(FOR) -c $< -o $@ $(FC_FLAGS)

$(BUILD_DIR)/subprograms/%.o: $(SRC_DIR)/subprograms/%.F90 makefile.inc
	@mkdir -p $(dir $@)
	@echo "Compiling Subprogram: $< -> $@"
	$(FOR) -c $< -o $@ $(FC_FLAGS)

UTILS_EXECS := $(BIN_DIR)/nc_nv_finder.x \
               $(BIN_DIR)/param_gen.x \
               $(BIN_DIR)/param_gen_vasp.x \
               $(BIN_DIR)/absorbance.x \
               $(BIN_DIR)/pce.x \
               $(BIN_DIR)/huckel2wtb.x

pp: $(UTILS_EXECS)
	@echo "--- Copying Python Scripts ---"
	@cp $(UTILS_DIR)/*.py $(BIN_DIR)

$(BIN_DIR)/nc_nv_finder.x: $(UTILS_DIR)/nc_nv_finder.F90 makefile.inc
	$(FOR) $< -o $@ $(L_FLAGS)

$(BIN_DIR)/param_gen.x: $(UTILS_DIR)/param_gen.F90 makefile.inc
	$(FOR) $< -o $@ $(L_FLAGS)

$(BIN_DIR)/param_gen_vasp.x: $(UTILS_DIR)/param_gen_vasp.F90 makefile.inc
	$(FOR) $< -o $@ $(L_FLAGS)

$(BIN_DIR)/absorbance.x: $(UTILS_DIR)/absorbance.F90 makefile.inc
	$(FOR) $< -o $@ $(L_FLAGS)

$(BIN_DIR)/pce.x: $(UTILS_DIR)/slme/pce-code.f90 $(UTILS_DIR)/slme/pce-subs.f90 makefile.inc
	$(FOR) $^ -o $@ $(L_FLAGS)

$(BIN_DIR)/huckel2wtb.x: $(UTILS_DIR)/huckel2wtb/src/overlaps_jc.f90 $(UTILS_DIR)/huckel2wtb/src/diagonalize.f90 $(UTILS_DIR)/huckel2wtb/src/Huckel_TB.f90 makefile.inc
	$(FOR) $^ -o $@ $(L_FLAGS)

clean:
	@echo "--- Cleaning build, bin, and .mod files ---"
	@rm -rf $(BUILD_DIR) $(BIN_DIR) ./*.mod
	@mkdir -p $(BUILD_DIR) $(BIN_DIR)

$(BUILD_DIR)/subroutines/coulomb_pot.o: $(BUILD_DIR)/subroutines/ei_spec_funct.o

.PHONY: all pp clean
