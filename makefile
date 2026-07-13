include makefile.inc

#IS_IFX := $(findstring ifx, $(FOR))

SRC_DIR   := ./src
BUILD_DIR := ./build
BIN_DIR   := $(DIR)
UTILS_DIR := ./utils

MOD_OUT_FLAG := -J$(BUILD_DIR)

#ifeq ($(IS_IFX), ifx)
#    MOD_OUT_FLAG := -module $(BUILD_DIR)
#else
#    MOD_OUT_FLAG := -J$(BUILD_DIR)
#endif

FC_FLAGS := $(OMP) $(COND) $(EXTRA) -I$(BUILD_DIR)
L_FLAGS  := $(OMP) $(LIBS) $(COND) $(EXTRA) -I$(BUILD_DIR)

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
    gw_subs \
    bcast_input_read 

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
    emission_PL \
    s_gw_mesh 

SUBROUTINE_OBJS := $(foreach name,$(SUBROUTINE_NAMES),$(BUILD_DIR)/subroutines/$(name).o)
SUBPROGRAM_OBJS := $(foreach name,$(SUBPROGRAM_NAMES),$(BUILD_DIR)/subprograms/$(name).o)
MODULE_OBJECTS  := $(SUBROUTINE_OBJS) $(SUBPROGRAM_OBJS)
INPUT_MODULE_OBJ := $(BUILD_DIR)/subroutines/module_input_read.o 
$(filter-out $(INPUT_MODULE_OBJ),$(MODULE_OBJECTS)): $(INPUT_MODULE_OBJ)

MAIN_SRC  := $(SRC_DIR)/wtb_main.F90
MAIN_EXEC := $(BIN_DIR)/wtb.x
LIB_FILE  := $(BUILD_DIR)/libwtb.a

PCE_SOURCES    := $(UTILS_DIR)/slme/pce-code.f90 $(UTILS_DIR)/slme/pce-subs.f90
HUCKEL_SOURCES := $(UTILS_DIR)/huckel2wtb/src/overlaps_jc.f90 $(UTILS_DIR)/huckel2wtb/src/diagonalize.f90 $(UTILS_DIR)/huckel2wtb/src/Huckel_TB.f90
SQ_SOURCES    := $(UTILS_DIR)/slme/sq-curve.f90 $(UTILS_DIR)/slme/pce-subs.f90

all: $(MAIN_EXEC) pp
	@echo "--- Build Complete ---"

main: $(MAIN_EXEC)
	@echo "--- Main Executable Build Complete ---"

$(MAIN_EXEC): $(LIB_FILE) $(MAIN_SRC) makefile.inc
	@mkdir -p $(BIN_DIR)
	@echo "--- Linking Main Executable: $@ ---"
	$(FOR) $(MAIN_SRC) -o $@ $(L_FLAGS) -L$(BUILD_DIR) -lwtb
	@cp $@ ./build/wtb.x

$(LIB_FILE): $(MODULE_OBJECTS)
	@echo "--- Creating Static Library: $@ ---"
	ar rcs $@ $(MODULE_OBJECTS)

$(BUILD_DIR)/subroutines/%.o: $(SRC_DIR)/subroutines/%.F90 makefile.inc
	@mkdir -p $(dir $@)
	@echo "Compiling Subroutine: $< -> $@"
	$(FOR) -c $< -o $@ $(FC_FLAGS) $(MOD_OUT_FLAG)

$(BUILD_DIR)/subprograms/%.o: $(SRC_DIR)/subprograms/%.F90 makefile.inc
	@mkdir -p $(dir $@)
	@echo "Compiling Subprogram: $< -> $@"
	$(FOR) -c $< -o $@ $(FC_FLAGS) $(MOD_OUT_FLAG)

UTILS_EXECS := $(BIN_DIR)/nc_nv_finder.x \
               $(BIN_DIR)/param_gen.x \
               $(BIN_DIR)/param_gen_vasp.x \
               $(BIN_DIR)/absorbance.x \
               $(BIN_DIR)/pce.x \
               $(BIN_DIR)/huckel2wtb.x \
               $(BIN_DIR)/sq_curve.x \
               
pp: $(UTILS_EXECS)
	@echo "--- Copying Python Scripts ---"
	@cp $(UTILS_DIR)/*.py $(BIN_DIR)

$(BIN_DIR)/nc_nv_finder.x: $(UTILS_DIR)/nc_nv_finder.F90 makefile.inc
	$(FOR) $< -o $@ $(L_FLAGS) $(MOD_OUT_FLAG)

$(BIN_DIR)/param_gen.x: $(UTILS_DIR)/param_gen.F90 makefile.inc
	$(FOR) $< -o $@ $(L_FLAGS) $(MOD_OUT_FLAG)

$(BIN_DIR)/param_gen_vasp.x: $(UTILS_DIR)/param_gen_vasp.F90 makefile.inc
	$(FOR) $< -o $@ $(L_FLAGS) $(MOD_OUT_FLAG)

$(BIN_DIR)/absorbance.x: $(UTILS_DIR)/absorbance.F90 makefile.inc
	$(FOR) $< -o $@ $(L_FLAGS) $(MOD_OUT_FLAG)

$(BIN_DIR)/pce.x: $(PCE_SOURCES) makefile.inc
	$(FOR) $(PCE_SOURCES) -o $@ $(L_FLAGS) $(MOD_OUT_FLAG)
	
$(BIN_DIR)/sq_curve.x: $(PCE_SOURCES) makefile.inc
	$(FOR) $(SQ_SOURCES) -o $@ $(L_FLAGS) $(MOD_OUT_FLAG)	

$(BIN_DIR)/huckel2wtb.x: $(HUCKEL_SOURCES) makefile.inc
	$(FOR) $(HUCKEL_SOURCES) -o $@ $(L_FLAGS) $(MOD_OUT_FLAG)

clean:
	@echo "--- Cleaning build, bin, and .mod files ---"
	@rm -rf $(BUILD_DIR) 
	@rm -rf $(BIN_DIR)/*.py
	@rm -rf $(BIN_DIR)/*.x
	@rm -f ./*.mod
	@mkdir -p $(BUILD_DIR) 

#$(BUILD_DIR)/subroutines/coulomb_pot.o: $(BUILD_DIR)/subroutines/ei_spec_funct.o

.PHONY: all pp clean main



