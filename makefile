include makefile.inc

FC_FLAGS := $(OMP) $(COND) $(EXTRA)
L_FLAGS  := $(OMP) $(LIBS) $(COND) $(EXTRA)

SRC_DIR   := ./src
BUILD_DIR := ./build
BIN_DIR   := $(DIR)
UTILS_DIR := ./utils

MAIN_SRC := $(SRC_DIR)/wtb_main.F90
MAIN_EXEC := $(BIN_DIR)wtb.x

# Find ALL .F90 files that are modules or subroutines
ALL_MODULE_SOURCES := $(shell find $(SRC_DIR)/subroutines -name '*.F90') \
                      $(shell find $(SRC_DIR)/subprograms -name '*.F90')

EXCLUDE_FILES := $(SRC_DIR)/subroutines/other.F90 \
                 $(SRC_DIR)/subroutines/old/hamiltonian_tb.f90 \
                 $(SRC_DIR)/subprograms/old/tdos-tool.F90 \
                 $(SRC_DIR)/subprograms/old/spin_txt.F90 \
                 $(SRC_DIR)/subprograms/old/pce-code.F90

MODULE_SOURCES := $(filter-out $(EXCLUDE_FILES), $(ALL_MODULE_SOURCES))

MODULE_OBJECTS := $(MODULE_SOURCES:$(SRC_DIR)/%.F90=$(BUILD_DIR)/%.o)

# Default target: build the main program and all utilities
all: $(MAIN_EXEC) pp
	@echo "--- Build Complete ---"

# Rule to link the main executable
$(MAIN_EXEC): $(MODULE_OBJECTS) $(MAIN_SRC)
	@mkdir -p $(BIN_DIR)
	@echo "--- Linking Main Executable: $@ ---"
	$(FOR) $(MAIN_SRC) $(MODULE_OBJECTS) -o $@ $(L_FLAGS)
	@cp $@ ./build/wtb.x

UTILS_EXECS := $(BIN_DIR)nc_nv_finder.x \
               $(BIN_DIR)param_gen.x \
               $(BIN_DIR)param_gen_vasp.x \
               $(BIN_DIR)absorbance.x \
               $(BIN_DIR)pce.x \
               $(BIN_DIR)huckel2wtb.x

pp: $(UTILS_EXECS)
	@echo "--- Copying Python Scripts ---"
	@cp $(UTILS_DIR)/*.py $(BIN_DIR)

# Rules for compiling utils
$(BIN_DIR)nc_nv_finder.x: $(UTILS_DIR)/nc_nv_finder.F90
	$(FOR) $< -o $@ $(L_FLAGS)

$(BIN_DIR)param_gen.x: $(UTILS_DIR)/param_gen.F90
	$(FOR) $< -o $@ $(L_FLAGS)

$(BIN_DIR)param_gen_vasp.x: $(UTILS_DIR)/param_gen_vasp.F90
	$(FOR) $< -o $@ $(L_FLAGS)

$(BIN_DIR)absorbance.x: $(UTILS_DIR)/absorbance.F90
	$(FOR) $< -o $@ $(L_FLAGS)

$(BIN_DIR)pce.x: $(UTILS_DIR)/slme/pce-code.f90 $(UTILS_DIR)/slme/pce-subs.f90
	$(FOR) $^ -o $@ $(L_FLAGS)

$(BIN_DIR)huckel2wtb.x: $(UTILS_DIR)/huckel2wtb/src/overlaps_jc.f90 $(UTILS_DIR)/huckel2wtb/src/diagonalize.f90 $(UTILS_DIR)/huckel2wtb/src/Huckel_TB.f90
	$(FOR) $^ -o $@ $(L_FLAGS)

# This rule teaches 'make' how to compile any .F90 file
# from 'src' into its corresponding .o file in 'build'.
#
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.F90
	@mkdir -p $(dir $@)
	@echo "Compiling Module: $< -> $@"
	$(FOR) -c $< -o $@ $(FC_FLAGS)


$(BUILD_DIR)/subroutines/coulomb_pot.o: $(BUILD_DIR)/subroutines/ei_spec_funct.o

# --- Clean Target ---
clean:
	@echo "--- Cleaning build, bin, and .mod files ---"
	@rm -rf $(BUILD_DIR) $(BIN_DIR) ./*.mod
	@mkdir -p $(BUILD_DIR) $(BIN_DIR)

.PHONY: all pp clean
