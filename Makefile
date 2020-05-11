CXX = g++

# path #
SRC_PATH = src
BUILD_PATH = build
BIN_PATH = $(BUILD_PATH)/bin

# executable # 
BIN_NAME = phoebe

# extensions #
SRC_EXT = cpp

# code lists #
# Find all source files in the source directory, sorted by
# most recently modified
#SOURCES = $(shell find $(SRC_PATH) -name '*.$(SRC_EXT)' | sort -k 1nr | cut -f2-)
SOURCES = src/constants/periodic_table.cpp src/algebra/utilities.cpp src/exceptions/exceptions.cpp src/io.cpp src/main.cpp src/statistics.cpp src/context.cpp src/crystal.cpp src/apps/app.cpp src/apps/phonon_transport_app.cpp src/harmonic/harmonic.cpp src/harmonic/phonon_h0.cpp src/points.cpp src/state.cpp src/bandstructure.cpp src/harmonic/electron_h0_fourier.cpp src/harmonic/electron_h0_wannier.cpp src/harmonic/window.cpp src/bte/vector_bte.cpp src/pugixml.cpp src/parser/qe_input_parser.cpp src/observable/observable.cpp src/bte/scattering.cpp src/apps/dos_app.cpp src/bte/drift.cpp src/delta_function/delta_function.cpp src/apps/bands_app.cpp

# Set the object file names, with the source directory stripped
# from the path, and the build path prepended in its place
OBJECTS = $(SOURCES:$(SRC_PATH)/%.$(SRC_EXT)=$(BUILD_PATH)/%.o)
# Set the dependency files that will be used to add header dependencies
DEPS = $(OBJECTS:.o=.d)

# flags #
COMPILE_FLAGS = -std=c++17 -Wall -Wextra -O3 -L./lib/spglib -g
INCLUDES = -I include -I /usr/local/include -I include/Eigen -I lib/pugixml-1.10/src
# Space-separated pkg-config libraries used by this project
LIBS = lib/libsymspg.a

.PHONY: default_target
default_target: release

.PHONY: release
release: export CXXFLAGS := $(CXXFLAGS) $(COMPILE_FLAGS)
release: dirs libs
	@$(MAKE) all

.PHONY: dirs
dirs:
	@echo "Creating directories"
	@mkdir -p $(dir $(OBJECTS))
	@mkdir -p $(BIN_PATH)

.PHONY: libs
libs:
	@echo "Creating libraries"
	(cd lib && make all)

.PHONY: doc
doc:
	(cd doc && doxygen)

.PHONY: clean
clean:
	@echo "Deleting $(BIN_NAME) symlink"
	@$(RM) $(BIN_NAME)
	@echo "Deleting directories"
	@$(RM) -r $(BUILD_PATH)
	@$(RM) -r $(BIN_PATH)
	(cd lib && make clean)
	@$(RM) -r doc/html doc/latex

# (cd ./lib && make clean)

# checks the executable and symlinks to the output
.PHONY: all
all: $(BIN_PATH)/$(BIN_NAME)
	@echo "Making symlink: $(BIN_NAME) -> $<"
	@$(RM) $(BIN_NAME)
	@ln -s $(BIN_PATH)/$(BIN_NAME) $(BIN_NAME)

# Creation of the executable
$(BIN_PATH)/$(BIN_NAME): $(OBJECTS)
	@echo "Linking: $@"
	$(CXX) $(OBJECTS) -o $@ $(LIBS)

# Add dependency files, if they exist
-include $(DEPS)

# Source file rules
# After the first compilation they will be joined with the rules from the
# dependency files to provide header dependencies
$(BUILD_PATH)/%.o: $(SRC_PATH)/%.$(SRC_EXT)
	@echo "Compiling: $< -> $@"
	$(CXX) $(CXXFLAGS) $(INCLUDES) -MP -MMD -c $< -o $@
