# --- Environment Paths ---
CONDA_PREFIX = /opt/conda/envs/qiime2-amplicon-2026.1
R_LIB_PATH   = $(CONDA_PREFIX)/lib/R/library

# Header locations
R_INC      = $(CONDA_PREFIX)/lib/R/include
RCPP_INC   = $(R_LIB_PATH)/Rcpp/include
RP_INC     = $(R_LIB_PATH)/RcppParallel/include

# Library locations (where libtbb.so and libR.so live)
TBB_LIB_PATH = $(R_LIB_PATH)/RcppParallel/lib
R_SHARED_LIB_PATH = $(CONDA_PREFIX)/lib/R/lib

# --- Compiler Settings ---
CXX = g++
# -DX64 handles architecture checks in the dada2 source
# -Isrc allows headers to be found relative to the src folder
CXXFLAGS = -g -O2 -std=c++11 -Wall -fPIC -DX64 -Isrc
INCLUDES = -I"$(R_INC)" -I"$(RCPP_INC)" -I"$(RP_INC)"

# --- Files and Directories ---
SRC_DIR = src
OBJ_DIR = obj
TARGET  = dada2_debug

# Find all .cpp files in src/ but ignore RcppExports (usually not needed for standalone)
# Includes Rmain.cpp which contains the dada_uniques logic
SRCS = $(shell find $(SRC_DIR) -name "*.cpp" ! -name "RcppExports.cpp")
MAIN = test_main.cpp
OBJS = $(SRCS:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)

# --- Linker Flags ---
# -L: Where to look during build
# -Wl,-rpath: Where the binary should look when it actually runs
LIBS = -L$(R_SHARED_LIB_PATH) \
       -L$(TBB_LIB_PATH) \
       -lR -ltbb -lpthread \
       -Wl,-rpath,$(R_SHARED_LIB_PATH) \
       -Wl,-rpath,$(TBB_LIB_PATH)

# --- Build Rules ---

all: $(TARGET)

# Link the final debug binary
$(TARGET): $(OBJS) $(MAIN)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(MAIN) $(OBJS) -o $(TARGET) $(LIBS)

# Compile source files into the obj directory
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Clean up build artifacts
clean:
	rm -rf $(OBJ_DIR) $(TARGET)

.PHONY: all clean