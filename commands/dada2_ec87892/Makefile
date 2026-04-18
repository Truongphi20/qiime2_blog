# --- Environment Paths ---
CONDA_PREFIX = /opt/conda/envs/qiime2-amplicon-2026.1
R_LIB_PATH   = $(CONDA_PREFIX)/lib/R/library

R_INC      = $(CONDA_PREFIX)/lib/R/include
RCPP_INC   = $(R_LIB_PATH)/Rcpp/include
RP_INC     = $(R_LIB_PATH)/RcppParallel/include

# --- Compiler Settings ---
CXX = g++
# -DX64 handles the architecture check in your source
# -Isrc allows #include "dada.h" to work from test_main.cpp
CXXFLAGS = -g -O2 -std=c++11 -Wall -fPIC -DX64 -Isrc
INCLUDES = -I"$(R_INC)" -I"$(RCPP_INC)" -I"$(RP_INC)"

# --- Files and Directories ---
SRC_DIR = src
OBJ_DIR = obj
TARGET  = dada2_debug

# Find all .cpp in src/ but ignore R-specific entry points
SRCS = $(shell find $(SRC_DIR) -name "*.cpp" ! -name "Rmain.cpp" ! -name "RcppExports.cpp")
MAIN = test_main.cpp
OBJS = $(SRCS:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)

# --- Linker Flags ---
LIBS = -L$(CONDA_PREFIX)/lib/R/lib -lR \
       -L$(R_LIB_PATH)/RcppParallel/lib -lRcppParallel

# --- Build Rules ---

all: $(TARGET)

# Link the final debug binary
$(TARGET): $(OBJS) $(MAIN)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(MAIN) $(OBJS) -o $(TARGET) $(LIBS)

# Compile source files into the obj directory
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -rf $(OBJ_DIR) $(TARGET)

.PHONY: all clean