# --- Environment Paths ---
CONDA_PREFIX = /opt/conda/envs/qiime2-amplicon-2026.1
R_LIB_PATH   = $(CONDA_PREFIX)/lib/R/library

# Header locations
R_INC      = $(CONDA_PREFIX)/lib/R/include
RCPP_INC   = $(R_LIB_PATH)/Rcpp/include
RP_INC     = $(R_LIB_PATH)/RcppParallel/include

# Library locations
TBB_LIB_PATH = $(R_LIB_PATH)/RcppParallel/lib
R_SHARED_LIB_PATH = $(CONDA_PREFIX)/lib/R/lib
RCPP_SO = $(R_LIB_PATH)/Rcpp/libs/Rcpp.so

# --- Compiler Settings ---
CXX = g++
# -Wno-cast-user-defined: Silences TBB warnings
# -DX64: Required for DADA2 architecture checks
CXXFLAGS = -g -O2 -std=c++11 -Wall -fPIC -DX64 -Isrc -Wno-cast-user-defined
INCLUDES = -I"$(R_INC)" -I"$(RCPP_INC)" -I"$(RP_INC)"

# --- Files and Directories ---
SRC_DIR = src
OBJ_DIR = obj
TARGET  = dada2_debug

# Exclude RcppExports.cpp to prevent duplicate symbol definitions
SRCS = $(shell find $(SRC_DIR) -name "*.cpp" ! -name "RcppExports.cpp")
MAIN = test_main.cpp
OBJS = $(SRCS:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)

# --- Linker Flags ---
# -ldl: Required for dlopen() in test_main.cpp
LIBS = -Wl,--no-as-needed \
       $(RCPP_SO) \
       -L$(R_SHARED_LIB_PATH) \
       -L$(TBB_LIB_PATH) \
       -lR -ltbb -lpthread -ldl \
       -Wl,-rpath,$(R_SHARED_LIB_PATH) \
       -Wl,-rpath,$(TBB_LIB_PATH) \
       -Wl,-rpath,$(dir $(RCPP_SO))

# --- Build Rules ---

all: $(TARGET)

$(TARGET): $(OBJS) $(MAIN)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(MAIN) $(OBJS) -o $(TARGET) $(LIBS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -rf $(OBJ_DIR) $(TARGET)

.PHONY: all clean