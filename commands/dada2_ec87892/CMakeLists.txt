cmake_minimum_required(VERSION 3.14)
project(DADA2_Standalone LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE Debug)

# -fPIC is often required when linking against R shared libraries
add_compile_options(-O0 -g -Wall -Wno-unknown-pragmas -fPIC)
add_compile_definitions(X64)
cmake_policy(SET CMP0135 NEW)

# --- Path Configuration ---
set(CONDA_PREFIX "/opt/conda/envs/qiime2-amplicon-2026.1")
set(R_LIB_PATH "${CONDA_PREFIX}/lib/R/library")
set(R_SHARED_LIB_PATH "${CONDA_PREFIX}/lib/R/lib")

# --- Dependencies ---
include(FetchContent)
FetchContent_Declare(zlib URL https://github.com/madler/zlib/archive/refs/tags/v1.3.1.tar.gz DOWNLOAD_EXTRACT_TIMESTAMP TRUE)
FetchContent_Declare(rds2cpp GIT_REPOSITORY https://github.com/LTLA/rds2cpp GIT_TAG master)
FetchContent_MakeAvailable(zlib rds2cpp)

# --- Include Directories ---
include_directories(
    "src"
    "${CONDA_PREFIX}/lib/R/include"
    "${R_LIB_PATH}/Rcpp/include"
    "${R_LIB_PATH}/RcppParallel/include"
)

# --- Source Management ---
file(GLOB_RECURSE DADA_SOURCES "src/*.cpp")

# FIX: Remove RcppExports.cpp but KEEP Rmain.cpp (where dada_uniques is defined)
list(FILTER DADA_SOURCES EXCLUDE REGEX "RcppExports.cpp")

# --- Executable ---
add_executable(dada_debug test_main.cpp ${DADA_SOURCES})

# --- Linking ---
# We need to link libR and the TBB library used by RcppParallel
target_link_libraries(dada_debug 
    PRIVATE 
    rds2cpp 
    zlibstatic 
    pthread
    dl
    "${R_SHARED_LIB_PATH}/libR.so"
    "${R_LIB_PATH}/RcppParallel/lib/libtbb.so"
    "${R_LIB_PATH}/RcppParallel/lib/libtbbmalloc.so"
)

# --- RPATH Configuration ---
# This ensures the executable finds libR.so and libtbb.so at runtime without LD_LIBRARY_PATH
set_target_properties(dada_debug PROPERTIES 
    BUILD_WITH_INSTALL_RPATH TRUE
    INSTALL_RPATH "${R_SHARED_LIB_PATH}:${R_LIB_PATH}/RcppParallel/lib"
)