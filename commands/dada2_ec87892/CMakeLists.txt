cmake_minimum_required(VERSION 3.14)
project(DADA2_Standalone LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE Debug)
add_compile_options(-O0 -g -Wall -Wno-unknown-pragmas -fPIC)
add_compile_definitions(X64)

set(CONDA_PREFIX "/opt/conda/envs/qiime2-amplicon-2026.1")
set(R_LIB_PATH "${CONDA_PREFIX}/lib/R/library")
set(R_SHARED_LIB_PATH "${CONDA_PREFIX}/lib/R/lib")

include_directories(
    "src"
    "${CONDA_PREFIX}/lib/R/include"
    "${R_LIB_PATH}/Rcpp/include"
    "${R_LIB_PATH}/RcppParallel/include"
)

file(GLOB_RECURSE DADA_SOURCES "src/*.cpp")
list(FILTER DADA_SOURCES EXCLUDE REGEX "RcppExports.cpp")

## denoise executor
add_executable(denoise_dada_debug denoise_main.cpp ${DADA_SOURCES})

target_link_libraries(denoise_dada_debug 
    PRIVATE 
    pthread
    dl
    "${R_SHARED_LIB_PATH}/libR.so"
    "${R_LIB_PATH}/RcppParallel/lib/libtbb.so"
    "${R_LIB_PATH}/RcppParallel/lib/libtbbmalloc.so"
)

set_target_properties(denoise_dada_debug PROPERTIES 
    BUILD_WITH_INSTALL_RPATH TRUE
    INSTALL_RPATH "${R_SHARED_LIB_PATH}:${R_LIB_PATH}/RcppParallel/lib:${CONDA_PREFIX}/lib"
)

## errlearn executor
add_executable(errlearn_dada_debug errlearn_main.cpp ${DADA_SOURCES})

target_link_libraries(errlearn_dada_debug 
    PRIVATE 
    pthread
    dl
    "${R_SHARED_LIB_PATH}/libR.so"
    "${R_LIB_PATH}/RcppParallel/lib/libtbb.so"
    "${R_LIB_PATH}/RcppParallel/lib/libtbbmalloc.so"
)

set_target_properties(errlearn_dada_debug PROPERTIES 
    BUILD_WITH_INSTALL_RPATH TRUE
    INSTALL_RPATH "${R_SHARED_LIB_PATH}:${R_LIB_PATH}/RcppParallel/lib:${CONDA_PREFIX}/lib"
)

## bimera executor
add_executable(bimera_dada_debug bimera_main.cpp ${DADA_SOURCES})

target_link_libraries(bimera_dada_debug 
    PRIVATE 
    pthread
    dl
    "${R_SHARED_LIB_PATH}/libR.so"
    "${R_LIB_PATH}/RcppParallel/lib/libtbb.so"
    "${R_LIB_PATH}/RcppParallel/lib/libtbbmalloc.so"
)

set_target_properties(bimera_dada_debug PROPERTIES 
    BUILD_WITH_INSTALL_RPATH TRUE
    INSTALL_RPATH "${R_SHARED_LIB_PATH}:${R_LIB_PATH}/RcppParallel/lib:${CONDA_PREFIX}/lib"
)