cmake_minimum_required(VERSION 3.10)
project(FastTree C)

# Set the C standard
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

# --- Options ---
option(USE_SINGLE "Use single precision floating point (default is double)" OFF)
option(NO_SSE "Turn off use of SSE3 instructions" OFF)
option(USE_NATIVE "Optimize for the local machine architecture (-march=native)" ON)
option(USE_UNSAFE_MATH "Use -funsafe-math-optimizations" ON)
option(TRACK_MEMORY "Enable detailed reports of memory usage (slower)" OFF)

# --- Dependencies ---
# Math library
find_library(MATH_LIBRARY m)

# OpenMP
find_package(OpenMP)

# --- Compiler Flags ---
# Recommended flags from documentation
set(BASE_FLAGS -O0 -Wall)

if(USE_NATIVE)
    list(APPEND BASE_FLAGS -march=native)
endif()

if(USE_UNSAFE_MATH)
    list(APPEND BASE_FLAGS -funsafe-math-optimizations)
endif()

# Check for OpenMP SIMD support
set(SIMD_FLAGS -fopenmp-simd)

# Definitions based on options
if(USE_SINGLE)
    add_definitions(-DUSE_SINGLE)
endif()

if(NO_SSE)
    add_definitions(-DNO_SSE)
endif()

if(TRACK_MEMORY)
    add_definitions(-DTRACK_MEMORY)
endif()

# --- Source Files ---
set(FASTTREE_SOURCES src/FastTree.c src/FastNJ.c src/hyper_parameters.c)

# --- Targets ---

# 1. FastTree (Serial)
add_executable(FastTree ${FASTTREE_SOURCES})
target_compile_options(FastTree PRIVATE ${BASE_FLAGS} ${SIMD_FLAGS})
target_include_directories(FastTree PRIVATE src)
if(MATH_LIBRARY)
    target_link_libraries(FastTree PRIVATE ${MATH_LIBRARY})
endif()

# 2. FastTreeMP (Multi-threaded)
if(OpenMP_C_FOUND)
    add_executable(FastTreeMP ${FASTTREE_SOURCES})
    # FastTreeMP requires -DOPENMP and -fopenmp
    target_compile_options(FastTreeMP PRIVATE ${BASE_FLAGS} ${SIMD_FLAGS} -DOPENMP)
    target_include_directories(FastTreeMP PRIVATE src)
    target_link_libraries(FastTreeMP PRIVATE OpenMP::OpenMP_C)
    if(MATH_LIBRARY)
        target_link_libraries(FastTreeMP PRIVATE ${MATH_LIBRARY})
    endif()
    message(STATUS "FastTreeMP will be built with OpenMP support")
else()
    message(STATUS "OpenMP not found. FastTreeMP will not be built.")
endif()

# --- Installation ---
install(TARGETS FastTree DESTINATION bin)
if(TARGET FastTreeMP)
    install(TARGETS FastTreeMP DESTINATION bin)
endif()
