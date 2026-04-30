cmake_minimum_required(VERSION 3.10)
project(mafft VERSION 7.525 LANGUAGES C CXX)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Options
option(ENABLE_MULTITHREAD "Enable multi-threading" ON)
option(ENABLE_ATOMIC "Enable atomic support" OFF)
option(BUILD_EXTENSIONS "Build RNA extensions" ON)
option(ENABLE_DASH_CLIENT "Enable dash_client (requires Go)" OFF)

if(ENABLE_MULTITHREAD)
    add_definitions(-Denablemultithread)
    find_package(Threads REQUIRED)
    set(THREADS_LIB Threads::Threads)
else()
    set(THREADS_LIB "")
endif()

if(ENABLE_ATOMIC)
    add_definitions(-Denableatomic)
    set(CMAKE_C_STANDARD 11)
endif()

# Installation paths
include(GNUInstallDirs)
set(MAFFT_LIBEXECDIR ${CMAKE_INSTALL_LIBEXECDIR}/mafft)

find_library(MATH_LIBRARY m)

# Core
add_subdirectory(core)

# Extensions
if(BUILD_EXTENSIONS)
    add_subdirectory(extensions)
endif()
