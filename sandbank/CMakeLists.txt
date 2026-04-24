cmake_minimum_required(VERSION 3.14)
project(RDSReaderProject LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find system ZLIB
find_package(ZLIB REQUIRED)

include(FetchContent)

# We need to ensure byteme (a dependency of rds2cpp) knows to use Zlib
FetchContent_Declare(
  rds2cpp
  GIT_REPOSITORY https://github.com/LTLA/rds2cpp
  GIT_TAG master
)

FetchContent_MakeAvailable(rds2cpp)

add_executable(rds_test readRdata.cpp)

# Link ZLIB and rds2cpp
target_link_libraries(rds_test PRIVATE rds2cpp ZLIB::ZLIB)

# This definition is CRITICAL to enable GzipFileReader in byteme
target_compile_definitions(rds_test PRIVATE BYTEME_USE_ZLIB)