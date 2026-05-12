cmake_minimum_required(VERSION 3.12)
project(DistanceWrap LANGUAGES C CXX)

# 1. Find Python3 Interpreter AND Development headers
find_package(Python3 REQUIRED COMPONENTS Interpreter Development.Module)

# 2. Find NumPy include directory
execute_process(
    COMMAND "${Python3_EXECUTABLE}" -c "import numpy; print(numpy.get_include())"
    OUTPUT_VARIABLE NUMPY_INCLUDE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE NUMPY_RESULT
)

if(NOT NUMPY_RESULT EQUAL 0)
    message(FATAL_ERROR "NumPy not found. Please install it in the active environment.")
endif()

# 3. Find pybind11
# First try finding it normally
find_package(pybind11 QUIET)
if(NOT pybind11_FOUND)
    # If not found, try to get the path from python
    execute_process(
        COMMAND "${Python3_EXECUTABLE}" -c "import pybind11; print(pybind11.get_cmake_dir())"
        OUTPUT_VARIABLE pybind11_DIR
        OUTPUT_STRIP_TRAILING_WHITESPACE
        RESULT_VARIABLE PYBIND11_RESULT
    )
    if(PYBIND11_RESULT EQUAL 0)
        find_package(pybind11 REQUIRED PATHS ${pybind11_DIR})
    else()
        message(FATAL_ERROR "pybind11 not found. Please install it.")
    endif()
endif()

# 4. Create the Python modules
# MODULE is the proper type for Python .so files (shared objects)
python3_add_library(_distance_wrap MODULE distance_wrap.c)

# Create the pybind11 module
pybind11_add_module(_distance_pybind distance_pybind.cpp)

# 5. Set include directories
target_include_directories(_distance_wrap PRIVATE
    ${NUMPY_INCLUDE_DIR}
    ${Python3_INCLUDE_DIRS}
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_include_directories(_distance_pybind PRIVATE
    ${NUMPY_INCLUDE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# 6. Force Debug mode and Disable Optimization for GDB Walkthrough
if(MSVC)
    set(DEBUG_FLAGS "/Od /Zi /RTC1")
    target_compile_options(_distance_wrap PRIVATE ${DEBUG_FLAGS})
    target_compile_options(_distance_pybind PRIVATE ${DEBUG_FLAGS})
else()
    # -O0 is critical: it prevents the compiler from optimizing out variables
    # -fno-inline ensures we can 'step' into every function call
    set(DEBUG_FLAGS "-O0 -g -fno-inline -fPIC")
    target_compile_options(_distance_wrap PRIVATE -O0 -g -fno-inline -fPIC)
    target_compile_options(_distance_pybind PRIVATE -O0 -g -fno-inline -fPIC)
endif()

# 7. Set output properties
set_target_properties(_distance_wrap PROPERTIES 
    PREFIX ""                       # Remove 'lib' prefix
    OUTPUT_NAME "_distance_wrap"    # Force exact name
)

set_target_properties(_distance_pybind PROPERTIES 
    PREFIX ""
    OUTPUT_NAME "_distance_pybind"
)

# 8. Default to Debug build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build." FORCE)
endif()

# Print status messages for verification
message(STATUS "--- Debugging Config ---")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Python Exec: ${Python3_EXECUTABLE}")
message(STATUS "NumPy include: ${NUMPY_INCLUDE_DIR}")
message(STATUS "pybind11 version: ${pybind11_VERSION}")
