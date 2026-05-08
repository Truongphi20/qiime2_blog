cmake_minimum_required(VERSION 3.12)
project(DistanceWrap LANGUAGES C)

# 1. Find Python3 Interpreter AND Development headers
# We need the Interpreter to run the NumPy check command
find_package(Python3 REQUIRED COMPONENTS Interpreter Development.Module)

# 2. Find NumPy include directory using the found Python3 executable
execute_process(
    COMMAND "${Python3_EXECUTABLE}" -c "import numpy; print(numpy.get_include())"
    OUTPUT_VARIABLE NUMPY_INCLUDE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE NUMPY_RESULT
)

if(NOT NUMPY_RESULT EQUAL 0)
    message(FATAL_ERROR "NumPy not found. Please install it in the active environment.")
endif()

# 3. Create the Python module
# MODULE is the proper type for Python .so files (shared objects)
python3_add_library(_distance_wrap MODULE distance_wrap.c)

# 4. Set include directories
# Include NumPy headers, Python headers, and current directory for distance_impl.h
target_include_directories(_distance_wrap PRIVATE
    ${NUMPY_INCLUDE_DIR}
    ${Python3_INCLUDE_DIRS}
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# 5. Force Debug mode and Disable Optimization for GDB Walkthrough
if(MSVC)
    target_compile_options(_distance_wrap PRIVATE /Od /Zi /RTC1)
else()
    # -O0 is critical: it prevents the compiler from optimizing out variables
    # -fno-inline ensures we can 'step' into every function call
    target_compile_options(_distance_wrap PRIVATE -O0 -g -fno-inline -fPIC)
endif()

# 6. Set output properties
set_target_properties(_distance_wrap PROPERTIES 
    PREFIX ""                       # Remove 'lib' prefix (e.g. _distance_wrap.so, not lib_distance_wrap.so)
    OUTPUT_NAME "_distance_wrap"    # Force exact name
)

# 7. Default to Debug build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build." FORCE)
endif()

# Print status messages for verification
message(STATUS "--- Debugging Config ---")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Python Exec: ${Python3_EXECUTABLE}")
message(STATUS "NumPy include: ${NUMPY_INCLUDE_DIR}")