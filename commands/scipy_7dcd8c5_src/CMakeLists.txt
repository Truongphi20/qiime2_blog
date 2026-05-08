cmake_minimum_required(VERSION 3.12)
project(DistanceWrap LANGUAGES C)

# Find Python3 and its development components
find_package(Python3 REQUIRED COMPONENTS Development)

# Find NumPy include directory
execute_process(
    COMMAND "${Python3_EXECUTABLE}" -c "import numpy; print(numpy.get_include())"
    OUTPUT_VARIABLE NUMPY_INCLUDE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE NUMPY_RESULT
)

if(NOT NUMPY_RESULT EQUAL 0)
    message(FATAL_ERROR "NumPy not found. Please install it (e.g., 'pip install numpy').")
endif()

# Create the Python module
# Using MODULE instead of SHARED is standard for Python extensions
python3_add_library(_distance_wrap MODULE distance_wrap.c)

# Set include directories
target_include_directories(_distance_wrap PRIVATE
    ${NUMPY_INCLUDE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# Force Debug mode and Disable Optimization
# This satisfies the requirement: "with debug mode, don't use optimizer"
if(MSVC)
    target_compile_options(_distance_wrap PRIVATE /Od /Zi /RTC1)
else()
    target_compile_options(_distance_wrap PRIVATE -O0 -g -fno-inline)
endif()

# Ensure the output name is exactly _distance_wrap.so (remove 'lib' prefix)
set_target_properties(_distance_wrap PROPERTIES PREFIX "")

# Default to Debug build if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build." FORCE)
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Python include: ${Python3_INCLUDE_DIRS}")
message(STATUS "NumPy include: ${NUMPY_INCLUDE_DIR}")
