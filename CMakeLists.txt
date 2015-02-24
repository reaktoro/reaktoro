#------------------------------------------------------------------------------
# Version Requirement
#------------------------------------------------------------------------------
cmake_minimum_required(VERSION 2.8)

#------------------------------------------------------------------------------
# Project Details
#------------------------------------------------------------------------------
project(Reaktor)

# Set the cmake module path of the project
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")

# Use cotire to speedup the compilation of the project
#include(cotire)

# Set the version of the project
set(REAKTOR_VERSION_MAJOR "0")
set(REAKTOR_VERSION_MINOR "1")
set(REAKTOR_VERSION_MICRO "0")
set(REAKTOR_VERSION "${REAKTOR_VERSION_MAJOR}.${REAKTOR_VERSION_MINOR}.${REAKTOR_VERSION_MICRO}")

#------------------------------------------------------------------------------
# Build Options
#------------------------------------------------------------------------------
option(BUILD_PYTHON_WRAPPERS  "Build wrappers for the python language." on)
option(BUILD_DEMOS            "Build demos." on)
option(BUILD_DOCUMENTATION    "Build documentation." on)
option(BUILD_TESTS            "Build tests." on)
option(BUILD_STATIC_LIBS      "Build static libraries." on)
option(BUILD_SHARED_LIBS      "Build shared libraries." on)
option(ENABLE_TESTING         "Enable testing." on)
option(BUILD_PHREEQC          "Build PHREEQC." on)
option(BUILD_GEMS3K           "Build GEMS3K." on)

#------------------------------------------------------------------------------
# Build Type
#------------------------------------------------------------------------------
# Set the default build type to RelWithDebInfo
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose the type of build: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif()

#------------------------------------------------------------------------------
# Compiler Flags
#------------------------------------------------------------------------------
# Check for support of the compiler flags
include(CheckCXXCompilerFlag)

check_cxx_compiler_flag(-std=c++11           HAVE_FLAG_STDCXX11)
check_cxx_compiler_flag(-fPIC                HAVE_FLAG_FPIC)
check_cxx_compiler_flag(-Wall                HAVE_FLAG_WALL)
check_cxx_compiler_flag(-pedantic            HAVE_FLAG_PEDANTIC)
check_cxx_compiler_flag(-Wno-variadic-macros HAVE_NO_VARIADIC_MACROS)

# Initialize the list of compiler flags
set(REAKTOR_CXX_FLAGS)

if(HAVE_FLAG_STDCXX11)
  set(REAKTOR_CXX_FLAGS ${REAKTOR_CXX_FLAGS} "-std=c++11")
endif()

if(HAVE_FLAG_FPIC)
  set(REAKTOR_CXX_FLAGS ${REAKTOR_CXX_FLAGS} "-fPIC")
endif()

if(HAVE_FLAG_WALL)
  set(REAKTOR_CXX_FLAGS ${REAKTOR_CXX_FLAGS} "-Wall")
endif()

if(HAVE_FLAG_PEDANTIC)
  set(REAKTOR_CXX_FLAGS ${REAKTOR_CXX_FLAGS} "-pedantic")
endif()

if(HAVE_NO_VARIADIC_MACROS)
  set(REAKTOR_CXX_FLAGS ${REAKTOR_CXX_FLAGS} "-Wno-variadic-macros")
endif()

# Armadillo flags
set(REAKTOR_CXX_FLAGS ${REAKTOR_CXX_FLAGS} "-DARMA_USE_BLAS")
set(REAKTOR_CXX_FLAGS ${REAKTOR_CXX_FLAGS} "-DARMA_USE_LAPACK")
set(REAKTOR_CXX_FLAGS ${REAKTOR_CXX_FLAGS} "-DARMA_DONT_USE_WRAPPER")

# Add the supported compiler flags
add_definitions(${REAKTOR_CXX_FLAGS})

#------------------------------------------------------------------------------
# Project Directories
#------------------------------------------------------------------------------
set(REAKTOR_DEMOS_DIR     ${CMAKE_SOURCE_DIR}/demos)
set(REAKTOR_DOCS_DIR      ${CMAKE_SOURCE_DIR}/docs)
set(REAKTOR_EXTERNAL_DIR  ${CMAKE_SOURCE_DIR}/external)
set(REAKTOR_RESOURCES_DIR ${CMAKE_SOURCE_DIR}/resources)
set(REAKTOR_SOURCE_DIR    ${CMAKE_SOURCE_DIR}/Reaktor)
set(REAKTOR_TESTS_DIR     ${CMAKE_SOURCE_DIR}/tests)
set(REAKTOR_RUNTIME_DIR   ${CMAKE_BINARY_DIR}/bin)
set(REAKTOR_LIBRARY_DIR   ${CMAKE_BINARY_DIR}/lib)
set(REAKTOR_ARQUIVE_DIR   ${CMAKE_BINARY_DIR}/lib)

#------------------------------------------------------------------------------
# Output Directories
#------------------------------------------------------------------------------
# Set the output path of the executable binaries produced by the project
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${REAKTOR_RUNTIME_DIR}"
  CACHE STRING "Single output directory for building all executables.")

# Set the output path of the shared library binaries produced by the project
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${REAKTOR_LIBRARY_DIR}"
  CACHE STRING "Single output directory for building all shared libraries.")

# Set the output path of the static library binaries produced by the project
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${REAKTOR_ARQUIVE_DIR}"
  CACHE STRING "Single output directory for building all static libraries.")

#------------------------------------------------------------------------------
# Output Prefixes and Suffixes
#------------------------------------------------------------------------------
# Set the suffix of the shared libraries
set(SUFFIX_SHARED_LIBS ""
  CACHE STRING "Suffix of the shared libraries.")

# Set the suffix of the static libraries
set(SUFFIX_STATIC_LIBS "-static"
  CACHE STRING "Suffix of the static libraries.")

#------------------------------------------------------------------------------
# Include Directories
#------------------------------------------------------------------------------
include_directories(${CMAKE_SOURCE_DIR})
include_directories(dependencies)
include_directories(utilities)

#------------------------------------------------------------------------------
# Build Project
#------------------------------------------------------------------------------
add_subdirectory(resources)
add_subdirectory(dependencies)
add_subdirectory(Reaktor)

if(BUILD_PYTHON_WRAPPERS)
    add_definitions("-DEIGEN_DEFAULT_TO_ROW_MAJOR")
    add_subdirectory(PyReaktor)
endif()

#------------------------------------------------------------------------------
# Build Demos
#------------------------------------------------------------------------------
if(BUILD_DEMOS)
  # Add the demos subdirectory to the build (not to the default target)
  add_subdirectory(demos EXCLUDE_FROM_ALL)

  # Add target "demos" for building demos
  add_custom_target(demos
  COMMAND ${CMAKE_MAKE_PROGRAM}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/demos")
endif()

#------------------------------------------------------------------------------
# Build Docs
#------------------------------------------------------------------------------
if(BUILD_DOCUMENTATION)
  # Add the docs subdirectory to the build (not to the default target)
  add_subdirectory(docs EXCLUDE_FROM_ALL)
endif()

#------------------------------------------------------------------------------
# Build Tests
#------------------------------------------------------------------------------
if(BUILD_TESTS)
  # Add the tests subdirectory to the build (not to the default target)
  add_subdirectory(tests EXCLUDE_FROM_ALL)

  # Add target "tests" for building tests
  add_custom_target(tests
      COMMAND ${CMAKE_MAKE_PROGRAM}
      WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/tests")
endif()