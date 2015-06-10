#------------------------------------------------------------------------------
# Version Requirement
#------------------------------------------------------------------------------
cmake_minimum_required(VERSION 2.8)

#------------------------------------------------------------------------------
# Project Details
#------------------------------------------------------------------------------
project(Reaktoro)

# Set the cmake module path of the project
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")

# Set the version of the project
set(REAKTORO_VERSION_MAJOR "0")
set(REAKTORO_VERSION_MINOR "1")
set(REAKTORO_VERSION_MICRO "0")
set(REAKTORO_VERSION "${REAKTORO_VERSION_MAJOR}.${REAKTORO_VERSION_MINOR}.${REAKTORO_VERSION_MICRO}")

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
option(BUILD_GEMS             "Build GEMS." on)

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
set(REAKTORO_CXX_FLAGS)

if(HAVE_FLAG_STDCXX11)
  set(REAKTORO_CXX_FLAGS ${REAKTORO_CXX_FLAGS} "-std=c++11")
endif()

if(HAVE_FLAG_FPIC)
  set(REAKTORO_CXX_FLAGS ${REAKTORO_CXX_FLAGS} "-fPIC")
endif()

if(HAVE_FLAG_WALL)
  set(REAKTORO_CXX_FLAGS ${REAKTORO_CXX_FLAGS} "-Wall")
endif()

if(HAVE_FLAG_PEDANTIC)
  set(REAKTORO_CXX_FLAGS ${REAKTORO_CXX_FLAGS} "-pedantic")
endif()

if(HAVE_NO_VARIADIC_MACROS)
  set(REAKTORO_CXX_FLAGS ${REAKTORO_CXX_FLAGS} "-Wno-variadic-macros")
endif()

# Add the supported compiler flags
add_definitions(${REAKTORO_CXX_FLAGS})

#------------------------------------------------------------------------------
# Project Directories
#------------------------------------------------------------------------------
set(REAKTORO_DEMOS_DIR     ${CMAKE_SOURCE_DIR}/demos)
set(REAKTORO_DOCS_DIR      ${CMAKE_SOURCE_DIR}/docs)
set(REAKTORO_EXTERNAL_DIR  ${CMAKE_SOURCE_DIR}/external)
set(REAKTORO_RESOURCES_DIR ${CMAKE_SOURCE_DIR}/resources)
set(REAKTORO_SOURCE_DIR    ${CMAKE_SOURCE_DIR}/Reaktoro)
set(REAKTORO_TESTS_DIR     ${CMAKE_SOURCE_DIR}/tests)
set(REAKTORO_RUNTIME_DIR   ${CMAKE_BINARY_DIR}/bin)
set(REAKTORO_LIBRARY_DIR   ${CMAKE_BINARY_DIR}/lib)
set(REAKTORO_ARQUIVE_DIR   ${CMAKE_BINARY_DIR}/lib)

#------------------------------------------------------------------------------
# Output Directories
#------------------------------------------------------------------------------
# Set the output path of the executable binaries produced by the project
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${REAKTORO_RUNTIME_DIR}"
  CACHE STRING "Single output directory for building all executables.")

# Set the output path of the shared library binaries produced by the project
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${REAKTORO_LIBRARY_DIR}"
  CACHE STRING "Single output directory for building all shared libraries.")

# Set the output path of the static library binaries produced by the project
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${REAKTORO_ARQUIVE_DIR}"
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
include_directories(dependencies/cvode/include)
include_directories(utilities)

#------------------------------------------------------------------------------
# Build Project
#------------------------------------------------------------------------------
add_subdirectory(resources)
add_subdirectory(dependencies)
add_subdirectory(Reaktoro)

if(BUILD_PYTHON_WRAPPERS)
    add_definitions("-DEIGEN_DEFAULT_TO_ROW_MAJOR")
    add_subdirectory(PyReaktoro)
    add_subdirectory(reaktoro)
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