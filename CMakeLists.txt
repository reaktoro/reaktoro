# Set cmake version requirement
cmake_minimum_required(VERSION 2.8)

# Set the name of the project
project(Reaktoro)

# Set the cmake module path of the project
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")

# Set the version of the project
set(REAKTORO_VERSION_MAJOR "0")
set(REAKTORO_VERSION_MINOR "1")
set(REAKTORO_VERSION_MICRO "0")
set(REAKTORO_VERSION "${REAKTORO_VERSION_MAJOR}.${REAKTORO_VERSION_MINOR}.${REAKTORO_VERSION_MICRO}")

# Define build options
option(BUILD_DEMOS           "Build demos." on)
option(BUILD_DOCS            "Build documentation." on)
option(BUILD_GEMS            "Build GEMS." off)
option(BUILD_PHREEQC         "Build PHREEQC." off)
option(BUILD_PYTHON_WRAPPERS "Build the python wrappers." off)
option(BUILD_SHARED_LIBS     "Build shared libraries." on)
option(BUILD_STATIC_LIBS     "Build static libraries." on)
option(BUILD_TESTS           "Build tests." on)
option(ENABLE_TESTING        "Enable testing." on)

# Set the default build type to RelWithDebInfo
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose the type of build: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif()

# Set the list of compiler flags for GNU compiler
if(${CMAKE_COMPILER_IS_GNUCXX})
  add_definitions("-std=c++11 -fPIC -Wall -pedantic -Wno-variadic-macros")
endif()

# Set the list of compiler flags for MSVC compiler
if(${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
  add_definitions("/W0 -D_SCL_SECURE_NO_WARNINGS /MP4")
endif()

# Set the project directory names
set(REAKTORO_DEMOS_DIR        ${CMAKE_SOURCE_DIR}/demos)
set(REAKTORO_DOCS_DIR         ${CMAKE_SOURCE_DIR}/docs)
set(REAKTORO_EXTERNAL_DIR     ${CMAKE_SOURCE_DIR}/external)
set(REAKTORO_RESOURCES_DIR    ${CMAKE_SOURCE_DIR}/resources)
set(REAKTORO_SOURCE_DIR       ${CMAKE_SOURCE_DIR}/Reaktoro)
set(REAKTORO_TESTS_DIR        ${CMAKE_SOURCE_DIR}/tests)
set(REAKTORO_THIRDPARTY_DIR   ${CMAKE_SOURCE_DIR}/thirdparty)
set(REAKTORO_UTILITIES_DIR    ${CMAKE_SOURCE_DIR}/utilities)

# Set the directories where libraries and binaries are kept
set(REAKTORO_RUNTIME_DIR ${CMAKE_BINARY_DIR}/bin)
set(REAKTORO_LIBRARY_DIR ${CMAKE_BINARY_DIR}/lib)
set(REAKTORO_ARQUIVE_DIR ${CMAKE_BINARY_DIR}/lib)

# Set the Boost include dir
set(BOOST_INCLUDE_DIR ${BOOST_ROOT})

# Set the include and library directories where the third party libraries were installed
set(THIRDPARTY_INCLUDE_DIR ${CMAKE_INSTALL_PREFIX}/include)
set(THIRDPARTY_LIBRARY_DIR ${CMAKE_INSTALL_PREFIX}/lib)

# Set the output path of the executable binaries produced by the project
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${REAKTORO_RUNTIME_DIR}"
  CACHE STRING "Single output directory for building all executables.")

# Set the output path of the shared library binaries produced by the project
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${REAKTORO_LIBRARY_DIR}"
  CACHE STRING "Single output directory for building all shared libraries.")

# Set the output path of the static library binaries produced by the project
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${REAKTORO_ARQUIVE_DIR}"
  CACHE STRING "Single output directory for building all static libraries.")

# Set the suffix of the shared libraries
set(SUFFIX_SHARED_LIBS ""
  CACHE STRING "Suffix of the shared libraries.")

# Set the suffix of the static libraries
set(SUFFIX_STATIC_LIBS "-static"
  CACHE STRING "Suffix of the static libraries.")

# Set the include directories
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${THIRDPARTY_INCLUDE_DIR})
include_directories(${BOOST_INCLUDE_DIR})

# Set the link directories
link_directories(${THIRDPARTY_LIBRARY_DIR})

# Build Project
add_subdirectory(${REAKTORO_THIRDPARTY_DIR})
add_subdirectory(${REAKTORO_SOURCE_DIR})

if(BUILD_PYTHON_WRAPPERS)
  add_subdirectory(python)
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
if(BUILD_DOCS)
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
