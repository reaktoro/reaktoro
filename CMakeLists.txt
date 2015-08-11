# Set cmake version requirement
cmake_minimum_required(VERSION 2.8)

# Set the name of the project
project(Reaktoro)

# Include the cmake variables with values for installation directories
include(GNUInstallDirs)

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

# Set the default build type to Release
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif()

# Set the list of compiler flags for GNU compiler
if(${CMAKE_COMPILER_IS_GNUCXX})
  add_definitions("-std=c++11 -Wall -pedantic -Wno-variadic-macros -Wno-deprecated")

  if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    add_definitions("-fPIC")
  endif()
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

# Set the directories where header files, libraries, binaries, etc. are installed
set(REAKTORO_INSTALL_INCLUDEDIR ${CMAKE_INSTALL_FULL_INCLUDEDIR})
set(REAKTORO_INSTALL_LIBDIR     ${CMAKE_INSTALL_FULL_LIBDIR})

# Set the Boost include dir
set(BOOST_INCLUDE_DIR ${BOOST_ROOT})

# Set the include and library directories where the third party libraries were installed
set(THIRDPARTY_INCLUDE_DIR ${CMAKE_INSTALL_FULL_INCLUDEDIR} ${CMAKE_INSTALL_PREFIX}/include)
set(THIRDPARTY_LIBRARY_DIR ${CMAKE_INSTALL_FULL_LIBDIR} ${CMAKE_INSTALL_PREFIX}/lib)

# Set the include directories
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${THIRDPARTY_INCLUDE_DIR})
include_directories(${BOOST_INCLUDE_DIR})

# Set the link directories
link_directories(${THIRDPARTY_LIBRARY_DIR})

# Build the third-party libraries
add_subdirectory(${REAKTORO_THIRDPARTY_DIR})

# Build the C++ lirary Reaktoro
add_subdirectory(${REAKTORO_SOURCE_DIR})

# Build the Python extension module PyReaktoro and the Python package reaktoro
if(BUILD_PYTHON_WRAPPERS)
  add_subdirectory(python)
endif()

# Build the demos
if(BUILD_DEMOS)
  # Add the demos subdirectory to the build (not to the default target)
  add_subdirectory(demos EXCLUDE_FROM_ALL)

  # Add target "demos" for building demos
  add_custom_target(demos
    COMMAND ${CMAKE_MAKE_PROGRAM}
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/demos")
endif()

# Build the documentation
if(BUILD_DOCS)
  # Add the docs subdirectory to the build (not to the default target)
  add_subdirectory(docs EXCLUDE_FROM_ALL)
endif()

# Build the tests
if(BUILD_TESTS)
  # Add the tests subdirectory to the build (not to the default target)
  add_subdirectory(tests EXCLUDE_FROM_ALL)

  # Add target "tests" for building tests
  add_custom_target(tests
    COMMAND ${CMAKE_MAKE_PROGRAM}
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/tests")
endif()
