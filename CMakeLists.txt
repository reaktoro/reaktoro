# Set cmake version requirement
cmake_minimum_required(VERSION 2.8)

# Set the name of the project
project(Reaktoro)

# Set the cmake module path of the project
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")

# Include the cmake variables with values for installation directories
include(GNUInstallDirs)

# Set the version of the project
set(REAKTORO_VERSION_MAJOR "0")
set(REAKTORO_VERSION_MINOR "1")
set(REAKTORO_VERSION_MICRO "0")
set(REAKTORO_VERSION "${REAKTORO_VERSION_MAJOR}.${REAKTORO_VERSION_MINOR}.${REAKTORO_VERSION_MICRO}")

# Define which components of Reaktoro to build
option(BUILD_ALL             "Build everything." off)
option(BUILD_DEMOS           "Build demos." off)
option(BUILD_DOCS            "Build documentation." off)
option(BUILD_GEMS            "Build GEMS." off)
option(BUILD_PHREEQC         "Build PHREEQC." off)
option(BUILD_PYTHON_WRAPPERS "Build the python wrappers." off)
option(BUILD_TESTS           "Build tests." off)

# Modify the BUILD_XXX variables accordingly to BUILD_ALL
set(BUILD_DEMOS           BUILD_ALL OR BUILD_DEMOS)
set(BUILD_DOCS            BUILD_ALL OR BUILD_DOCS)
set(BUILD_GEMS            BUILD_ALL OR BUILD_GEMS)
set(BUILD_PHREEQC         BUILD_ALL OR BUILD_PHREEQC)
set(BUILD_PYTHON_WRAPPERS BUILD_ALL OR BUILD_PYTHON_WRAPPERS)
set(BUILD_TESTS           BUILD_ALL OR BUILD_TESTS)

# Define which types of libraries to build
option(BUILD_SHARED_LIBS "Build shared libraries." on)
option(BUILD_STATIC_LIBS "Build static libraries." on)

# Define custom options
option(ENABLE_TESTING "Enable testing." on)

# Set the default build type to Release
if(NOT CMAKE_BUILD_TYPE)
    message(STATUS "Setting build type to 'Release' as none was specified.")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
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

# Set the Boost include dir
set(BOOST_INCLUDE_DIR ${BOOST_ROOT})

# Set the directory where the third-party libraries are installed after build
set(THIRDPARTY_DIR ${CMAKE_BINARY_DIR}/thirdparty/install)

# Set the include and library directories where the third party libraries were installed
set(THIRDPARTY_INCLUDE_DIR ${THIRDPARTY_DIR}/${CMAKE_INSTALL_INCLUDEDIR})
set(THIRDPARTY_LIBRARY_DIR ${THIRDPARTY_DIR}/${CMAKE_INSTALL_LIBDIR} ${THIRDPARTY_DIR}/lib)

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

# Configure the Reaktoro package target
set(CPACK_PACKAGE_NAME "Reaktoro")
set(CPACK_PACKAGE_VENDOR "Allan Leal (allan.leal@ethz.ch)")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "A unified framework for modeling chemically reactive systems")
set(CPACK_PACKAGE_VERSION ${REAKTORO_VERSION})
set(CPACK_PACKAGE_VERSION_MAJOR ${REAKTORO_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${REAKTORO_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${REAKTORO_VERSION_MICRO})
set(CPACK_PACKAGE_INSTALL_DIRECTORY "Reaktoro")

# Configure the package components (CPACK_COMPONENT_${COMPONENT_NAME_ALL_CAPS}_GROUP).
set(CPACK_COMPONENTS_ALL applications development interfaces)

# More descriptive names for each of the components, and component groups
set(CPACK_COMPONENT_APPLICATIONS_DISPLAY_NAME "Applications")
set(CPACK_COMPONENT_APPLICATIONS_DESCRIPTION "The applications developed using Reaktoro")

set(CPACK_COMPONENT_DEVELOPMENT_DISPLAY_NAME "Development")
set(CPACK_COMPONENT_DEVELOPMENT_DESCRIPTION "The Reaktoro C++ headers and libraries for application development")

set(CPACK_COMPONENT_INTERFACES_DISPLAY_NAME "Interfaces")
set(CPACK_COMPONENT_INTERFACES_DESCRIPTION "The interfaces for accessing Reaktoro from programming languages other than C++")

# Set the license file for output during the installation
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/COPYING")

set(CPACK_NSIS_ENABLE_UNINSTALL_BEFORE_INSTALL ON)
set(CPACK_NSIS_MODIFY_PATH ON)

# Set the icon used inside the installer wizard and the icons for the installer and un-installer
# set(CPACK_PACKAGE_ICON "${CMAKE_SOURCE_DIR}/hw.bmp")
# set(CPACK_NSIS_MUI_ICON "${CMAKE_SOURCE_DIR}/hw.ico")
# set(CPACK_NSIS_MUI_UNIICON "${CMAKE_SOURCE_DIR}/hw.ico")

# Set the programs displayed under the Start menu shortcut
# We do not use CPACK_PACKAGE_EXECUTABLES because it expects the binaries to be in a
# directory named "bin" (CMake bug tracker ID: 0007829).
set(CPACK_NSIS_MENU_LINKS "ireaktoro" "iReaktoro")

include(CPack)
