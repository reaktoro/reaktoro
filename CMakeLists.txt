# Set cmake version requirement
cmake_minimum_required(VERSION 2.8)

# Set the name of the project
project(Reaktoro)

# Set the cmake module path of the project
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Include the cmake variables with values for installation directories
include(GNUInstallDirs)

# Set the version of the project
set(REAKTORO_VERSION_MAJOR "1")
set(REAKTORO_VERSION_MINOR "0")
set(REAKTORO_VERSION_MICRO "0")
set(REAKTORO_VERSION "${REAKTORO_VERSION_MAJOR}.${REAKTORO_VERSION_MINOR}.${REAKTORO_VERSION_MICRO}")

# Define which components of Reaktoro to build
option(BUILD_ALL         "Build everything." OFF)
option(BUILD_DEMOS       "Build demos." OFF)
option(BUILD_DOCS        "Build documentation." OFF)
option(BUILD_INTERPRETER "Build the reacktoro interpreter package." OFF)
option(BUILD_PYTHON      "Build the python wrappers and python package reaktoro." OFF)
option(BUILD_TESTS       "Build tests." OFF)

# Options for linking Reaktoro with other codes
option(LINK_GEMS    "Link Reaktoro with GEMS." ON)
option(LINK_PHREEQC "Link Reaktoro with PHREEQC." ON)

# Check if GEMS is to be linked with Reaktoro
if(LINK_GEMS)
    add_definitions(-DLINK_GEMS)
endif()

# Check if PHREEQC is to be linked with Reaktoro
if(LINK_PHREEQC)
    add_definitions(-DLINK_PHREEQC)
endif()

# Modify the BUILD_XXX variables accordingly to BUILD_ALL
if(BUILD_ALL MATCHES ON)
    set(BUILD_DEMOS       ON)
    set(BUILD_DOCS        ON)
    set(BUILD_INTERPRETER ON)
    set(BUILD_PYTHON      ON)
    set(BUILD_TESTS       ON)
endif()

# Find Boost
find_package(Boost REQUIRED)

# Find the Python interpreter
find_package(PythonInterp)

# Find PyInstaller
find_package(PyInstaller)

# Check if python was found
if(NOT PYTHONINTERP_FOUND)
    set(BUILD_PYTHON OFF)
    message(WARNING "Python executable was not found. "
        "Cannot build python wrappers and packages for Reaktoro.")
endif()

# Check if PyInstaller was found
if(NOT PYINSTALLER_FOUND)
    message(WARNING "PyInstaller was not found. "
        "Cannot generate an executable application for Reaktoro's interpreter.")
endif()

# Ensure Eigen uses row-major matrices for compatibility with the python wrappers
if(BUILD_PYTHON)
    add_definitions("-DEIGEN_DEFAULT_TO_ROW_MAJOR")
endif()

# Set the output directories of the built libraries and binaries
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Define which types of libraries to build
option(BUILD_SHARED_LIBS "Build shared libraries." ON)
option(BUILD_STATIC_LIBS "Build static libraries." ON)

# Reaktoro currently is not setup to produce a dynamic library using MSVC, only static
if(${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
    set(BUILD_SHARED_LIBS OFF)
endif()

# Define custom options
option(ENABLE_TESTING "Enable testing." on)

# Set the default build type to Release
if(NOT CMAKE_BUILD_TYPE)
    message(STATUS "Setting build type to 'Release' as none was specified.")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Set the list of compiler flags for GNU compiler
if(${CMAKE_CXX_COMPILER_ID} MATCHES "GNU" OR ${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
    add_definitions("-std=c++11 -Wall -pedantic -Wno-variadic-macros -Wno-deprecated")

    if(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Windows")
        add_definitions("-fPIC")
    endif()
endif()

# Set the list of compiler flags for MSVC compiler
if(${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
    add_definitions("/W0 -D_SCL_SECURE_NO_WARNINGS /MP4")
endif()

# Set the project directory names
set(REAKTORO_DATABASES_DIR  ${CMAKE_SOURCE_DIR}/databases)
set(REAKTORO_DEMOS_DIR      ${CMAKE_SOURCE_DIR}/demos)
set(REAKTORO_DOCS_DIR       ${CMAKE_SOURCE_DIR}/docs)
set(REAKTORO_EXTERNAL_DIR   ${CMAKE_SOURCE_DIR}/external)
set(REAKTORO_RESOURCES_DIR  ${CMAKE_SOURCE_DIR}/resources)
set(REAKTORO_SOURCE_DIR     ${CMAKE_SOURCE_DIR}/Reaktoro)

# Set the local directory where the third-party libraries are installed
set(THIRDPARTY_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/thirdparty/install)

# Set the include directories
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${BOOST_INCLUDE_DIR})
include_directories(${THIRDPARTY_INSTALL_PREFIX}/include)

# Set the link directories
link_directories(${THIRDPARTY_INSTALL_PREFIX}/lib)

# Build the third-party libraries
add_subdirectory(thirdparty)

# Build the C++ lirary Reaktoro
add_subdirectory(Reaktoro)

# Build the Python extension module PyReaktoro and the Python packages reaktoro
if(BUILD_PYTHON)
    add_subdirectory(ReaktoroPython)
endif()

# Build the Python package ireaktoro, the Reaktoro interpreter
if(BUILD_INTERPRETER)
    add_subdirectory(ReaktoroInterpreter)
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

# Create an install target for all built libraries
install(DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/
    DESTINATION "lib" COMPONENT development)

# Create an install target for all built binaries
if(BUILD_INTERPRETER)
    install(DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/
        DESTINATION "bin" COMPONENT applications)
endif()

# Create an install target for Reaktoro header files
install(DIRECTORY ${REAKTORO_SOURCE_DIR}
    DESTINATION "include" COMPONENT development
    FILES_MATCHING PATTERN "*.hpp" PATTERN "*.hxx")

# Package Reaktoro
include(PackageReaktoro)
