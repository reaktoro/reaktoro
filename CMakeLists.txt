# Set cmake version requirement
cmake_minimum_required(VERSION 3.6)

# Set the name of the project
project(Reaktoro)

# Set the cmake module path of the project
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Include cotire module to speed compilation times
include(cotire)

# Include this module to enable testing with CTest and dashboard submissions to CDash
include(CTest)

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
option(BUILD_INTERPRETER "Build the interpreter executable reaktoro." ON)
option(BUILD_PYTHON      "Build the python wrappers and python package reaktoro." OFF)
option(BUILD_TESTS       "Build tests." OFF)

# Set the CTest option BUILD_TESTING to the value of BUILD_TESTS
set(BUILD_TESTING BUILD_TEST)

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

# Set libraries to be compiled with position independent code mode (i.e., fPIC option in GNU compilers)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Set the list of compiler flags for GNU compiler
if(${CMAKE_CXX_COMPILER_ID} MATCHES "GNU")
    add_compile_options(-std=c++11)
endif()

# Set the list of compiler flags for Clang compiler
if(${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
    add_compile_options(-std=c++11)
endif()

# Set the list of compiler flags for Intel compiler
if(${CMAKE_CXX_COMPILER_ID} MATCHES "Intel")
    add_compile_options(-std=c++11)
endif()

# Set the list of compiler flags for MSVC compiler
if(${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
    add_compile_options("/W0 -D_SCL_SECURE_NO_WARNINGS /MP4")
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

# Build the C++ library Reaktoro
add_subdirectory(Reaktoro)

# Build the Python extension module PyReaktoro and the Python packages reaktoro
if(BUILD_PYTHON)
    add_subdirectory(python)
endif()

# Build the Python package ireaktoro, the Reaktoro interpreter
if(BUILD_INTERPRETER)
    add_subdirectory(interpreter)
endif()

# Build the demonstration applications
if(BUILD_DEMOS)
    add_subdirectory(demos)
else()
    add_subdirectory(demos EXCLUDE_FROM_ALL)
endif()

# Build the project documentation
if(BUILD_DOCS)
    add_subdirectory(docs)
else()
    add_subdirectory(docs EXCLUDE_FROM_ALL)
endif()

# Build the tests
if(BUILD_TESTS)
    add_subdirectory(tests)
else()
    add_subdirectory(tests EXCLUDE_FROM_ALL)
endif()

# Build the utilities
if(BUILD_UTILITIES)
    add_subdirectory(utilities)
else()
    add_subdirectory(utilities EXCLUDE_FROM_ALL)
endif()

# Add target "python" for manual building of python wrappers, as `make python`, if BUILD_PYTHON is OFF
add_custom_target(python
    COMMAND ${CMAKE_MAKE_PROGRAM}
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/python")

# Add target "interpreter" for manual building of interpreter, as `make interpreter`, if BUILD_INTERPRETER is OFF
add_custom_target(interpreter
    COMMAND ${CMAKE_MAKE_PROGRAM}
    WORKING_DIRECTORY "${CURRENT_BINARY_DIR}/interpreter")

# Add target "demos" for manual building of demos, as `make demos`, if BUILD_DEMOS is OFF
add_custom_target(demos
    COMMAND ${CMAKE_MAKE_PROGRAM}
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/demos")

# Add target "tests" for manual building of tests, as `make tests`, if BUILD_TESTS is OFF
add_custom_target(tests
    COMMAND ${CMAKE_MAKE_PROGRAM}
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/tests")

# Add target "utilities" for manual building of utilities, as `make utilities`, if BUILD_UTILITIES is OFF
add_custom_target(utilities
    COMMAND ${CMAKE_MAKE_PROGRAM}
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/utilities")

# Copy database files to binary directory
add_custom_target(databases ALL COMMENT "Copy database files to binary directory")
add_custom_command(
    TARGET databases POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
            ${CMAKE_SOURCE_DIR}/databases
            ${CMAKE_BINARY_DIR}/databases)

# Package Reaktoro
include(PackageReaktoro)
