# Set cmake version requirement
cmake_minimum_required(VERSION 3.9)

# Set the name of the project
project(Reaktoro VERSION 1.0.0 LANGUAGES CXX)

# Set the cmake module path of the project
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

# Check if a conda environment is active
include(CondaRecognizer)

# Include this module to enable testing with CTest and dashboard submissions to CDash
include(CTest)

# Include the cmake variables with values for installation directories
include(GNUInstallDirs)

# Set the output directories of the built libraries and binaries
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})

# Define which components of Reaktoro to build
option(REAKTORO_BUILD_ALL         "Build everything." OFF)
option(REAKTORO_BUILD_DEMOS       "Build demos." OFF)
option(REAKTORO_BUILD_DOCS        "Build documentation." OFF)
option(REAKTORO_BUILD_INTERPRETER "Build the interpreter executable reaktoro." ON)
option(REAKTORO_BUILD_PYTHON      "Build the python wrappers and python package reaktoro." ON)
option(REAKTORO_BUILD_TESTS       "Build tests." OFF)

# Modify the REAKTORO_BUILD_* variables accordingly to BUILD_ALL
if(REAKTORO_BUILD_ALL MATCHES ON)
    set(REAKTORO_BUILD_DEMOS       ON)
    set(REAKTORO_BUILD_DOCS        ON)
    set(REAKTORO_BUILD_INTERPRETER ON)
    set(REAKTORO_BUILD_PYTHON      ON)
    set(REAKTORO_BUILD_TESTS       ON)
endif()

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

# Define which types of libraries to build
option(BUILD_SHARED_LIBS "Build shared libraries." ON)
option(BUILD_STATIC_LIBS "Build static libraries." OFF)

# Define custom options
option(ENABLE_TESTING "Enable testing." on)

# Set the default build type to Release
if(NOT CMAKE_BUILD_TYPE)
    # The build type selection for the project
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the build type for ${PROJECT_NAME}." FORCE)

    # The build type options for the project
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release MinSizeRel RelWithDebInfo)
endif()

# Set libraries to be compiled with position independent code mode (i.e., fPIC option in GNU compilers)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set the list of compiler flags for MSVC compiler
if(${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
    add_compile_options(
        /D_SCL_SECURE_NO_WARNINGS
        /D_CRT_SECURE_NO_WARNINGS=1
        /MP4
        /EHsc
        /D_SILENCE_TR1_NAMESPACE_DEPRECATION_WARNING
        /DNOMINMAX
    )
endif()

# Set the local directory where the third-party libraries are installed
set(REAKTORO_THIRDPARTY_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/thirdparty/install)

# Set path variables for the include and lib directories of the installed dependencies
set(REAKTORO_THIRDPARTY_INCLUDE_PATH ${REAKTORO_THIRDPARTY_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR})
set(REAKTORO_THIRDPARTY_LIBRARY_PATH ${REAKTORO_THIRDPARTY_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})

# Add the root path of this project to the cmake include directories
include_directories(${PROJECT_SOURCE_DIR})

# Add the include path of the installed external dependencies to the cmake include directories
include_directories(${REAKTORO_THIRDPARTY_INCLUDE_PATH})

# Add the library path of the installed external dependencies to the cmake link directories
link_directories(${REAKTORO_THIRDPARTY_LIBRARY_PATH})

# Find the Boost library 
find_package(Boost REQUIRED)

# Add the boost include path to the cmake include directories
include_directories(${Boost_INCLUDE_DIRS})

# Build the third-party libraries
add_subdirectory(thirdparty)

# Build the C++ library Reaktoro
add_subdirectory(Reaktoro)

# Build the Python extension module PyReaktoro and the Python packages reaktoro
if(REAKTORO_BUILD_PYTHON)
    add_subdirectory(python)
else()
    add_subdirectory(python EXCLUDE_FROM_ALL)
endif()

# Build the Python package ireaktoro, the Reaktoro interpreter
if(REAKTORO_BUILD_INTERPRETER)
    add_subdirectory(interpreter)
endif()

# Build the demonstration applications
if(REAKTORO_BUILD_DEMOS)
    add_subdirectory(demos)
else()
    add_subdirectory(demos EXCLUDE_FROM_ALL)
endif()

# Build the project documentation
if(REAKTORO_BUILD_DOCS)
    add_subdirectory(docs)
else()
    add_subdirectory(docs EXCLUDE_FROM_ALL)
endif()

# Build the tests
if(REAKTORO_BUILD_TESTS)
    add_subdirectory(tests)
else()
    add_subdirectory(tests EXCLUDE_FROM_ALL)
endif()

# Build the utilities
add_subdirectory(utilities EXCLUDE_FROM_ALL)

# Add target "python" for manual building of python wrappers, as `make python`, if REAKTORO_BUILD_PYTHON is OFF
add_custom_target(python
    COMMAND ${CMAKE_MAKE_PROGRAM}
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/python")
    
# Add target "interpreter" for manual building of interpreter, as `make interpreter`, if REAKTORO_BUILD_INTERPRETER is OFF
add_custom_target(interpreter
    COMMAND ${CMAKE_MAKE_PROGRAM}
    WORKING_DIRECTORY "${CURRENT_BINARY_DIR}/interpreter")

# Add target "demos" for manual building of demos, as `make demos`, if REAKTORO_BUILD_DEMOS is OFF
add_custom_target(demos
    COMMAND ${CMAKE_MAKE_PROGRAM}
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/demos")

# Add target "tests" for manual building of tests, as `make tests`, if REAKTORO_BUILD_TESTS is OFF
add_custom_target(tests
    COMMAND ${CMAKE_MAKE_PROGRAM}
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/tests")

# Add target "utilities" for manual building of utilities, as `make utilities`, if REAKTORO_BUILD_UTILITIES is OFF
add_custom_target(utilities
    COMMAND ${CMAKE_MAKE_PROGRAM}
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/utilities")

# Copy database files to binary directory
add_custom_target(databases ALL COMMENT "Copy database files to binary directory")
add_custom_command(
    TARGET databases POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
            ${PROJECT_SOURCE_DIR}/databases
            ${CMAKE_BINARY_DIR}/databases)

# Package Reaktoro
include(PackageReaktoro)
