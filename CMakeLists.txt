# Set cmake version requirement
cmake_minimum_required(VERSION 3.17)

# Set the cmake module path of the project
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Use ccache to speed up repeated compilations
include(CCache)

# Use CMakeRC to enable easy and convenient access to resources
include(CMakeRC)

# Set the name of the projectw
project(Reaktoro VERSION 2.0 LANGUAGES CXX)

# Check if a conda environment is active
include(CondaAware)

# Include the cmake variables with values for installation directories
include(GNUInstallDirs)

# Define which components of Reaktoro to build
option(REAKTORO_BUILD_ALL         "Build everything." OFF)
option(REAKTORO_BUILD_EXAMPLES    "Build examples." OFF)
option(REAKTORO_BUILD_DOCS        "Build documentation." OFF)
option(REAKTORO_BUILD_PYTHON      "Build the python wrappers and python package reaktoro." ON)
option(REAKTORO_BUILD_TESTS       "Build tests." ON)

# Modify the REAKTORO_BUILD_* variables accordingly to BUILD_ALL
if(REAKTORO_BUILD_ALL MATCHES ON)
    set(REAKTORO_BUILD_DOCS        ON)
    set(REAKTORO_BUILD_EXAMPLES    ON)
    set(REAKTORO_BUILD_INTERPRETER ON)
    set(REAKTORO_BUILD_PYTHON      ON)
    set(REAKTORO_BUILD_TESTS       ON)
endif()

# Define is Reaktoro should be built linking against openlibm instead of system's default libm
option(REAKTORO_ENABLE_OPENLIBM "Build linking with openlibm." OFF)

# Define if shared library should be build instead of static.
option(BUILD_SHARED_LIBS "Build shared libraries." ON)

# Set the default build type to Release
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to Release as none was specified.")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the build type for ${PROJECT_NAME}." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release MinSizeRel RelWithDebInfo)
endif()

# Enable parallel build if MSVC is used
add_compile_options($<$<CXX_COMPILER_ID:MSVC>:/MP>)

# Set the list of compiler flags for MSVC compiler
if(${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
    add_compile_options(
        /D_SCL_SECURE_NO_WARNINGS
        /D_CRT_SECURE_NO_WARNINGS=1
        /EHsc
        /D_SILENCE_TR1_NAMESPACE_DEPRECATION_WARNING
        /DNOMINMAX
        /Zc:__cplusplus  # This is needed to ensure __cplusplus is replaced with a correct value (e.g. 201703L) instead of fixed 199711L (see more on https://docs.microsoft.com/bs-cyrl-ba/cpp/build/reference/zc-cplusplus?view=vs-2019)
    )
endif()

# The relative path where site-packages directory is located for current platform
execute_process(
    COMMAND "${PYTHON_EXECUTABLE}" -c "if True:
        from distutils import sysconfig as sc
        print(sc.get_python_lib(prefix='', plat_specific=True))"
    OUTPUT_VARIABLE PYTHON_SITE_PACKAGES_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE)

# The directory where the python package reaktoro is built
set(REAKTORO_PYTHON_PACKAGE_DIR ${CMAKE_BINARY_DIR}/python/package)

# Define variables REAKTORO_PATH and REAKTORO_PYTHONPATH used to set environment variables PATH and PYTHONPATH.
if(WIN32)
    set(REAKTORO_PATH "$<TARGET_FILE_DIR:Reaktoro>")
    set(REAKTORO_PYTHONPATH "${REAKTORO_PYTHON_PACKAGE_DIR}")
    file(TO_NATIVE_PATH ${REAKTORO_PATH} REAKTORO_PATH)
    file(TO_NATIVE_PATH ${REAKTORO_PYTHONPATH} REAKTORO_PYTHONPATH)
    set(REAKTORO_PATH "${REAKTORO_PATH}\;$ENV{PATH}")
    set(REAKTORO_PYTHONPATH "${REAKTORO_PYTHONPATH}\;$ENV{PYTHONPATH}")
else()
    set(REAKTORO_PATH "$ENV{PATH}")
    set(REAKTORO_PYTHONPATH "${REAKTORO_PYTHON_PACKAGE_DIR}:$ENV{PYTHONPATH}")
endif()

# Find all Reaktoro dependencies
include(ReaktoroFindDeps)

# Build the C++ resource library containing embedded databases
add_subdirectory(databases)

# Build the C++ library Reaktoro
add_subdirectory(Reaktoro)

# Build the Python package reaktoro
if(REAKTORO_BUILD_PYTHON)
    add_subdirectory(python)
endif()

# Build the examples
if(REAKTORO_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

# Build the project documentation
if(REAKTORO_BUILD_DOCS)
    add_subdirectory(docs)
endif()

# Build the tests
if(REAKTORO_BUILD_TESTS)
    add_subdirectory(tests)
endif()

# Process sub-directory scripts
add_subdirectory(scripts)

# Build the utilities
add_subdirectory(utilities)

# Package Reaktoro
include(PackageReaktoro)

# Install the cmake config files that permit users to use find_package(Reaktoro)
include(ReaktoroInstallCMakeConfigFiles)
