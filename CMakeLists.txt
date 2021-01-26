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

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)

# Check if a conda environment is active
include(CondaAware)

# Include the cmake variables with values for installation directories
include(GNUInstallDirs)

# Set the output directories of the built libraries and binaries
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})

# Define which components of Reaktoro to build
option(REAKTORO_BUILD_ALL         "Build everything." OFF)
option(REAKTORO_BUILD_EXAMPLES    "Build examples." OFF)
option(REAKTORO_BUILD_DOCS        "Build documentation." OFF)
option(REAKTORO_BUILD_INTERPRETER "Build the interpreter executable reaktoro." OFF)
option(REAKTORO_BUILD_PYTHON      "Build the python wrappers and python package reaktoro." ON)
option(REAKTORO_BUILD_TESTS       "Build tests." OFF)

# Modify the REAKTORO_BUILD_* variables accordingly to BUILD_ALL
if(REAKTORO_BUILD_ALL MATCHES ON)
    set(REAKTORO_BUILD_DOCS        ON)
    set(REAKTORO_BUILD_EXAMPLES    ON)
    set(REAKTORO_BUILD_INTERPRETER ON)
    set(REAKTORO_BUILD_PYTHON      ON)
    set(REAKTORO_BUILD_TESTS       ON)
endif()

# Define is Reaktoro should be built linking against openlibm instead of system's default libm
option(REAKTORO_USE_OPENLIBM "Build linking with openlibm." OFF)

# Define if shared library should be build instead of static.
option(BUILD_SHARED_LIBS "Build shared libraries." ON)

# Define custom options
option(ENABLE_TESTING "Enable testing." on)

# Set the default build type to Release
if(NOT CMAKE_BUILD_TYPE)
    # The build type selection for the project
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the build type for ${PROJECT_NAME}." FORCE)

    # The build type options for the project
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

# The path where the external dependency libraries are installed
set(REAKTORO_DEPS_PATH ${PROJECT_SOURCE_DIR}/deps/build/install)

# The paths where the public and private external dependency libraries are installed
set(REAKTORO_DEPS_PATH_PUBLIC ${REAKTORO_DEPS_PATH}/public)
set(REAKTORO_DEPS_PATH_PRIVATE ${REAKTORO_DEPS_PATH}/private)

# The paths where the include directory of the external dependencies are installed
set(REAKTORO_INCLUDE_PATH_DEPS
    ${REAKTORO_DEPS_PATH_PUBLIC}/include
    ${REAKTORO_DEPS_PATH_PRIVATE}/include)

# The paths where the lib directory of the external dependencies are installed
set(REAKTORO_DEPS_LIBRARY_PATH
    ${REAKTORO_DEPS_PATH_PUBLIC}/lib
    ${REAKTORO_DEPS_PATH_PRIVATE}/lib)

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

# Build the Python package ireaktoro, the Reaktoro interpreter
if(REAKTORO_BUILD_INTERPRETER)
    add_subdirectory(interpreter)
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

# Build the utilities
add_subdirectory(utilities EXCLUDE_FROM_ALL)

# Process sub-directory scripts
add_subdirectory(scripts)

# Add target "python" for manual building of python wrappers, as `make python`, if REAKTORO_BUILD_PYTHON is OFF
# add_custom_target(python
#     COMMAND ${CMAKE_MAKE_PROGRAM}
#     WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/python")

# Add target "interpreter" for manual building of interpreter, as `make interpreter`, if REAKTORO_BUILD_INTERPRETER is OFF
add_custom_target(interpreter
    COMMAND ${CMAKE_MAKE_PROGRAM}
    WORKING_DIRECTORY "${CURRENT_BINARY_DIR}/interpreter")

# Add target "examples" for manual building of examples, as `make examples`, if REAKTORO_BUILD_EXAMPLES is OFF
add_custom_target(examples
    COMMAND ${CMAKE_MAKE_PROGRAM}
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/demos")

# Add target "utilities" for manual building of utilities, as `make utilities`, if REAKTORO_BUILD_UTILITIES is OFF
add_custom_target(utilities
    COMMAND ${CMAKE_MAKE_PROGRAM}
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/utilities")

# Add target "deps" for updating the external dependencies
add_custom_target(deps
    COMMAND ${CMAKE_COMMAND} -P deps/install
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")

# Add target "optima" for building and installing Optima
add_custom_target(optima
    COMMAND ${CMAKE_COMMAND} --build . --parallel
    COMMAND ${CMAKE_COMMAND} --install .
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/deps/build/etc/src/OPTIMA-build")

# Package Reaktoro
include(PackageReaktoro)

# Install the cmake config files that permit users to use find_package(Reaktoro)
include(ReaktoroInstallCMakeConfigFiles)
