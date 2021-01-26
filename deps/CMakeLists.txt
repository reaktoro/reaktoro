cmake_minimum_required(VERSION 3.17)

# The project that builds the dependencies of Reaktoro
project(deps)

# Set the cmake module path of the dependency build project
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../cmake")

# Use ccache to speed up repeated compilations
include(CCache)

# CMake commands to add external projects
include(ExternalProject)

# Set the list of compiler flags for the external projects
if(WIN32)
    set(CFLAGS "${CMAKE_C_FLAGS} /MP")
    set(CXXFLAGS "${CMAKE_CXX_FLAGS} /MP") # Enable parallel build if MSVC is used with /MP
else()
    set(CFLAGS "${CFLAGS} -fPIC")
    set(CXXFLAGS "${CXXFLAGS} -fPIC")
endif()

if(NOT ${REAKTORO_DEPS_EXTRA_BUILD_ARGS} STREQUAL "")
    separate_arguments(REAKTORO_DEPS_EXTRA_BUILD_ARGS NATIVE_COMMAND "${REAKTORO_DEPS_EXTRA_BUILD_ARGS}")
endif()

# Set the common cmake arguments to all external projects
set(REAKTORO_DEPS_COMMON_INSTALL_ARGS
    -DCMAKE_GENERATOR=${CMAKE_GENERATOR}
        $<$<BOOL:CMAKE_GENERATOR_PLATFORM>:-DCMAKE_GENERATOR_PLATFORM=${CMAKE_GENERATOR_PLATFORM}>
    -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DCMAKE_C_COMPILER_LAUNCHER=${CMAKE_C_COMPILER_LAUNCHER}
    -DCMAKE_CXX_COMPILER_LAUNCHER=${CMAKE_CXX_COMPILER_LAUNCHER}
    -DCMAKE_C_FLAGS=${CFLAGS}
    -DCMAKE_CXX_FLAGS=${CXXFLAGS}
    -DCMAKE_INSTALL_INCLUDEDIR=include
    -DCMAKE_INSTALL_LIBDIR=lib
    -DCMAKE_INSTALL_BINDIR=bin
    -DBUILD_SHARED_LIBS=OFF
    -DCMAKE_CXX_STANDARD=17
    -DCMAKE_CXX_STANDARD_REQUIRED=ON
    -DCMAKE_INSTALL_MESSAGE=LAZY # Print "Installing" but not "Up-to-date messages" during installation step.
    --no-warn-unused-cli # Disables the cmake warning message `Manually-specified variables were not used by the project` because of unused arguments above
    ${REAKTORO_DEPS_EXTRA_BUILD_ARGS} # Allows the specification of extra args via command-line
)

# Ensure the dependencies handled by git submodules are initialized
execute_process(COMMAND git submodule update --init --recursive
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

# Build Eigen library
ExternalProject_Add(EIGEN
    PREFIX etc
    URL ${CMAKE_CURRENT_SOURCE_DIR}/eigen/eigen-master-9a4d04c05fbaad6e6de750b13cfeab9edeb71c90.tar.gz
    CMAKE_ARGS
        ${REAKTORO_DEPS_COMMON_INSTALL_ARGS}
        -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/install/public
        -DBUILD_TESTING=FALSE
)

# Build Optima library
ExternalProject_Add(OPTIMA
    PREFIX etc
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/optima
    DEPENDS EIGEN
    BUILD_ALWAYS ON
    CMAKE_ARGS
        ${REAKTORO_DEPS_COMMON_INSTALL_ARGS}
        -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/install/public
        -DOPTIMA_BUILD_PYTHON=ON
        -DOPTIMA_BUILD_DEMOS=OFF
)

# Build autodiff library
ExternalProject_Add(AUTODIFF
    PREFIX etc
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/autodiff
    DEPENDS EIGEN
    BUILD_ALWAYS ON
    CMAKE_ARGS
        ${REAKTORO_DEPS_COMMON_INSTALL_ARGS}
        -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/install/public
        -DAUTODIFF_BUILD_PYTHON=ON
        -DAUTODIFF_BUILD_EXAMPLES=OFF
        -DAUTODIFF_BUILD_TESTS=OFF
)

# Build and install the cvode library
ExternalProject_Add(CVODE
    PREFIX etc
    URL ${CMAKE_CURRENT_SOURCE_DIR}/cvode/cvode-2.8.2.tar.gz
    CMAKE_ARGS
        ${REAKTORO_DEPS_COMMON_INSTALL_ARGS}
        -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/install/private
        -DEXAMPLES_ENABLE=OFF
        -DEXAMPLES_INSTALL=OFF
)

# Build and install the Phreeqc library
ExternalProject_Add(PHREEQC
    PREFIX etc
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/phreeqc
    CMAKE_ARGS
        ${REAKTORO_DEPS_COMMON_INSTALL_ARGS}
        -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/install/private
)

# # Build and install the Gems library
# ExternalProject_Add(GEMS
#     PREFIX etc
#     SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/gems
#     CMAKE_ARGS
#         ${REAKTORO_DEPS_COMMON_INSTALL_ARGS}
#         -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/install/private
# )
