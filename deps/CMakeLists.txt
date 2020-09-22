cmake_minimum_required(VERSION 3.9)

# The project that builds the dependencies of Reaktoro
project(deps)

# Set the cmake module path of the dependency build project
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../cmake")

# Use ccache to speed up repeated compilations
include(CCache)

# CMake commands to add external projects
include(ExternalProject)

# Set the list of compiler flags for the external projects
if(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set(CFLAGS "${CFLAGS} -fPIC")
    set(CXXFLAGS "${CXXFLAGS} -fPIC")
else()
    set(CFLAGS "${CMAKE_C_FLAGS}")
    set(CXXFLAGS "${CMAKE_CXX_FLAGS}")
endif()

if(NOT ${REAKTORO_DEPS_EXTRA_BUILD_ARGS} STREQUAL "")
    separate_arguments(REAKTORO_DEPS_EXTRA_BUILD_ARGS NATIVE_COMMAND "${REAKTORO_DEPS_EXTRA_BUILD_ARGS}")
endif()

# Set the common cmake arguments to all external projects
set(REAKTORO_DEPS_COMMON_INSTALL_ARGS
    -DCMAKE_GENERATOR=${CMAKE_GENERATOR}
    $<$<BOOL:CMAKE_GENERATOR_PLATFORM>:-DCMAKE_GENERATOR_PLATFORM=${CMAKE_GENERATOR_PLATFORM}>
    # -DCMAKE_INSTALL_PREFIX=${REAKTORO_DEPS_INSTALL_PREFIX}
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
    -DCMAKE_CXX_STANDARD=14
    -DCMAKE_CXX_STANDARD_REQUIRED=ON
    # Disables the cmake warning message `Manually-specified
     # variables were not used by the project` because of
     # unused arguments above
    --no-warn-unused-cli
    # Allows the specification of extra args via command-line
    ${REAKTORO_DEPS_EXTRA_BUILD_ARGS}
)

# Build Eigen library
ExternalProject_Add(EIGEN
    PREFIX         build
    SOURCE_DIR     eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG        master
    GIT_SHALLOW    TRUE
    GIT_PROGRESS   TRUE
    CMAKE_ARGS
        -DBUILD_TESTING=FALSE
        -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/install
)

# Build Optima library
ExternalProject_Add(OPTIMA
    PREFIX         build
    SOURCE_DIR     optima
    GIT_REPOSITORY https://github.com/allanleal/optima
    GIT_TAG        develop
    GIT_SHALLOW    TRUE
    GIT_PROGRESS   TRUE
    CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/install
        -DOPTIMA_BUILD_PYTHON=ON
        -DOPTIMA_BUILD_DEMOS=OFF
)

# Build autodiff library
ExternalProject_Add(AUTODIFF
    PREFIX         build
    SOURCE_DIR     autodiff
    GIT_REPOSITORY https://github.com/allanleal/autodiff
    GIT_TAG        develop
    GIT_SHALLOW    TRUE
    GIT_PROGRESS   TRUE
    CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/install
        -DAUTODIFF_BUILD_PYTHON=ON
        -DAUTODIFF_BUILD_EXAMPLES=OFF
        -DAUTODIFF_BUILD_TESTS=OFF
)

# Build and install the cvode library
ExternalProject_Add(CVODE
    PREFIX build
    URL ${CMAKE_CURRENT_SOURCE_DIR}/cvode/cvode-2.8.2.tar.gz
    CMAKE_ARGS ${REAKTORO_DEPS_COMMON_INSTALL_ARGS}
        -DEXAMPLES_ENABLE=OFF
        -DEXAMPLES_INSTALL=OFF)

# Build and install the Phreeqc library
ExternalProject_Add(PHREEQC
    PREFIX build
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/phreeqc
    CMAKE_ARGS ${REAKTORO_DEPS_COMMON_INSTALL_ARGS})

# Build and install the Gems library
ExternalProject_Add(GEMS
    PREFIX build
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/gems
    CMAKE_ARGS ${REAKTORO_DEPS_COMMON_INSTALL_ARGS})