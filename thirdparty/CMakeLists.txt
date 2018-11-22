# Set cmake version requirement
cmake_minimum_required(VERSION 3.6)

# Set the name of the project
project(Reaktoro-thirdparty)

# Include the ExternalProject module
include(ExternalProject)

# Include the cmake variables with values for installation directories
include(GNUInstallDirs)

# Set the default build type to Release for the external dependencies
if(NOT CMAKE_BUILD_TYPE)
    # The build type selection for the external dependencies
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the common build type for the external dependencies." FORCE)

    # The build type options for the external dependencies
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release MinSizeRel RelWithDebInfo)
endif()

# Output a message detailing the build type used
message(STATUS "Using build type `${CMAKE_BUILD_TYPE}` for all external dependencies.")

# Set the path where the external dependencies should be installed
set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR})

# Set the list of compiler flags for the external projects
if(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set(CFLAGS "${CFLAGS} -fPIC")
    set(CXXFLAGS "${CXXFLAGS} -fPIC")
else()
    set(CFLAGS "${CMAKE_C_FLAGS}")
    set(CXXFLAGS "${CMAKE_CXX_FLAGS}")
endif()
 
separate_arguments(REAKTORO_THIRDPARTY_INSTALL_ARGS NATIVE_COMMAND "${REAKTORO_THIRDPARTY_INSTALL_ARGS}")

# Set the common cmake arguments to all external projects
set(REAKTORO_THIRDPARTY_COMMON_INSTALL_ARGS
    -DCMAKE_GENERATOR=${CMAKE_GENERATOR}
    $<$<BOOL:CMAKE_GENERATOR_PLATFORM>:-DCMAKE_GENERATOR_PLATFORM=${CMAKE_GENERATOR_PLATFORM}>
    -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DCMAKE_C_FLAGS=${CFLAGS}
    -DCMAKE_CXX_FLAGS=${CXXFLAGS}
    -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
    -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}
    -DCMAKE_INSTALL_BINDIR=${CMAKE_INSTALL_BINDIR}
    -DBUILD_SHARED_LIBS=OFF
    -DCMAKE_CXX_STANDARD=14
    -DCMAKE_CXX_STANDARD_REQUIRED=ON
    # Disables the cmake warning message `Manually-specified
     # variables were not used by the project` because of
     # unused arguments above
    --no-warn-unused-cli
    # Allows the specification of extra args via command-line
    ${REAKTORO_THIRDPARTY_INSTALL_ARGS}
    )

# Build and install the pugixml library
ExternalProject_Add(PUGIXML
    PREFIX build
    BUILD_ALWAYS TRUE
    URL ${CMAKE_CURRENT_SOURCE_DIR}/pugixml/pugixml-1.9.tar.gz
    CMAKE_ARGS ${REAKTORO_THIRDPARTY_COMMON_INSTALL_ARGS})

# Build and install the json library
ExternalProject_Add(JSON
    PREFIX build
    BUILD_ALWAYS TRUE
    URL https://github.com/nlohmann/json/releases/download/v3.1.2/json.hpp
    URL_HASH SHA256=fbdfec4b4cf63b3b565d09f87e6c3c183bdd45c5be1864d3fcb338f6f02c1733
    DOWNLOAD_DIR download/nlohmann
    DOWNLOAD_NO_EXTRACT TRUE
    DOWNLOAD_NO_PROGRESS TRUE
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ${CMAKE_COMMAND} -E copy_directory 
        ${CMAKE_CURRENT_BINARY_DIR}/download/nlohmann 
            ${CMAKE_INSTALL_PREFIX}/include/nlohmann
    CMAKE_ARGS ${REAKTORO_THIRDPARTY_COMMON_INSTALL_ARGS}
)

# Build and install the miniz library
ExternalProject_Add(MINIZ
    PREFIX build
    BUILD_ALWAYS TRUE
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/miniz
    CMAKE_ARGS ${REAKTORO_THIRDPARTY_COMMON_INSTALL_ARGS})

# Build and install the cvode library
ExternalProject_Add(CVODE
    PREFIX build
    BUILD_ALWAYS TRUE
    URL ${CMAKE_CURRENT_SOURCE_DIR}/cvode/cvode-2.8.2.tar.gz
    BUILD_ALWAYS TRUE
    CMAKE_ARGS ${REAKTORO_THIRDPARTY_COMMON_INSTALL_ARGS}
        -DEXAMPLES_ENABLE=OFF
        -DEXAMPLES_INSTALL=OFF)

# Build and install the Phreeqc library
ExternalProject_Add(PHREEQC
    PREFIX build
    BUILD_ALWAYS TRUE
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/phreeqc
    CMAKE_ARGS ${REAKTORO_THIRDPARTY_COMMON_INSTALL_ARGS})

# Build and install the Gems library
ExternalProject_Add(GEMS
    PREFIX build
    BUILD_ALWAYS TRUE
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/gems
    BUILD_ALWAYS TRUE
    CMAKE_ARGS ${REAKTORO_THIRDPARTY_COMMON_INSTALL_ARGS})

# Build and install the pybind11 library
ExternalProject_Add(PYBIND11
    PREFIX build
    BUILD_ALWAYS TRUE
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
    GIT_SHALLOW TRUE
    CMAKE_ARGS ${REAKTORO_THIRDPARTY_COMMON_INSTALL_ARGS}
        -DPYBIND11_TEST=OFF
)
