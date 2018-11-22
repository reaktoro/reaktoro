# Include the ExternalProject module
include(ExternalProject)

# Set the list of compiler flags for the external projects
if(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set(CFLAGS "${CFLAGS} -fPIC")
    set(CXXFLAGS "${CXXFLAGS} -fPIC")
else()
    set(CFLAGS "${CMAKE_C_FLAGS}")
    set(CXXFLAGS "${CMAKE_CXX_FLAGS}")
endif()

separate_arguments(REAKTORO_THIRDPARTY_INSTALL_ARGS NATIVE_COMMAND "${REAKTORO_THIRDPARTY_INSTALL_ARGS}")
if (NOT ${THIRDPARTY_EXTRA_BUILD_ARGS} STREQUAL "")
    separate_arguments(THIRDPARTY_EXTRA_BUILD_ARGS NATIVE_COMMAND "${THIRDPARTY_EXTRA_BUILD_ARGS}")
endif()

# Set the common cmake arguments to all external projects
set(REAKTORO_THIRDPARTY_COMMON_INSTALL_ARGS
    -DCMAKE_GENERATOR=${CMAKE_GENERATOR}
    $<$<BOOL:CMAKE_GENERATOR_PLATFORM>:-DCMAKE_GENERATOR_PLATFORM=${CMAKE_GENERATOR_PLATFORM}>
    -DCMAKE_INSTALL_PREFIX=${REAKTORO_THIRDPARTY_INSTALL_PREFIX}
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
    ${THIRDPARTY_COMMON_ARGS}
    )

# Build and install the pugixml library
ExternalProject_Add(PUGIXML
    PREFIX build
    URL ${CMAKE_CURRENT_SOURCE_DIR}/pugixml/pugixml-1.9.tar.gz
    CMAKE_ARGS ${REAKTORO_THIRDPARTY_COMMON_INSTALL_ARGS})

# Build and install the json library
ExternalProject_Add(JSON
    PREFIX build
    URL https://github.com/nlohmann/json/releases/download/v3.1.2/json.hpp
    URL_HASH SHA256=fbdfec4b4cf63b3b565d09f87e6c3c183bdd45c5be1864d3fcb338f6f02c1733
    DOWNLOAD_DIR download/nlohmann
    DOWNLOAD_NO_EXTRACT TRUE
    DOWNLOAD_NO_PROGRESS TRUE
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_BINARY_DIR}/download/nlohmann
            ${REAKTORO_THIRDPARTY_INSTALL_PREFIX}/include/nlohmann
    CMAKE_ARGS ${REAKTORO_THIRDPARTY_COMMON_INSTALL_ARGS}
)

# Build and install the miniz library
ExternalProject_Add(MINIZ
    PREFIX build
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/miniz
    CMAKE_ARGS ${REAKTORO_THIRDPARTY_COMMON_INSTALL_ARGS})

# Build and install the cvode library
ExternalProject_Add(CVODE
    PREFIX build
    URL ${CMAKE_CURRENT_SOURCE_DIR}/cvode/cvode-2.8.2.tar.gz
    CMAKE_ARGS ${REAKTORO_THIRDPARTY_COMMON_INSTALL_ARGS}
        -DEXAMPLES_ENABLE=OFF
        -DEXAMPLES_INSTALL=OFF)

# Build and install the Phreeqc library
ExternalProject_Add(PHREEQC
    PREFIX build
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/phreeqc
    CMAKE_ARGS ${REAKTORO_THIRDPARTY_COMMON_INSTALL_ARGS})

# Build and install the Gems library
ExternalProject_Add(GEMS
    PREFIX build
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/gems
    CMAKE_ARGS ${REAKTORO_THIRDPARTY_COMMON_INSTALL_ARGS})
