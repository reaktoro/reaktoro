# Include the ExternalProject module
include(ExternalProject)

# Set the list of compiler flags for the external projects
if(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set(CFLAGS "-fPIC")
    set(CXXFLAGS "-fPIC")
endif()

# Set the common cmake arguments to all external projects
set(THIRDPARTY_COMMON_ARGS
    -DCMAKE_GENERATOR=${CMAKE_GENERATOR}
    -DCMAKE_INSTALL_PREFIX=${THIRDPARTY_INSTALL_PREFIX}
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DCMAKE_C_FLAGS=${CFLAGS}
    -DCMAKE_CXX_FLAGS=${CXXFLAGS}
    -DCMAKE_INSTALL_INCLUDEDIR=include
    -DCMAKE_INSTALL_LIBDIR=lib
    -DCMAKE_INSTALL_BINDIR=bin
    -DBUILD_SHARED_LIBS=OFF
     # Disables the cmake warning message `Manually-specified
     # variables were not used by the project` because of
     # unused arguments above
    --no-warn-unused-cli
    )

# Build and install the pugixml library
ExternalProject_Add(PUGIXML
    PREFIX build
    URL ${CMAKE_CURRENT_SOURCE_DIR}/pugixml/pugixml-1.6.tar.gz
    SOURCE_DIR build/src/pugixml
    BINARY_DIR build/src/pugixml-build
    STAMP_DIR build/src/pugixml-stamp
    CONFIGURE_COMMAND ${CMAKE_COMMAND} <SOURCE_DIR>/scripts ${THIRDPARTY_COMMON_ARGS})

# Build and install the yaml-cpp library
ExternalProject_Add(YAMLCPP
    PREFIX build
    URL ${CMAKE_CURRENT_SOURCE_DIR}/yaml-cpp/yaml-cpp-0.5.3.zip
    SOURCE_DIR build/src/yaml-cpp
    BINARY_DIR build/src/yaml-cpp-build
    STAMP_DIR build/src/yaml-cpp-stamp
    CMAKE_ARGS ${THIRDPARTY_COMMON_ARGS}
        -DYAML_CPP_BUILD_TOOLS=OFF
        -DYAML_CPP_BUILD_CONTRIB=OFF)

# Build and install the miniz library
ExternalProject_Add(MINIZ
    PREFIX build
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/miniz
    CMAKE_ARGS ${THIRDPARTY_COMMON_ARGS})

# Build and install the cvode library
ExternalProject_Add(CVODE
    PREFIX build
    URL ${CMAKE_CURRENT_SOURCE_DIR}/cvode/cvode-2.8.2.tar.gz
    SOURCE_DIR build/src/cvode
    BINARY_DIR build/src/cvode-build
    STAMP_DIR build/src/cvode-stamp
    CMAKE_ARGS ${THIRDPARTY_COMMON_ARGS}
        -DEXAMPLES_ENABLE=OFF
        -DEXAMPLES_INSTALL=OFF)

# Build and install the Phreeqc library
if(LINK_PHREEQC)
    ExternalProject_Add(PHREEQC
        PREFIX build
        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/phreeqc
        CMAKE_ARGS ${THIRDPARTY_COMMON_ARGS})
endif()

# Build and install the Gems library
if(LINK_GEMS)
    ExternalProject_Add(GEMS
        PREFIX build
        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/gems
        CMAKE_ARGS ${THIRDPARTY_COMMON_ARGS})
endif()
