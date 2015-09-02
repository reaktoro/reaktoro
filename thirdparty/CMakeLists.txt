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
    -DCMAKE_INSTALL_PREFIX=${THIRDPARTY_DIR}
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DCMAKE_C_FLAGS=${CFLAGS}
    -DCMAKE_CXX_FLAGS=${CXXFLAGS}
    -DCMAKE_INSTALL_INCLUDEDIR=include
    -DCMAKE_INSTALL_LIBDIR=lib
    -DCMAKE_INSTALL_BINDIR=bin
    -DBUILD_SHARED_LIBS=OFF)

# Download and install the pugixml library
ExternalProject_Add(PUGIXML
    PREFIX thirdparty
    GIT_REPOSITORY https://github.com/zeux/pugixml.git
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ${CMAKE_COMMAND} <SOURCE_DIR>/scripts ${THIRDPARTY_COMMON_ARGS})

# Download and install the Eigen library
ExternalProject_Add(EIGEN
    PREFIX thirdparty
    URL https://bitbucket.org/eigen/eigen/get/3.2.5.tar.bz2
    CMAKE_ARGS ${THIRDPARTY_COMMON_ARGS}
    PATCH_COMMAND patch -p1 < ${CMAKE_SOURCE_DIR}/thirdparty/patches/Eigen.patch)

# Build and install the cvode library
ExternalProject_Add(CVODE
    PREFIX thirdparty
    URL ${CMAKE_SOURCE_DIR}/thirdparty/cvode/cvode-2.8.0.tar.gz
    SOURCE_DIR thirdparty/src/cvode
    BINARY_DIR thirdparty/src/cvode-build
    STAMP_DIR thirdparty/src/cvode-stamp
    CMAKE_ARGS ${THIRDPARTY_COMMON_ARGS})

# Build and install the Phreeqc library
if(BUILD_PHREEQC)
    ExternalProject_Add(PHREEQC
        PREFIX thirdparty
        SOURCE_DIR ${CMAKE_SOURCE_DIR}/thirdparty/phreeqc
        CMAKE_ARGS ${THIRDPARTY_COMMON_ARGS})
endif()

# Build and install the Gems library
if(BUILD_GEMS)
    ExternalProject_Add(GEMS
        PREFIX thirdparty
        SOURCE_DIR ${CMAKE_SOURCE_DIR}/thirdparty/gems
        CMAKE_ARGS ${THIRDPARTY_COMMON_ARGS})
endif()

# Install the headers-only third-party library Eigen
install(DIRECTORY ${THIRDPARTY_DIR}/include/eigen3
    DESTINATION "include"
    COMPONENT development)