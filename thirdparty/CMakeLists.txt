# Include the ExternalProject module
include(ExternalProject)

# Ignore all warning from the third-party libraries
if(CMAKE_COMPILER_IS_GNUCXX)
  add_definitions("-w")
endif()

# Download and install the pugixml library
ExternalProject_Add(PUGIXML
    PREFIX thirdparty
    GIT_REPOSITORY https://github.com/zeux/pugixml.git
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ${CMAKE_COMMAND} <SOURCE_DIR>/scripts
        -DCMAKE_GENERATOR=${CMAKE_GENERATOR}
        -DCMAKE_INSTALL_PREFIX=${THIRDPARTY_DIR}
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
        -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}
        -DCMAKE_INSTALL_BINDIR=${CMAKE_INSTALL_LIBDIR}
        -DBUILD_SHARED_LIBS=ON
)

# Download and install the Eigen library
ExternalProject_Add(EIGEN
    PREFIX thirdparty
    URL https://bitbucket.org/eigen/eigen/get/3.2.5.tar.bz2
    CMAKE_ARGS -DCMAKE_GENERATOR=${CMAKE_GENERATOR}
               -DCMAKE_INSTALL_PREFIX=${THIRDPARTY_DIR}
               -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    PATCH_COMMAND patch -p1 < ${CMAKE_SOURCE_DIR}/thirdparty/patches/Eigen.patch
)

# Build and install the cvode library
ExternalProject_Add(CVODE
    PREFIX thirdparty
    URL ${CMAKE_SOURCE_DIR}/thirdparty/cvode/cvode-2.8.0.tar.gz
    SOURCE_DIR thirdparty/src/cvode
    BINARY_DIR thirdparty/src/cvode-build
    STAMP_DIR thirdparty/src/cvode-stamp
    CMAKE_ARGS -DCMAKE_GENERATOR=${CMAKE_GENERATOR}
               -DCMAKE_INSTALL_PREFIX=${THIRDPARTY_DIR}
               -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
)

# Build and install the Phreeqc library
if(BUILD_PHREEQC)
    ExternalProject_Add(PHREEQC
        PREFIX thirdparty
        SOURCE_DIR ${CMAKE_SOURCE_DIR}/thirdparty/phreeqc
        CMAKE_ARGS -DCMAKE_GENERATOR=${CMAKE_GENERATOR}
                   -DCMAKE_INSTALL_PREFIX=${THIRDPARTY_DIR}
                   -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    )
endif(BUILD_PHREEQC)

# Build and install the Gems library
if(BUILD_GEMS)
    ExternalProject_Add(GEMS
        PREFIX thirdparty
        SOURCE_DIR ${CMAKE_SOURCE_DIR}/thirdparty/gems
        CMAKE_ARGS -DCMAKE_GENERATOR=${CMAKE_GENERATOR}
                   -DCMAKE_INSTALL_PREFIX=${THIRDPARTY_DIR}
                   -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    )
endif(BUILD_GEMS)

# Install the third-party dependencies (headers and libraries)
install(DIRECTORY ${THIRDPARTY_DIR}/include ${THIRDPARTY_DIR}/lib
    DESTINATION "/"
    COMPONENT development
    )
