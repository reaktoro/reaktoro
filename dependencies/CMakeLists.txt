# Include the ExternalProject module
include(ExternalProject)

# Ignore all warning from the dependencies
if(CMAKE_COMPILER_IS_GNUCXX)
  add_definitions("-w")
endif()

# Download and install the pugixml library
ExternalProject_Add(PUGIXML
    PREFIX thirdparty
    GIT_REPOSITORY https://github.com/zeux/pugixml.git
    UPDATE_COMMAND ""
    INSTALL_DIR ${CMAKE_INSTALL_PREFIX}
    CONFIGURE_COMMAND ${CMAKE_COMMAND} <SOURCE_DIR>/scripts
        -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_PREFIX}/include
        -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_PREFIX}/lib
        -DBUILD_SHARED_LIBS=ON
)

# Download and install the Eigen library
ExternalProject_Add(EIGEN
    PREFIX thirdparty
    URL https://bitbucket.org/eigen/eigen/get/3.2.5.tar.bz2
    INSTALL_DIR ${CMAKE_INSTALL_PREFIX}
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
               -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
)

# Build and install the cvode library
ExternalProject_Add(CVODE
    PREFIX thirdparty
    URL ${CMAKE_SOURCE_DIR}/dependencies/cvode/cvode-2.8.0.tar.gz
    SOURCE_DIR thirdparty/cvode
    BINARY_DIR thirdparty/cvode-build
    STAMP_DIR thirdparty/cvode-stamp
    INSTALL_DIR ${CMAKE_INSTALL_PREFIX}
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
               -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
)

# Build and install the Phreeqc library
if(BUILD_PHREEQC)
    ExternalProject_Add(PHREEQC
        PREFIX thirdparty
        SOURCE_DIR ${CMAKE_SOURCE_DIR}/dependencies/phreeqc
        INSTALL_DIR ${CMAKE_INSTALL_PREFIX}
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
                   -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    )
endif(BUILD_PHREEQC)

# Build and install the Gems library
if(BUILD_GEMS)
    ExternalProject_Add(GEMS
        PREFIX thirdparty
        SOURCE_DIR ${CMAKE_SOURCE_DIR}/dependencies/gems
        INSTALL_DIR ${CMAKE_INSTALL_PREFIX}
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
                   -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    )
endif(BUILD_GEMS)
