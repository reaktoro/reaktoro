###############################################################################
# Description: This cmake script install Reaktoro library and its python
#   wrappers and packages.
# Author: Allan Leal
# Date: 24 August 2015
###############################################################################

# Install all built libraries
install(DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/
    DESTINATION "lib"
    COMPONENT development)

# Install all built runtime binaries
install(DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/
    DESTINATION "bin"
    COMPONENT applications)

# Install the header-only third-party library Eigen
install(DIRECTORY ${THIRDPARTY_DIR}/include/eigen3
    DESTINATION "include"
    COMPONENT development)

# Install the header files preserving the directory hierarchy
install(DIRECTORY ${REAKTORO_SOURCE_DIR}
    DESTINATION "include"
    COMPONENT development
    FILES_MATCHING PATTERN "*.hpp" PATTERN "*.hxx")