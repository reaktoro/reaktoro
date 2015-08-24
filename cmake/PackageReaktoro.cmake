###############################################################################
# Description: This cmake script creates a package for Reaktoro.
# Author: Allan Leal
# Date: 24 August 2015
###############################################################################

# Set the default CPack generators for Linux
if(NOT CPACK_GENERATOR AND ${CMAKE_SYSTEM_NAME} MATCHES Linux)
    set(CPACK_GENERATOR "TGZ;DEB;RPM")
endif()

# Configure the Reaktoro package target
set(CPACK_PACKAGE_NAME "Reaktoro")
set(CPACK_PACKAGE_CONTACT "Allan Leal (allan.leal@erdw.ethz.ch)")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "A unified framework for modeling chemically reactive systems")
set(CPACK_PACKAGE_VERSION ${REAKTORO_VERSION})
set(CPACK_PACKAGE_VERSION_MAJOR ${REAKTORO_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${REAKTORO_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${REAKTORO_VERSION_MICRO})
set(CPACK_PACKAGE_INSTALL_DIRECTORY "Reaktoro")

# Configure the package components (CPACK_COMPONENT_${COMPONENT_NAME_ALL_CAPS}_GROUP).
set(CPACK_COMPONENTS_ALL applications development)

# More descriptive names for each of the components, and component groups
set(CPACK_COMPONENT_APPLICATIONS_DISPLAY_NAME "Applications")
set(CPACK_COMPONENT_APPLICATIONS_DESCRIPTION "The applications developed using Reaktoro.")

set(CPACK_COMPONENT_DEVELOPMENT_DISPLAY_NAME "Development")
set(CPACK_COMPONENT_DEVELOPMENT_DESCRIPTION "The Reaktoro C++ headers and libraries for application development. This also install the Python wrappers for Reaktoro.")

# Set the license file for output during the installation
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/COPYING")

###############################################################################
# Debian Options
###############################################################################
set(CPACK_DEBIAN_PACKAGE_HOMEPAGE "www.reaktoro.org")

###############################################################################
# NSIS Options
###############################################################################
# Allow Reaktoro to be uninstalled before installed
set(CPACK_NSIS_ENABLE_UNINSTALL_BEFORE_INSTALL ON)

# Allow installation of Reaktoro to change the environment variable PATH
set(CPACK_NSIS_MODIFY_PATH ON)

# Set the icon used inside the installer wizard and the icons for the installer and un-installer
# set(CPACK_PACKAGE_ICON "${CMAKE_SOURCE_DIR}/hw.bmp")
# set(CPACK_NSIS_MUI_ICON "${CMAKE_SOURCE_DIR}/hw.ico")
# set(CPACK_NSIS_MUI_UNIICON "${CMAKE_SOURCE_DIR}/hw.ico")

# Set the programs displayed under the Start menu shortcut
# We do not use CPACK_PACKAGE_EXECUTABLES because it expects the binaries to be in a
# directory named "bin" (CMake bug tracker ID: 0007829).
set(CPACK_NSIS_MENU_LINKS "ireaktoro" "iReaktoro")

include(CPack)
