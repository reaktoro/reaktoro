# The path where cmake config files are installed
set(REAKTORO_INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake)

install(EXPORT ReaktoroTargets
    FILE ReaktoroTargets.cmake
    NAMESPACE Reaktoro::
    DESTINATION ${REAKTORO_INSTALL_CONFIGDIR}
    COMPONENT cmake)

include(CMakePackageConfigHelpers)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/ReaktoroConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion)

configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/ReaktoroConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/ReaktoroConfig.cmake
    INSTALL_DESTINATION ${REAKTORO_INSTALL_CONFIGDIR}
    PATH_VARS REAKTORO_INSTALL_CONFIGDIR)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/ReaktoroConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/ReaktoroConfigVersion.cmake
    DESTINATION ${REAKTORO_INSTALL_CONFIGDIR})
