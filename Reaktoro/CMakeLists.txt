# Recursively collect all header files from the current directory
file(GLOB_RECURSE HEADER_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.hpp *.hxx)

# Recursively collect all source files from the current directory
file(GLOB_RECURSE SOURCE_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

# The name of the shared and static libraries
set(REAKTORO_SHARED_LIB ${PROJECT_NAME}${SUFFIX_SHARED_LIBS})
set(REAKTORO_STATIC_LIB ${PROJECT_NAME}${SUFFIX_STATIC_LIBS})

# Set the list of names of the third-party targets and libraries
set(THIRDPARTY_TARGETS PUGIXML MINIZ CVODE JSON)
set(THIRDPARTY_LIBS pugixml miniz sundials_cvode sundials_nvecserial)

if(LINK_PHREEQC)
    set(THIRDPARTY_TARGETS ${THIRDPARTY_TARGETS} PHREEQC)
    set(THIRDPARTY_LIBS ${THIRDPARTY_LIBS} phreeqc)
endif()

if(LINK_GEMS)
    set(THIRDPARTY_TARGETS ${THIRDPARTY_TARGETS} GEMS)
    set(THIRDPARTY_LIBS ${THIRDPARTY_LIBS} gems)
endif()

# Enable automatic creation of a module definition (.def) file for a SHARED library on Windows.
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS TRUE)

# Add the library path of the external dependencies to link directories
# Note: cmake v3.13 is needed to use target_link_directories command!
link_directories(${REAKTORO_THIRDPARTY_LIBRARY_PATH})

# Build a Reaktoro C++ library
add_library(Reaktoro ${HEADER_FILES} ${SOURCE_FILES})

# Specify the external dependencies of Reaktoro to ensure proper build sequence
add_dependencies(Reaktoro ${THIRDPARTY_TARGETS})

# Add the include paths to the Reaktoro target
target_include_directories(Reaktoro
    PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    PRIVATE ${REAKTORO_THIRDPARTY_INCLUDE_PATH})

# Link Reaktoro library against external dependencies
target_link_libraries(Reaktoro
    PRIVATE ${THIRDPARTY_LIBS}
    PUBLIC Boost::boost)

# Install Reaktoro C++ library
install(TARGETS Reaktoro
    EXPORT ReaktoroTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT libraries
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT libraries
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT libraries)

# Install Reaktoro header files
install(DIRECTORY ${PROJECT_SOURCE_DIR}/Reaktoro
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} COMPONENT headers
    PATTERN "*.cpp" EXCLUDE
    PATTERN "CMakeLists.txt" EXCLUDE)
