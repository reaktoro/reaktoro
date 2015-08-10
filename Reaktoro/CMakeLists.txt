# Require a certain version of cmake
cmake_minimum_required(VERSION 2.8)

# Recursively collect all header files from the current directory
file(GLOB_RECURSE HEADER_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.hpp)

# Recursively collect all source files from the current directory
file(GLOB_RECURSE SOURCE_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

# Remove Phreeqc files if not to be built
if(NOT BUILD_PHREEQC)
    set(PHREEQC_HEADER_FILES "Interfaces/Phreeqc.hpp" "Interfaces/internal/PhreeqcUtils.hpp")
    set(PHREEQC_SOURCE_FILES "Interfaces/Phreeqc.cpp" "Interfaces/internal/PhreeqcUtils.cpp")
    list(REMOVE_ITEM HEADER_FILES ${PHREEQC_HEADER_FILES})
    list(REMOVE_ITEM SOURCE_FILES ${PHREEQC_SOURCE_FILES})
endif()

# Remove Gems files if not to be built
if(NOT BUILD_GEMS)
    list(REMOVE_ITEM HEADER_FILES "Interfaces/Gems.hpp")
    list(REMOVE_ITEM SOURCE_FILES "Interfaces/Gems.cpp")
endif()

# The name of the shared and static libraries
set(REAKTORO_SHARED_LIB ${PROJECT_NAME}${SUFFIX_SHARED_LIBS})
set(REAKTORO_STATIC_LIB ${PROJECT_NAME}${SUFFIX_STATIC_LIBS})

# Create the shared library using the collected source files
add_library(${REAKTORO_SHARED_LIB} SHARED ${HEADER_FILES} ${SOURCE_FILES})

# Specify the dependencies of the Reaktoro library to ensure proper build steps
add_dependencies(${REAKTORO_SHARED_LIB} PUGIXML EIGEN CVODE)

# Link the shared library with the third party libraries
target_link_libraries(${REAKTORO_SHARED_LIB} pugixml)
target_link_libraries(${REAKTORO_SHARED_LIB} sundials_cvode)
target_link_libraries(${REAKTORO_SHARED_LIB} sundials_nvecserial)

# Link the shared library with the optional Phreeqc library
if(BUILD_PHREEQC)
    add_dependencies(${REAKTORO_SHARED_LIB} PHREEQC)
    target_link_libraries(${REAKTORO_SHARED_LIB} phreeqc)
endif()

# Link the shared library with the optional Gems library
if(BUILD_GEMS)
    add_dependencies(${REAKTORO_SHARED_LIB} GEMS)
    target_link_libraries(${REAKTORO_SHARED_LIB} gems)
endif()

# Install the target shared library
install(TARGETS ${REAKTORO_SHARED_LIB} DESTINATION lib)

# Install the header files preserving the directory hierarchy
install(DIRECTORY ${REAKTORO_SOURCE_DIR} DESTINATION include FILES_MATCHING PATTERN "*.hpp" PATTERN "*.hxx")
