# Require a certain version of cmake
cmake_minimum_required(VERSION 2.8)

# Recursively collect all header files from the current directory
file(GLOB_RECURSE HEADER_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.hpp)

# Recursively collect all source files from the current directory
file(GLOB_RECURSE SOURCE_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

# The name of the shared and static libraries
set(REAKTORO_SHARED_LIB ${PROJECT_NAME}${SUFFIX_SHARED_LIBS})
set(REAKTORO_STATIC_LIB ${PROJECT_NAME}${SUFFIX_STATIC_LIBS})

# Set the list of names of the third-party targets and libraries
set(THIRDPARTY_TARGETS PUGIXML EIGEN CVODE)
set(THIRDPARTY_LIBS pugixml sundials_cvode sundials_nvecserial)

if(BUILD_PHREEQC)
    set(THIRDPARTY_TARGETS ${THIRDPARTY_TARGETS} PHREEQC)
    set(THIRDPARTY_LIBS ${THIRDPARTY_LIBS} phreeqc)
endif()

if(BUILD_GEMS)
    set(THIRDPARTY_TARGETS ${THIRDPARTY_TARGETS} GEMS)
    set(THIRDPARTY_LIBS ${THIRDPARTY_LIBS} gems)
endif()

# Compile Reaktoro into object files
add_library(ReaktoroObject OBJECT ${HEADER_FILES} ${SOURCE_FILES})

# Specify the third-party dependencies of Reaktoro to ensure proper build sequence
add_dependencies(ReaktoroObject ${THIRDPARTY_TARGETS})

# Create both static and shared libraries out of the object files
add_library(ReaktoroShared SHARED $<TARGET_OBJECTS:ReaktoroObject>)
add_library(ReaktoroStatic STATIC $<TARGET_OBJECTS:ReaktoroObject>)

# Link the shared library with the third party libraries
target_link_libraries(ReaktoroStatic ${THIRDPARTY_LIBS})
target_link_libraries(ReaktoroShared ${THIRDPARTY_LIBS})

# Set the names of the static and shared Reaktoro libraries
set_target_properties(ReaktoroStatic PROPERTIES OUTPUT_NAME Reaktoro)
set_target_properties(ReaktoroShared PROPERTIES OUTPUT_NAME Reaktoro)

