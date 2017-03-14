# Require a certain version of cmake
cmake_minimum_required(VERSION 2.8)

# Recursively collect all header files from the current directory
file(GLOB_RECURSE HEADER_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.hpp *.hxx)

# Recursively collect all source files from the current directory
file(GLOB_RECURSE SOURCE_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

# The name of the shared and static libraries
set(REAKTORO_SHARED_LIB ${PROJECT_NAME}${SUFFIX_SHARED_LIBS})
set(REAKTORO_STATIC_LIB ${PROJECT_NAME}${SUFFIX_STATIC_LIBS})

# Set the list of names of the third-party targets and libraries
set(THIRDPARTY_TARGETS PUGIXML YAMLCPP MINIZ CVODE)
set(THIRDPARTY_LIBS pugixml yaml-cpp miniz sundials_cvode sundials_nvecserial)

if(LINK_PHREEQC)
    set(THIRDPARTY_TARGETS ${THIRDPARTY_TARGETS} PHREEQC)
    set(THIRDPARTY_LIBS ${THIRDPARTY_LIBS} phreeqc)
endif()

if(LINK_GEMS)
    set(THIRDPARTY_TARGETS ${THIRDPARTY_TARGETS} GEMS)
    set(THIRDPARTY_LIBS ${THIRDPARTY_LIBS} gems)
endif()

# Compile Reaktoro into object files
add_library(ReaktoroObject OBJECT ${HEADER_FILES} ${SOURCE_FILES})

# Apply cotire to ReaktotoObject target to improve compilation time
#cotire(ReaktoroObject)

# Specify the third-party dependencies of Reaktoro to ensure proper build sequence
add_dependencies(ReaktoroObject ${THIRDPARTY_TARGETS})

# Check if a shared library for Reaktoro must be built
if(BUILD_SHARED_LIBS)
    add_library(ReaktoroShared SHARED $<TARGET_OBJECTS:ReaktoroObject>)
    target_link_libraries(ReaktoroShared ${THIRDPARTY_LIBS})
    set_target_properties(ReaktoroShared PROPERTIES OUTPUT_NAME Reaktoro)
    install(TARGETS ReaktoroShared DESTINATION "lib" COMPONENT libraries)
endif()

# Check if a static library for Reaktoro must be built
if(BUILD_STATIC_LIBS)
    add_library(ReaktoroStatic STATIC $<TARGET_OBJECTS:ReaktoroObject>)
    target_link_libraries(ReaktoroStatic ${THIRDPARTY_LIBS})
    set_target_properties(ReaktoroStatic PROPERTIES OUTPUT_NAME Reaktoro)
    install(TARGETS ReaktoroStatic DESTINATION "lib" COMPONENT libraries)
endif()

# Check if an interpreter executable for Reaktoro must be built
if(BUILD_INTERPRETER)
    add_executable(ReaktoroInterpreter $<TARGET_OBJECTS:ReaktoroObject>)
    if(BUILD_STATIC_LIBS)
        message(STATUS "Interpreter will be linked against `static library`.")
        target_link_libraries(ReaktoroInterpreter ReaktoroStatic)
    else()
        message(STATUS "Interpreter will be linked against `shared library`.")
        target_link_libraries(ReaktoroInterpreter ReaktoroShared)
    endif()
    set_target_properties(ReaktoroInterpreter PROPERTIES OUTPUT_NAME reaktoro)
    install(TARGETS ReaktoroInterpreter DESTINATION "bin" COMPONENT applications)
endif()

# Create an install target for Reaktoro header files
install(DIRECTORY ${REAKTORO_SOURCE_DIR}
    DESTINATION "include" COMPONENT headers
    FILES_MATCHING PATTERN "*.hpp" PATTERN "*.hxx")

# Create an install target for Eigen header files
install(DIRECTORY ${REAKTORO_SOURCE_DIR}/Math/Eigen
    DESTINATION "include/Reaktoro/Math" COMPONENT headers)
