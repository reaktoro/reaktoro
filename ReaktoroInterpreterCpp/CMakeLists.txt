# Require a certain version of cmake
cmake_minimum_required(VERSION 2.8)

# Recursively collect all header files from the current directory
file(GLOB_RECURSE HEADER_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.hpp *.hxx)

# Recursively collect all source files from the current directory
file(GLOB_RECURSE SOURCE_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

# Compile ReaktoroInterpreter into object files
add_library(ReaktoroInterpreterObject OBJECT ${HEADER_FILES} ${SOURCE_FILES})

# Compile Reaktoro interpreter into a shared library
add_library(ReaktoroInterpreter SHARED $<TARGET_OBJECTS:ReaktoroInterpreterObject>)

# Compile Reaktoro interpreter into an executable
add_executable(reaktoro $<TARGET_OBJECTS:ReaktoroInterpreterObject>)

# Specify that Reaktoro must be compiled first before the interpreter
add_dependencies(ReaktoroInterpreterObject ReaktoroObject)
add_dependencies(reaktoro ReaktoroObject)

# Link the library and executable against Reaktoro
target_link_libraries(ReaktoroInterpreter ReaktoroShared)
target_link_libraries(reaktoro ReaktoroStatic)