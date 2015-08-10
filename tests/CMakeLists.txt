#------------------------------------------------------------------------------
# Version Requirement
#------------------------------------------------------------------------------
cmake_minimum_required(VERSION 2.8)

#------------------------------------------------------------------------------
# Project Details
#------------------------------------------------------------------------------
project(tests)

#------------------------------------------------------------------------------
# Header and Source Files
#------------------------------------------------------------------------------
# Recursively collect all header files from the current directory
file(GLOB_RECURSE HEADER_FILES *.hpp)

# Recursively collect all source files from the current directory
file(GLOB_RECURSE SOURCE_FILES *.cpp)

#------------------------------------------------------------------------------
# Output Directories
#------------------------------------------------------------------------------
# Update the output path of the executable binaries produced by the project
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/tests)

#------------------------------------------------------------------------------
# Build Tests
#------------------------------------------------------------------------------
# Create an executable for each source file within the current directory
#foreach(sourcefile ${SOURCE_FILES})
#  get_filename_component(filename ${sourcefile} NAME_WE)
#  add_executable(${filename} ${sourcefile})
#  target_link_libraries(${filename} ${CMAKE_PROJECT_NAME})
#  add_test(NAME ${filename} COMMAND ${filename})
#endforeach()

#add_executable(ReaktoroTester ${SOURCE_FILES})
#target_link_libraries(ReaktoroTester ${CMAKE_PROJECT_NAME} ${THIRD_PARTY_LIBRARIES})
#add_test(NAME ReaktoroTester COMMAND ReaktoroTester)