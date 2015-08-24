# Require a certain version of cmake
cmake_minimum_required(VERSION 2.8)

# Set the name of the project
project(gems)

# Recursively collect all header files from the current directory
file(GLOB_RECURSE HEADER_FILES *.h)

# Recursively collect all source files from the current directory
file(GLOB_RECURSE SOURCE_FILES *.cpp)

# Set some necessary definitions
add_definitions(-DIPMGEMPLUGIN)
add_definitions(-DNOPARTICLEARRAY)

# Create the shared library using the collected source files
add_library(gems STATIC ${HEADER_FILES} ${SOURCE_FILES})

# Install the target shared library
install(TARGETS gems DESTINATION lib)

# Install the header files preserving the directory hierarchy
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} DESTINATION include FILES_MATCHING PATTERN "*.h")
