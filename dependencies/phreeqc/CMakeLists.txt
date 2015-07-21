# Require a certain version of cmake
cmake_minimum_required(VERSION 2.8)

# Set the name of the project
project(phreeqc)

add_definitions(-DPhreeqc=PHREEQC)

# Recursively collect all header files from the current directory
file(GLOB_RECURSE HEADER_FILES *.h)

# Recursively collect all source files from the current directory
file(GLOB_RECURSE SOURCE_FILES *.cpp *.cxx)

# Create the shared library using the collected source files
add_library(phreeqc OBJECT ${HEADER_FILES} ${SOURCE_FILES})
