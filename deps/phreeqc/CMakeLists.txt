# Require a certain version of cmake
cmake_minimum_required(VERSION 2.8)

# Set the name of the project
project(phreeqc)

# Replace the word Phreeqc by PHREEQC
add_definitions(-DPhreeqc=PHREEQC)

# Enable std::unique_ptr with Visual Studio 2019 (CXX17).
add_compile_definitions($<$<CXX_COMPILER_ID:MSVC>:_HAS_AUTO_PTR_ETC=1>)

# Recursively collect all header files from the current directory
file(GLOB_RECURSE HEADER_FILES *.h)

# Recursively collect all source files from the current directory
file(GLOB_RECURSE SOURCE_FILES *.cpp *.cxx)

# Include directories for proper compilation
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/common)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/PhreeqcKeywords)

# Create the shared library using the collected source files
add_library(phreeqc STATIC ${HEADER_FILES} ${SOURCE_FILES})

# Install the target shared library
install(TARGETS phreeqc DESTINATION lib)

# Install the header files preserving the directory hierarchy
install(FILES ${HEADER_FILES} DESTINATION include/phreeqc)
