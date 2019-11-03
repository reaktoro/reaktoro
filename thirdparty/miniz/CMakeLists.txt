# Require a certain version of cmake
cmake_minimum_required(VERSION 2.8)

# Set the name of the project
project(miniz)

# Recursively collect all header files from the current directory
set(HEADER_FILES miniz.h zip_file.hpp)

# Recursively collect all source files from the current directory
set(SOURCE_FILES miniz.c zip_file.cpp)

# Set the c++11 as the standard for GNU and clang compilers
if(${CMAKE_CXX_COMPILER_ID} MATCHES "GNU" OR ${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c11 -Wno-implicit-function-declaration")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif()

# Create the shared library using the collected source files
add_library(miniz STATIC ${HEADER_FILES} ${SOURCE_FILES})

# Install the target shared library
install(TARGETS miniz DESTINATION lib)

# Install the header files preserving the directory hierarchy
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} DESTINATION include
    FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")
