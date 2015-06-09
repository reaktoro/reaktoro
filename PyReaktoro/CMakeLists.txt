project(PyReaktoro)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/reaktoro/base)

find_package(PythonLibs 2.7 REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

find_package(Numpy 2.7 REQUIRED QUIET)
include_directories(${NUMPY_INCLUDE_DIR})

find_package(Boost COMPONENTS python REQUIRED QUIET)
include_directories(${Boost_INCLUDE_DIRS})

file(GLOB_RECURSE HEADER_FILES *.hpp)
file(GLOB_RECURSE SOURCE_FILES *.cpp)

# Suppress pedantic warnings related to ISO C++
add_definitions(-Wno-pedantic)

add_library(PyReaktoro SHARED ${HEADER_FILES} ${SOURCE_FILES})
target_link_libraries(PyReaktoro Reaktoro ${Boost_LIBRARIES})
set_target_properties(PyReaktoro PROPERTIES PREFIX "")