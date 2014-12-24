project(pyreaktor)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

find_package(PythonLibs 2.7 REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

find_package(Boost COMPONENTS python REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

file(GLOB_RECURSE HEADER_FILES *.hpp)
file(GLOB_RECURSE SOURCE_FILES *.cpp)

add_library(pyreaktor SHARED ${HEADER_FILES} ${SOURCE_FILES})
target_link_libraries(pyreaktor Reaktor ${Boost_LIBRARIES})
set_target_properties(pyreaktor PROPERTIES PREFIX "")