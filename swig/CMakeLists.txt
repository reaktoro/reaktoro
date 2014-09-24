find_package(SWIG 3.0 REQUIRED)

include(${SWIG_USE_FILE})

find_package(PythonLibs)

include_directories(${PYTHON_INCLUDE_PATH})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR})

set(CMAKE_CXX_FLAGS "-std=c++11")
set(CMAKE_SWIG_FLAGS "")

set_source_files_properties(Reaktor.i PROPERTIES CPLUSPLUS ON)

swig_add_module(reaktor python Reaktor.i ${SOURCE_FILES})
swig_link_libraries(reaktor ${PYTHON_LIBRARIES})