# Ensure Eigen uses row-major matrices for compatibility with the python wrappers
add_definitions("-DEIGEN_DEFAULT_TO_ROW_MAJOR")

# Include the parent directory of PyReaktoro
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/..)

find_package(PythonLibs 2.7 REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

find_package(Numpy 2.7 REQUIRED)
include_directories(${NUMPY_INCLUDE_DIR})

find_package(Boost COMPONENTS python REQUIRED QUIET)
include_directories(${Boost_INCLUDE_DIRS})

file(GLOB_RECURSE HEADER_FILES *.hpp)
file(GLOB_RECURSE SOURCE_FILES *.cpp)

# Suppress pedantic warnings related to ISO C++ for GNU compiler
if(${CMAKE_COMPILER_IS_GNUCXX})
    add_definitions(-Wno-pedantic)
endif()

add_library(PyReaktoro SHARED ${HEADER_FILES} ${SOURCE_FILES})
target_link_libraries(PyReaktoro Reaktoro ${Boost_LIBRARIES})
set_target_properties(PyReaktoro PROPERTIES PREFIX "")

# Install the target shared library
install(TARGETS PyReaktoro DESTINATION lib)

# Copy the Python extension module PyReaktoro to the reaktoro/core directory
add_custom_command(TARGET PyReaktoro POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:PyReaktoro>
    ${CMAKE_BINARY_DIR}/python/reaktoro/core/$<TARGET_FILE_NAME:PyReaktoro>
    COMMENT "Copying the Python extension module PyReaktoro to the reaktoro/core directory."
)