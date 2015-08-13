# Ensure Eigen uses row-major matrices for compatibility with the python wrappers
add_definitions("-DEIGEN_DEFAULT_TO_ROW_MAJOR")

# Include the parent directory of PyReaktoro
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/..)

# Find Python, Numpy and Boost.Python libraries
find_package(PythonLibs 2.7 REQUIRED)
find_package(Numpy 2.7 REQUIRED)
find_package(Boost COMPONENTS python REQUIRED)

# Specify Python, Numpy and Boost include directories
include_directories(${PYTHON_INCLUDE_DIRS})
include_directories(${NUMPY_INCLUDE_DIR})
include_directories(${Boost_INCLUDE_DIRS})

# Collect the header and source files of PyReaktoro
file(GLOB_RECURSE HEADER_FILES *.hpp)
file(GLOB_RECURSE SOURCE_FILES *.cpp)

# Suppress pedantic warnings related to ISO C++ for GNU compiler
if(${CMAKE_COMPILER_IS_GNUCXX})
    add_definitions(-Wno-pedantic)
endif()

# Create the PyReaktoro library
add_library(PyReaktoro SHARED ${HEADER_FILES} ${SOURCE_FILES})

# Link PyReaktoro against Boost.Python and Python libraries
target_link_libraries(PyReaktoro Reaktoro ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})

# Remove the `lib` prefix of the PyReaktoro library file name
set_target_properties(PyReaktoro PROPERTIES PREFIX "")

# Check if PyReaktoro library is output as `.dll` (Windows) and change it to `.pyd`
if(${CMAKE_SHARED_LIBRARY_SUFFIX} MATCHES ".dll")
    set_target_properties(PyReaktoro PROPERTIES SUFFIX ".pyd")
endif()

# Check if PyReaktoro library is output as `.dylib` (MacOS) and change it to `.so`
if(${CMAKE_SHARED_LIBRARY_SUFFIX} MATCHES ".dylib")
    set_target_properties(PyReaktoro PROPERTIES SUFFIX ".so")
endif()

# Install the target shared library
install(TARGETS PyReaktoro DESTINATION ${REAKTORO_INSTALL_LIBDIR})

# Copy the Python extension module PyReaktoro to the reaktoro/core directory
add_custom_command(TARGET PyReaktoro POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:PyReaktoro>
    ${CMAKE_BINARY_DIR}/python/reaktoro/core/$<TARGET_FILE_NAME:PyReaktoro>
    COMMENT "Copying the Python extension module PyReaktoro to the reaktoro/core directory."
)