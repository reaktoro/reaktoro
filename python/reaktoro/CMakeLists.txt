# Copy the reaktoro directory to the binary directory
file(COPY . DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Define a variable with the name of the compiled python module file as defined by pybind11 (which is used in the file setup.py.in)
set(REAKTORO_PYTHON_MODULE_FILENAME PyReaktoro${PYTHON_MODULE_EXTENSION})

# Configure the setup.py file
set(SETUP_PY_IN ${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in)
set(SETUP_PY ${CMAKE_CURRENT_BINARY_DIR}/setup.py)

# Configure the setup.py script with the provided CMake variables
configure_file(${SETUP_PY_IN} ${SETUP_PY})

# Create a custom target to build reaktoro python package during build stage
add_custom_target(reaktoro ALL
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:PyReaktoro> 
        ${CMAKE_CURRENT_BINARY_DIR}/$<TARGET_FILE_NAME:PyReaktoro>
    COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY} --quiet install --prefix=${CMAKE_BINARY_DIR}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Set dependencies of reaktoro target
add_dependencies(reaktoro PyReaktoro)

# Install the reaktoro python package using setuptools
install(CODE
"
    execute_process(
        COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY} install --prefix=\${CMAKE_INSTALL_PREFIX} 
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
"
)
