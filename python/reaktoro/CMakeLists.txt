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

# Set the path where the python package is installed to CMAKE_INSTALL_PREFIX if not given
if(NOT DEFINED REAKTORO_PYTHON_INSTALL_PREFIX)
    set(REAKTORO_PYTHON_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
endif()

# Set the path where the dist's source is installed to CMAKE_SOURCE_DIR/dist if not given
if(NOT DEFINED REAKTORO_PYTHON_DIST_DIR)
    set(REAKTORO_PYTHON_DIST_DIR ${REAKTORO_PYTHON_INSTALL_PREFIX}/dist)
endif()

# If the path is already in Windows format (with backslashes), it can't be added directly
# to the string below, otherwise CMake will later complain about "Invalid escape sequence".
file(TO_CMAKE_PATH "${REAKTORO_PYTHON_INSTALL_PREFIX}" REAKTORO_PYTHON_INSTALL_PREFIX)
file(TO_CMAKE_PATH "${REAKTORO_PYTHON_DIST_DIR}" REAKTORO_PYTHON_DIST_DIR)

# Install the reaktoro python package using setuptools
install(CODE
"
    file(TO_NATIVE_PATH \"${REAKTORO_PYTHON_INSTALL_PREFIX}\" REAKTORO_PYTHON_INSTALL_PREFIX_NATIVE)
    file(TO_NATIVE_PATH \"${REAKTORO_PYTHON_DIST_DIR}\" REAKTORO_PYTHON_DIST_DIR_NATIVE)

    if(EXISTS ${CMAKE_CURRENT_BINARY_DIR}/../../lib/PyReaktoro.pdb)
        string(REPLACE .pyd .pdb REAKTORO_PDB_FILENAME \"${REAKTORO_PYTHON_MODULE_FILENAME}\")

        execute_process(
            COMMAND \${CMAKE_COMMAND} -E copy ../../lib/PyReaktoro.pdb \${REAKTORO_PDB_FILENAME}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    endif()

    execute_process(
        COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY}
            sdist --dist-dir=\${REAKTORO_PYTHON_DIST_DIR_NATIVE}
            bdist_wheel --dist-dir=\${REAKTORO_PYTHON_DIST_DIR_NATIVE}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

    execute_process(
        COMMAND ${PYTHON_EXECUTABLE} -m pip
            install reaktoro3 --prefix=\${REAKTORO_PYTHON_INSTALL_PREFIX_NATIVE}
            --find-links=\${REAKTORO_PYTHON_DIST_DIR_NATIVE}
        WORKING_DIRECTORY ${REAKTORO_PYTHON_INSTALL_PREFIX})
"
)
