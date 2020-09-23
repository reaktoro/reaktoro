# The name of the compiled python module file as defined by pybind11 used in setup.py.in
set(REAKTORO_PYTHON_MODULE_FILENAME reaktoro4py${PYTHON_MODULE_EXTENSION})

# The path in the binary dir where the python package is assembled
set(REAKTORO_PYTHON_PACKAGE_PATH ${CMAKE_CURRENT_BINARY_DIR}/reaktoro)

# Configure the setup.py file
set(SETUP_PY_IN ${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in)
set(SETUP_PY ${REAKTORO_PYTHON_PACKAGE_PATH}/setup.py) # ensure setup.py is inside the package dir
configure_file(${SETUP_PY_IN} ${SETUP_PY})

# Create a custom target to build the python package during build stage
add_custom_target(reaktoro-setuptools ALL
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/reaktoro ${REAKTORO_PYTHON_PACKAGE_PATH}
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:reaktoro4py> ${REAKTORO_PYTHON_PACKAGE_PATH}
    COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY} --quiet install --force --prefix=${CMAKE_BINARY_DIR}
    WORKING_DIRECTORY ${REAKTORO_PYTHON_PACKAGE_PATH})

# Set dependencies of the python package target
add_dependencies(reaktoro-setuptools reaktoro4py)

# Set the path where the python package is installed to CMAKE_INSTALL_PREFIX if not given
if(NOT DEFINED REAKTORO_PYTHON_INSTALL_PREFIX)
    set(REAKTORO_PYTHON_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
endif()

# If the path is already in Windows format (with backslashes), it can't be
# added directly to the string below. CMake complains about "Invalid escape sequence".
file(TO_CMAKE_PATH "${REAKTORO_PYTHON_INSTALL_PREFIX}" REAKTORO_PYTHON_INSTALL_PREFIX)

# Install the python package using setuptools
install(CODE
"
    file(TO_NATIVE_PATH \"${REAKTORO_PYTHON_INSTALL_PREFIX}\" REAKTORO_PYTHON_INSTALL_PREFIX_NATIVE)

    if(EXISTS ${CMAKE_CURRENT_BINARY_DIR}/../../lib/reaktoro4py.pdb)
        string(REPLACE .pyd .pdb REAKTORO_PDB_FILENAME \"${REAKTORO_PYTHON_MODULE_FILENAME}\")

        execute_process(
            COMMAND \${CMAKE_COMMAND} -E copy ../../lib/reaktoro4py.pdb \${REAKTORO_PDB_FILENAME}
            WORKING_DIRECTORY ${REAKTORO_PYTHON_PACKAGE_PATH})
    endif()

    execute_process(
        COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY} install --prefix=\${REAKTORO_PYTHON_INSTALL_PREFIX_NATIVE}
        WORKING_DIRECTORY ${REAKTORO_PYTHON_PACKAGE_PATH})
")
