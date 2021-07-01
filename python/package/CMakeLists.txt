# The name of the compiled python module file as defined by pybind11 used in setup.py.in
set(REAKTORO_PYTHON_MODULE_FILENAME reaktoro4py${PYTHON_MODULE_EXTENSION})

# The path in the binary dir where the python package is assembled before it is built
set(REAKTORO_PYTHON_PACKAGE_PATH ${CMAKE_CURRENT_BINARY_DIR}/reaktoro)

# Configure the setup.py file
set(SETUP_PY_IN ${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in)
set(SETUP_PY ${CMAKE_CURRENT_BINARY_DIR}/setup.py)
configure_file(${SETUP_PY_IN} ${SETUP_PY})

# This is needed in Windows so that the next step produces the python package in CMAKE_BINARY_DIR\Lib\site-packages\reaktoro
file(TO_NATIVE_PATH ${CMAKE_BINARY_DIR} CMAKE_BINARY_DIR_NATIVE)

# Create a custom target to build the python package during build stage
add_custom_target(reaktoro-setuptools ALL
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/reaktoro ${CMAKE_CURRENT_BINARY_DIR}/reaktoro
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:reaktoro4py> ${CMAKE_CURRENT_BINARY_DIR}/reaktoro
    COMMAND ${CMAKE_COMMAND} -E env
        "PATH=${REAKTORO_PATH}"
        "PYTHONPATH=${REAKTORO_PYTHONPATH}"
            pybind11-stubgen --ignore-invalid=all --no-setup-py --log-level=ERROR reaktoro
    COMMAND ${CMAKE_COMMAND} -E copy_directory stubs/reaktoro-stubs/ reaktoro
    COMMAND ${CMAKE_COMMAND} -E rm -rf ${CMAKE_CURRENT_BINARY_DIR}/build  # remove build dir created by previous `python setup.py install` commands (see next) to ensure fresh rebuild since changed python files are not overwritten even with --force option
    COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY} --quiet install --force --no-compile --prefix=installed/$<CONFIG>
    COMMAND ${CMAKE_COMMAND} -E rm ${CMAKE_CURRENT_BINARY_DIR}/reaktoro/$<TARGET_FILE_NAME:reaktoro4py>
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Set dependencies of the python package target
add_dependencies(reaktoro-setuptools reaktoro4py)

# Ensure the path where the python package is installed is not empty
if(NOT DEFINED REAKTORO_PYTHON_INSTALL_PREFIX)
    set(REAKTORO_PYTHON_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
endif()

# Create an install target for the reaktoro python package
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/installed/$<CONFIG>/
    DESTINATION ${REAKTORO_PYTHON_INSTALL_PREFIX})
