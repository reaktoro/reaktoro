# Configure the setup.py.in and version.py.in files
configure_file(setup.py.in ${CMAKE_CURRENT_BINARY_DIR}/setup.py)
configure_file(version.py.in ${CMAKE_CURRENT_BINARY_DIR}/reaktoro/version.py)

# Create pybind11-stubgen command string depending if this program has been found or not
if(PYBIND11_STUBGEN)
    set(STUBGEN_COMMAND1 ${CMAKE_COMMAND} -E env "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}" ${PYBIND11_STUBGEN} --ignore-all-errors reaktoro)
    set(STUBGEN_COMMAND2 ${CMAKE_COMMAND} -E copy_directory stubs/reaktoro/ reaktoro)
else()
    set(STUBGEN_COMMAND1 "")  # do nothing when it it's time to generate python stubs for reaktoro in the target below
    set(STUBGEN_COMMAND2 "")  # do nothing when it it's time to copy the generated stub directory
endif()

# Create a custom target to build the python package during build stage
add_custom_target(reaktoro-setuptools ALL
    DEPENDS reaktoro4py
    COMMAND ${CMAKE_COMMAND} -E rm -rf build  # remove any previously created build dir (due to `python setup.py install` below) to ensure fresh rebuild since any modified python files are not overwritten even when --force option is used
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/reaktoro reaktoro
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:reaktoro4py> reaktoro
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:Reaktoro> reaktoro
    COMMAND ${STUBGEN_COMMAND1}
    COMMAND ${STUBGEN_COMMAND2}
    COMMAND ${Python_EXECUTABLE} setup.py --quiet build --force
    COMMAND ${CMAKE_COMMAND} -E rm reaktoro/$<TARGET_FILE_NAME:reaktoro4py>
    COMMAND ${CMAKE_COMMAND} -E rm reaktoro/$<TARGET_FILE_NAME:Reaktoro>
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Ensure the path where the python package is installed is not empty
if(NOT DEFINED REAKTORO_PYTHON_INSTALL_PREFIX)
    file(TO_NATIVE_PATH ${CMAKE_INSTALL_PREFIX} REAKTORO_PYTHON_INSTALL_PREFIX)
endif()

# Create an install target for the reaktoro python package
install(CODE "execute_process(COMMAND pip install . --no-deps --force-reinstall --ignore-installed --no-compile --prefix=${REAKTORO_PYTHON_INSTALL_PREFIX} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})")
