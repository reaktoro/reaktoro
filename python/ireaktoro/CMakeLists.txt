# Copy the reaktoro directory to the binary directory
file(COPY . DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Configure the setup.py file
set(SETUP_PY_IN "${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in")
set(SETUP_PY    "${CMAKE_CURRENT_BINARY_DIR}/setup.py")

# Configure the setup.py script with the provided CMake variables
configure_file(${SETUP_PY_IN} ${SETUP_PY})

# Collect all files in the ireactoro directory
file(GLOB_RECURSE FILES ${CMAKE_CURRENT_SOURCE_DIR} *)

# Find PyInstaller, which is used to create an executable for ireaktoro
find_package(PyInstaller)

# Check if PyInstaller was found
if(PYINSTALLER_FOUND)

    # Set the string `pythonX.Y` where X and Y are the major and minor versions of python
    set(PYTHONXX "python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}")

    # Set the directory where the Reaktoro python packages have been built
    set(REAKTORO_PYTHON_PACKAGES_DIR
        "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${PYTHONXX}/site-packages")

    # Set the path to the file `ireaktoro.py`, the main file of the Reaktoro interpreter project
    set(IREAKTORO_PY "${CMAKE_CURRENT_SOURCE_DIR}/ireaktoro/ireaktoro.py")

    # Create a custom command that is executed whenever a file in FILES is changed
    add_custom_command(
        DEPENDS ${FILES}
        OUTPUT timestamp
        COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY} install --prefix=${CMAKE_BINARY_DIR}
        COMMAND ${PYINSTALLER} ${IREAKTORO_PY} --paths=${REAKTORO_PYTHON_PACKAGES_DIR} --hidden-import=numpy --exclude=dolfin --clean -y --distpath=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
        COMMAND ${CMAKE_COMMAND} -E touch timestamp
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

    # Create the ireaktoro target that depends on the stamp generated in the previous command
    add_custom_target(ireaktoro ALL DEPENDS timestamp)

endif()