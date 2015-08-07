# Copy the reaktoro directory to the binary directory
file(COPY . DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Configure the setup.py file
set(SETUP_PY_IN "${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in")
set(SETUP_PY    "${CMAKE_CURRENT_BINARY_DIR}/setup.py")

# Configure the setup.py script with the provided CMake variables
configure_file(${SETUP_PY_IN} ${SETUP_PY})

# Find python executable
find_program(PYTHON "python")

# Install the Python reaktoro project
if(PYTHON)
    install(CODE "execute_process(COMMAND ${PYTHON} ${SETUP_PY} install --prefix=${CMAKE_INSTALL_PREFIX})")
endif()