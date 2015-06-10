# Copy the reaktoro directory to the binary directory
file(COPY reaktoro DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Configure the setup.py file
set(SETUP_PY_IN "${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in")
set(SETUP_PY    "${CMAKE_CURRENT_BINARY_DIR}/setup.py")

configure_file(${SETUP_PY_IN} ${SETUP_PY})

# Find python executable
find_program(PYTHON "python")

if(PYTHON)

    add_custom_target(reaktoro
        COMMENT "Building the Python package reaktoro."
        COMMAND ${PYTHON} ${SETUP_PY} build)

    install(CODE "execute_process(COMMAND ${PYTHON} ${SETUP_PY} install)")
endif()