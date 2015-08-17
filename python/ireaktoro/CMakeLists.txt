# Copy the reaktoro directory to the binary directory
file(COPY . DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Configure the setup.py file
set(SETUP_PY_IN "${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in")
set(SETUP_PY    "${CMAKE_CURRENT_BINARY_DIR}/setup.py")

# Configure the setup.py script with the provided CMake variables
configure_file(${SETUP_PY_IN} ${SETUP_PY})

# Find python executable
find_program(PYTHON "python")

# Check if python was found
if(${PYTHON_NOTFOUND})
    message(FATAL_ERROR "Python executable is needed to build ireaktoro target.")
endif()

# Collect all files in the ireactoro directory
file(GLOB_RECURSE FILES ${CMAKE_CURRENT_SOURCE_DIR} *)

# Create a custom command that is executed whenever a file in FILES is changed
add_custom_command(
    DEPENDS ${FILES}
    OUTPUT timestamp
    COMMAND ${PYTHON} ${SETUP_PY} install --prefix=${CMAKE_CURRENT_BINARY_DIR}/install
    COMMAND ${CMAKE_COMMAND} -E touch timestamp
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Create the ireaktoro target that depends on the stamp generated in the previous command
add_custom_target(ireaktoro ALL
    DEPENDS timestamp)

# Install target for the Python ireaktoro package
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/install/
    DESTINATION "/"
    COMPONENT interfaces)

