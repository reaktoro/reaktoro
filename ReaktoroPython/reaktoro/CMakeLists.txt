# Copy the reaktoro directory to the binary directory
file(COPY . DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Configure the setup.py file
set(SETUP_PY_IN "${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in")
set(SETUP_PY    "${CMAKE_CURRENT_BINARY_DIR}/setup.py")

# Configure the setup.py script with the provided CMake variables
configure_file(${SETUP_PY_IN} ${SETUP_PY})

# Collect all files in the reactoro directory
file(GLOB_RECURSE FILES ${CMAKE_CURRENT_SOURCE_DIR} *)

# Create a custom command that is executed whenever a file in FILES is changed
#add_custom_command(
#    DEPENDS ${FILES}
#    OUTPUT timestamp
#    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:PyReaktoro> 
#        ${CMAKE_CURRENT_BINARY_DIR}/$<TARGET_FILE_NAME:PyReaktoro>
#    COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY} build
#    COMMAND ${CMAKE_COMMAND} -E touch timestamp
#    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
#
#add_custom_command(
#    DEPENDS ${FILES}
#    OUTPUT timestamp
#    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:PyReaktoro> 
#        ${CMAKE_CURRENT_BINARY_DIR}/$<TARGET_FILE_NAME:PyReaktoro>
#    COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY} install --prefix=${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
#    COMMAND ${CMAKE_COMMAND} -E touch timestamp
#    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})


add_custom_command(
    DEPENDS ${FILES} DEPENDS PyReaktoro
    OUTPUT timestamp
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:PyReaktoro> 
        ${CMAKE_CURRENT_BINARY_DIR}/$<TARGET_FILE_NAME:PyReaktoro>
    COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY} install --prefix=${CMAKE_BINARY_DIR}
    COMMAND ${CMAKE_COMMAND} -E touch timestamp
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    
# Create the reaktoro target that depends on the stamp generated in the previous command
add_custom_target(reaktoro ALL DEPENDS timestamp)

# Set dependencies of reaktoro target
add_dependencies(reaktoro PyReaktoro)

# Define the install target
install(CODE "
    execute_process(
        COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY} install --prefix=${CMAKE_INSTALL_PREFIX} 
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    ")