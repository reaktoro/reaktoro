# Copy the reaktoro directory to the binary directory
file(COPY . DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Configure the setup.py file
set(SETUP_PY_IN "${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in")
set(SETUP_PY    "${CMAKE_CURRENT_BINARY_DIR}/setup.py")

# Configure the setup.py script with the provided CMake variables
configure_file(${SETUP_PY_IN} ${SETUP_PY})

# Collect all files in the ireactoro directory
file(GLOB_RECURSE FILES ${CMAKE_CURRENT_SOURCE_DIR} *)

# Create a custom command that builds the python package ireaktoro
add_custom_command(
    DEPENDS ${FILES}
    OUTPUT timestamp-build
    COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY} install --prefix=${CMAKE_BINARY_DIR}
    COMMAND ${CMAKE_COMMAND} -E touch timestamp-build
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Create the ireaktoro target that depends on the stamp generated in the previous command
add_custom_target(ireaktoro ALL DEPENDS timestamp-build)

# Set dependencies of ireaktoro target
add_dependencies(ireaktoro reaktoro)

# Find PyInstaller
find_package(PyInstaller)

include(FindPythonModule)

# Check if required python packages are installed
find_python_module(tabulate REQUIRED)

# Check if PyInstaller was found
if(NOT PYINSTALLER_FOUND)
    message(WARNING "PyInstaller was not found. "
        "Cannot generate an executable application for Reaktoro's interpreter. "
        "Install PyInstaller with `pip install pyinstaller`.")
endif()

# Check if PyInstaller was found
if(PYINSTALLER_FOUND)
    # Set the path to the installed dir of the Reaktoro python packages ireaktoro and reaktoro
    if(${CMAKE_SYSTEM_NAME} MATCHES Windows)
        set(REAKTORO_PYTHON_PACKAGES_DIR "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/Lib/site-packages")
    else()
        set(REAKTORO_PYTHON_PACKAGES_DIR "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages")
    endif()
    
    # Set the path to the file `ireaktoro.py`, the main file of the Reaktoro interpreter project
    set(IREAKTORO_PY "${REAKTORO_PYTHON_PACKAGES_DIR}/ireaktoro/ireaktoro.py")

    # Set the path to Reaktoro's icon
    set(REAKTORO_ICON "${REAKTORO_RESOURCES_DIR}/icons/reaktoro.ico")
    
    # Create a custom command that is executed whenever a file in FILES is changed
    add_custom_command(
        DEPENDS ${FILES} DEPENDS PyReaktoro DEPENDS reaktoro
        OUTPUT timestamp-interpreter
        COMMAND ${PYINSTALLER} ${IREAKTORO_PY}
            --paths=${REAKTORO_PYTHON_PACKAGES_DIR}
            --hidden-import=numpy
            --exclude=dolfin
            --clean -y
            --debug
            --distpath=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
            --icon=${REAKTORO_ICON}
        COMMAND ${CMAKE_COMMAND} -E copy_directory
            ${REAKTORO_DATABASES_DIR}
            ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/ireaktoro/databases
        COMMAND ${CMAKE_COMMAND} -E touch timestamp-interpreter
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

    # Create the target that depends on the stamp generated in the previous command
    add_custom_target(PyInstallerCommand ALL DEPENDS timestamp-interpreter)
    
    # Set the above target to depend on ireaktoro target
    add_dependencies(PyInstallerCommand ireaktoro)

    # Define the interpreter directory that is used on the next install target
    if(${CMAKE_SYSTEM_NAME} MATCHES Windows)
        # This causes all interpreter files/libraries to be located at bin
        set(REAKTORO_INTERPRETER_DIR ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/ireaktoro/)
    else()
        # This causes all interpreter files/libraries to be located at bin/ireaktoro
        set(REAKTORO_INTERPRETER_DIR ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/ireaktoro)
    endif()
endif()