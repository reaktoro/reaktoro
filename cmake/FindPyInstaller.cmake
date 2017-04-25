# - Find the PyInstaller libraries
# This module finds if PyInstaller is installed, and sets the following variables:
#
#  PYINSTALLER                     - the PyInstaller executable
#  PYINSTALLER_FOUND               - the flag that indicates if PyInstaller was found
#  PYINSTALLER_VERSION             - the version of PyInstaller as a string
#  PYINSTALLER_VERSION_MAJOR       - the major version number of PyInstaller
#  PYINSTALLER_VERSION_MINOR       - the minor version number of PyInstaller
#  PYINSTALLER_VERSION_PATCH       - the patch version number of PyInstaller

# Find PyInstaller, which is used to create an executable out of ireaktoro project
if(NOT ${CMAKE_SYSTEM_NAME} MATCHES Windows)
    find_program(PYINSTALLER pyinstaller)
else()
    if(CMAKE_SIZEOF_VOID_P EQUAL 8)
        set(PYTHON_DIR "C:/Python27-w64/Scripts")
    else()
        set(PYTHON_DIR "C:/Python27-w32/Scripts")
    endif()
    find_program(PYINSTALLER pyinstaller PATH ${PYTHON_DIR} NO_DEFAULT_PATH)
endif()

# Check if PyInstaller was found
if(PYINSTALLER)
    set(PYINSTALLER_FOUND TRUE)
else()
    set(PYINSTALLER_FOUND FALSE)

    if(${PyInstaller_FIND_REQUIRED})
        message(FATAL_ERROR "Could not find required package PyInstaller.")
    else()
        message(STATUS "Could not find PyInstaller.")
    endif()

    message(STATUS "On Windows, make sure PyInstaller can be found at "
        "C:\\Python27-w32\\Scripts (for Win32 version) or "
        "C:\\Python27-w64\\Scripts (for Win64 version).")
endif()

# Set the PyInstaller version variables
if(PYINSTALLER_FOUND)
    # Get the version of PyInstaller
#    execute_process(COMMAND ${PYINSTALLER} "--version"
#                    OUTPUT_VARIABLE PYINSTALLER_VERSION)
#
    # Remove the endline character
#    string(REPLACE "\n" "" PYINSTALLER_VERSION ${PYINSTALLER_VERSION})
#
    # Replace dots with `;` to create a list of strings
#    string(REPLACE "." ";" PYINSTALLER_LIST ${PYINSTALLER_VERSION})
#
    # Append a 0 to the list to ensure that patch is 0 if not provided (e.g., 2.1)
#    list(APPEND PYINSTALLER_LIST 0)
#
    # Set the major, minor and patch version numbers
#    list(GET PYINSTALLER_LIST 0 PYINSTALLER_VERSION_MAJOR)
#    list(GET PYINSTALLER_LIST 1 PYINSTALLER_VERSION_MINOR)
#    list(GET PYINSTALLER_LIST 2 PYINSTALLER_VERSION_PATCH)
#
    # Print a message to indicate that PyInstaller was successfuly found
#    message(STATUS "Found PyInstaller version ${PYINSTALLER_VERSION} at ${PYINSTALLER}.")
    message(STATUS "Found PyInstaller.")
endif()
