# - Find the Numpy libraries
# This module finds if Numpy is installed, and sets the following variables:
#
#  NUMPY_FOUND               - the flag that indicates if Numpy was found
#  NUMPY_INCLUDE_DIR         - path to the Numpy include directory
#  NUMPY_VERSION             - the version of Numpy as a string `major.minor.patch`
#  NUMPY_VERSION_MAJOR       - the major version number of Numpy
#  NUMPY_VERSION_MINOR       - the minor version number of Numpy
#  NUMPY_VERSION_PATCH       - the patch version number of Numpy

# Make sure the variable `PYTHON_EXECUTABLE` is initialised before proceeding with finding numpy
if(${Numpy_FIND_REQUIRED})
    find_package(PythonInterp REQUIRED)
else()
    find_package(PythonInterp)
endif()

# Execute python and print the include directory of numpy and its version
if(PYTHONINTERP_FOUND)
    execute_process(COMMAND ${PYTHON_EXECUTABLE}
      "-c" "import numpy; print numpy.get_include(); print numpy.version.version; print numpy.version.version.split('.')[0]; print numpy.version.version.split('.')[1]; print numpy.version.version.split('.')[2];"
      OUTPUT_VARIABLE NUMPY_OUTPUT_STRING
      RESULT_VARIABLE NUMPY_NOT_FOUND
      ERROR_VARIABLE  NUMPY_ERROR_VALUE)
else()
    message(FATAL_ERROR "Searching for Numpy requires Python executable, which was not found. Please install Python.")
endif()

# Check if numpy was found and act approapriately
if(${NUMPY_NOT_FOUND} STREQUAL "0")
    set(NUMPY_FOUND TRUE)
else()
    set(NUMPY_FOUND FALSE)
    set(NUMPY_INCLUDE_DIR -NOTFOUND)
    if(${Numpy_FIND_REQUIRED})
        message(FATAL_ERROR "Searching for Numpy resulted in the failure:\n${NUMPY_ERROR_VALUE}")
    endif()
    return()
endif()

# Replace newline symbols to `;` to create a list of strings
string(REGEX REPLACE "\n" ";" NUMPY_LIST ${NUMPY_OUTPUT_STRING})

# For every specific entry of the list `NUMPY_LIST`, set the approapriate variable
list(GET NUMPY_LIST 0 NUMPY_INCLUDE_DIR)
list(GET NUMPY_LIST 1 NUMPY_VERSION)
list(GET NUMPY_LIST 2 NUMPY_VERSION_MAJOR)
list(GET NUMPY_LIST 3 NUMPY_VERSION_MINOR)
list(GET NUMPY_LIST 4 NUMPY_VERSION_PATCH)

# Print a nice message to indicate that numpy was successfuly found
message(STATUS "Found Numpy version: ${NUMPY_VERSION}")

mark_as_advanced(NUMPY_INCLUDE_DIR)
