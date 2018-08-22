# Find pybind11 package, which is also going to find the python libraries and executable
find_package(pybind11 REQUIRED)

if(NOT pybind11_FOUND)
    set(BUILD_PYTHON OFF)
    message(WARNING "pybind11 was not found. "
        "Cannot build python wrappers and packages for Optima.")
endif()
    
# Check if python was found
if(NOT PYTHONINTERP_FOUND)
    set(BUILD_PYTHON OFF)
    message(WARNING "Python executable was not found. "
        "Cannot build python wrappers and packages for Reaktoro.")
endif()

# Build the Python wrappers of Reaktoro library
add_subdirectory(PyReaktoro)

# Build the Python package for Reaktoro
add_subdirectory(reaktoro)
