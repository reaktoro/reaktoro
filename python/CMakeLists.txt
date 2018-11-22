# Find pybind11 among the installed external dependencies
find_package(pybind11 REQUIRED
    NO_DEFAULT_PATH # prevent the use of pybind11 installed in the system
    PATHS ${REAKTORO_DEPENDENCIES_PATH})
    
# Build the Python wrappers of Reaktoro library
add_subdirectory(PyReaktoro)

# Build the Python package for Reaktoro
add_subdirectory(reaktoro)
