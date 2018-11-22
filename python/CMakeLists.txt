# Find pybind11 among the installed external dependencies
find_package(pybind11 REQUIRED)

# Build the Python wrappers of Reaktoro library
add_subdirectory(PyReaktoro)

# Build the Python package for Reaktoro
add_subdirectory(reaktoro)
