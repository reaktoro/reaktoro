# This cmake module determines if a conda environment is active.
# If so, CMAKE_INSTALL_PREFIX is set to the environmental variable CONDA_PREFIX.
# This automatic behavior can be overriden by manually specifying a different CMAKE_INSTALL_PREFIX.

# Check if user didn't override CMAKE_INSTALL_PREFIX and a conda environment is active
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT and DEFINED ENV{CONDA_PREFIX})
    message(STATUS "Conda environment recognized in $ENV{CONDA_PREFIX}")
    set(CMAKE_INSTALL_PREFIX $ENV{CONDA_PREFIX})
endif()
