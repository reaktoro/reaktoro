# This cmake module determines if a conda environment is active.
# Note: CMAKE_INSTALL_PREFIX is set to the environment variable CONDA_PREFIX if a conda environment is found active.
# This automatic behavior can be overriden by manually specifying a different CMAKE_INSTALL_PREFIX.

# Check if a conda environment is active
if(DEFINED ENV{CONDA_PREFIX})
    # Show that conda env has been recognized
    message(STATUS "Conda environment recognized in $ENV{CONDA_PREFIX}")

    # Set CMAKE_INSTALL_PREFIX to CONDA_PREFIX if not specified by the user
    if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
        set(CMAKE_INSTALL_PREFIX $ENV{CONDA_PREFIX})
    endif()

    # Set BOOST_ROOT to CONDA_PREFIX to ensure conda's boost version is used
    set(BOOST_ROOT $ENV{CONDA_PREFIX})
endif()
