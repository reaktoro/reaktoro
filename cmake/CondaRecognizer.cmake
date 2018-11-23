# This cmake module determines if a conda environment is active.
# If so, CMAKE_INSTALL_PREFIX is set to the environmental variable CONDA_PREFIX.
# This automatic behavior can prevented by passing the option -DCMAKE_CONDA_IGNORE=TRUE


# Keep the previous value of CMAKE_INSTALL_PREFIX before a possible change to CONDA_PREFIX
set(CMAKE_INSTALL_PREFIX_BEFORE_CONDA_PREFIX ${CMAKE_INSTALL_PREFIX})

# Check if a conda environment is active
if(DEFINED ENV{CONDA_PREFIX} AND NOT CMAKE_CONDA_IGNORE)
    message(STATUS "Conda environment recognized in $ENV{CONDA_PREFIX}")
    set(CMAKE_INSTALL_PREFIX $ENV{CONDA_PREFIX})
endif()

# If conda environment should be ignored, then set CMAKE_INSTALL_PREFIX to its previous value
if(CMAKE_CONDA_IGNORE)
    set(CMAKE_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX_BEFORE_CONDA_PREFIX})
endif()