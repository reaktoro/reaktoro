# Ignore all warning from the dependencies
add_definitions(-w)

# Build the pugixml library
add_subdirectory(pugixml)

# Build the cvode library
add_subdirectory(cvode)

# Build the PHREEQC library
if(BUILD_PHREEQC)
    add_subdirectory(phreeqc)
endif(BUILD_PHREEQC)

# Build the GEMS library
if(BUILD_GEMS)
    add_subdirectory(gems)
endif(BUILD_GEMS)
