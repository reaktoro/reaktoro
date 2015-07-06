# Ignore all warning from the dependencies
if(CMAKE_COMPILER_IS_GNUCXX)
  add_definitions("-w")
endif()

# Build the pugixml library
add_subdirectory(pugixml)

# Build the tsolmod library
add_subdirectory(tsolmod)

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
