# Ignore all warning from the dependencies
add_definitions(-w)

# Build the pugixml library
add_subdirectory(Reaktor/pugixml)

# Build the PHREEQC library
if(BUILD_PHREEQC)
    add_subdirectory(Reaktor/phreeqc)
endif(BUILD_PHREEQC)

# Build the GEMS library
if(BUILD_GEMS3K)
    add_subdirectory(Reaktor/gems3k)
endif(BUILD_GEMS3K)
