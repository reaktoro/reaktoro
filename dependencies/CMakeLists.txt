# Add the subdirectories to be built
add_subdirectory(Reaktor/pugixml)

if(BUILD_PHREEQC)
    add_subdirectory(Reaktor/phreeqc)
endif(BUILD_PHREEQC)

if(BUILD_GEMS3K)
    add_subdirectory(Reaktor/gems3k)
endif(BUILD_GEMS3K)
