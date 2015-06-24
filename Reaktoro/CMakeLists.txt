# Require a certain version of cmake
cmake_minimum_required(VERSION 2.8)

# Recursively collect all header files from the current directory
file(GLOB_RECURSE HEADER_FILES *.hpp)

# Recursively collect all source files from the current directory
file(GLOB_RECURSE SOURCE_FILES *.cpp)

# The name of the shared and static libraries
set(REAKTORO_SHARED_LIB ${PROJECT_NAME}${SUFFIX_SHARED_LIBS})
set(REAKTORO_STATIC_LIB ${PROJECT_NAME}${SUFFIX_STATIC_LIBS})

# Create a list of target objects
# Link the shared library with the pugixml library
set(REAKTORO_TARGET_OBJECTS ${REAKTORO_TARGET_OBJECTS} $<TARGET_OBJECTS:pugixml>)

# Link the shared library with the optional PHREEQC library
if(BUILD_PHREEQC)
    set(REAKTORO_TARGET_OBJECTS ${REAKTORO_TARGET_OBJECTS} $<TARGET_OBJECTS:phreeqc>)
endif()

# Link the shared library with the optional GEMS3K library
if(BUILD_GEMS)
    set(REAKTORO_TARGET_OBJECTS ${REAKTORO_TARGET_OBJECTS} $<TARGET_OBJECTS:gems>)
endif()

# Create the shared library using the collected source files
add_library(${REAKTORO_SHARED_LIB} SHARED ${HEADER_FILES} ${SOURCE_FILES} ${REAKTORO_TARGET_OBJECTS})

# Link the shared library with the CVODE library
target_link_libraries(${REAKTORO_SHARED_LIB} sundials_cvode_shared)
target_link_libraries(${REAKTORO_SHARED_LIB} sundials_nvecserial_shared)

# Install the target shared library
install(TARGETS ${REAKTORO_SHARED_LIB} DESTINATION lib)

# Install the header files preserving the directory hierarchy
install(DIRECTORY ${REAKTORO_SOURCE_DIR} DESTINATION include FILES_MATCHING PATTERN "*.hpp" PATTERN "*.hxx")
