# Require a certain version of cmake
cmake_minimum_required(VERSION 2.8)

# Recursively collect all header files from the current directory
file(GLOB_RECURSE HEADER_FILES *.hpp)

# Recursively collect all source files from the current directory
file(GLOB_RECURSE SOURCE_FILES *.cpp)

# The name of the shared library
set(REAKTOR_SHARED_LIB ${PROJECT_NAME}${SUFFIX_SHARED_LIBS})

# Create the shared library using the collected source files
add_library(${REAKTOR_SHARED_LIB} SHARED ${HEADER_FILES} ${SOURCE_FILES})

# Link the shared library with the required third-party libraries
target_link_libraries(${REAKTOR_SHARED_LIB} ${THIRD_PARTY_LIBRARIES} gems3k phreeqc)

# Apply the cotire function to the library target
cotire(${REAKTOR_SHARED_LIB})

# Install the target shared library
install(TARGETS ${REAKTOR_SHARED_LIB} DESTINATION lib)

# Install the header files preserving the directory hierarchy
install(DIRECTORY ${REAKTOR_SOURCE_DIR} DESTINATION include FILES_MATCHING PATTERN "*.hpp" PATTERN "*.hxx")