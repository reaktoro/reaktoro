#------------------------------------------------------------------------------
# Version Requirement
#------------------------------------------------------------------------------
cmake_minimum_required(VERSION 2.8)

#------------------------------------------------------------------------------
# Header and Source Files
#------------------------------------------------------------------------------
# Recursively collect all header files from the current directory
file(GLOB_RECURSE HEADER_FILES *.hpp)

# Recursively collect all source files from the current directory
file(GLOB_RECURSE SOURCE_FILES *.cpp)

#------------------------------------------------------------------------------
# Build and Install Steps
#------------------------------------------------------------------------------
# add_library(reaktor OBJECT ${HEADER_FILES} ${SOURCE_FILES})

# if(BUILD_STATIC_LIBS)
#   # The name of the static library
#   set(REAKTOR_STATIC_LIB ${PROJECT_NAME}${SUFFIX_STATIC_LIBS})

#   # Create the static library using the collected source files
#   add_library(${REAKTOR_STATIC_LIB} STATIC $<TARGET_OBJECTS:reaktor>)

#   # Link the static library with the required third-party libraries
#   target_link_libraries(${REAKTOR_STATIC_LIB} ${THIRD_PARTY_LIBRARIES})

#   # Apply the cotire function to the library target
#   cotire(${REAKTOR_STATIC_LIB})

#   # Install the target libraries
#   install(TARGETS ${REAKTOR_STATIC_LIB} DESTINATION lib)
# endif()

# if(BUILD_SHARED_LIBS)
#   # The name of the shared library
#   set(REAKTOR_SHARED_LIB ${PROJECT_NAME}${SUFFIX_SHARED_LIBS})

#   # Create the shared library using the collected source files
#   add_library(${REAKTOR_SHARED_LIB} SHARED $<TARGET_OBJECTS:reaktor>)

#   # Link the shared library with the required third-party libraries
#   target_link_libraries(${REAKTOR_SHARED_LIB} ${THIRD_PARTY_LIBRARIES})

#   # Apply the cotire function to the library target
#   cotire(${REAKTOR_SHARED_LIB})

#   # Install the target shared library
#   install(TARGETS ${REAKTOR_SHARED_LIB} DESTINATION lib)
# endif()

# The name of the shared library
set(REAKTOR_SHARED_LIB ${PROJECT_NAME}${SUFFIX_SHARED_LIBS})

# Create the shared library using the collected source files
add_library(${REAKTOR_SHARED_LIB} SHARED ${HEADER_FILES} ${SOURCE_FILES})

# Link the shared library with the required third-party libraries
target_link_libraries(${REAKTOR_SHARED_LIB} ${THIRD_PARTY_LIBRARIES})

# Apply the cotire function to the library target
cotire(${REAKTOR_SHARED_LIB})

# Install the target shared library
install(TARGETS ${REAKTOR_SHARED_LIB} DESTINATION lib)

# Install the header files preserving the directory hierarchy
install(DIRECTORY ${REAKTOR_SOURCE_DIR} DESTINATION include FILES_MATCHING PATTERN "*.hpp" PATTERN "*.hxx")