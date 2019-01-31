# Collect all python files in the demos/python directory
file(GLOB_RECURSE PYFILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.py)

# Copy all collected python files to demos/python in the binary directory
file(COPY ${PYFILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Copy the databases directory to demos/python/databases in the binary directory
file(COPY ${PROJECT_SOURCE_DIR}/databases DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
