# Copy resources dir to the binary build dir along with the executable demos
file(COPY resources DESTINATION ${CMAKE_BINARY_DIR}/bin)

# Copy json dir to the binary build dir along with the executable demos
file(COPY json DESTINATION ${CMAKE_BINARY_DIR}/bin)

# Build the python and cpp demos
add_subdirectory(python)
add_subdirectory(cpp)
add_subdirectory(cpp-aux-tests)
