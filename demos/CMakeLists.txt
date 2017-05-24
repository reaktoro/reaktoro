# Copy resources dir for the demos to the binary build dir along with the demos
file(COPY resources DESTINATION .)

# Build the cpp demos
add_subdirectory(cpp)