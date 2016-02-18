# Update the output path of the executable binaries produced by the project
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/tests/ReaktoroInterpreter)

# Recursively collect all header files from the current directory
file(GLOB_RECURSE HEADER_FILES *.hpp)

# Recursively collect all source files from the current directory
file(GLOB_RECURSE SOURCE_FILES *.cpp)

# Create an executable for each source file within the current directory
foreach(CPPFILE ${SOURCE_FILES})
  get_filename_component(filename ${CPPFILE} NAME_WE)
  add_executable(${filename} ${CPPFILE})
  target_link_libraries(${filename} ReaktoroShared)
  add_test(NAME ${filename} COMMAND ${filename})
endforeach()

