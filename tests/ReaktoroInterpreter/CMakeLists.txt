# Update the output path of the executable binaries produced by the project
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/tests/ReaktoroInterpreter)

# Recursively collect all header files from the current directory
file(GLOB_RECURSE HEADER_FILES *.hpp)

# Recursively collect all source files from the current directory
file(GLOB_RECURSE SOURCE_FILES *.cpp)

#foreach(CPPFILE ${CPPFILES})
#    get_filename_component(CPPNAME ${CPPFILE} NAME_WE)
#    add_executable(${CPPNAME} ${CPPFILE})
#    target_link_libraries(${CPPNAME} ReaktoroShared)
#endforeach()

# Create an executable for each source file within the current directory
foreach(sourcefile ${SOURCE_FILES})
  get_filename_component(filename ${sourcefile} NAME_WE)
  add_executable(${filename} ${sourcefile})
  target_link_libraries(${filename} ReaktoroInterpreter)
  add_test(NAME ${filename} COMMAND ${filename})
endforeach()

#add_executable(ReaktoroTester ${SOURCE_FILES})
#target_link_libraries(ReaktoroTester ${CMAKE_PROJECT_NAME} ${THIRD_PARTY_LIBRARIES})
#add_test(NAME ReaktoroTester COMMAND ReaktoroTester)