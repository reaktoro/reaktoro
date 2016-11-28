# Copy and configure the files to binary dir
configure_file(Doxyfile ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)

# Add target "docs" for building docs
add_custom_target(docs
  COMMAND doxygen Doxyfile
  COMMAND cmake -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/css ${CMAKE_CURRENT_BINARY_DIR}/html/css
  COMMAND cmake -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/fonts ${CMAKE_CURRENT_BINARY_DIR}/html/fonts
  COMMAND cmake -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/js ${CMAKE_CURRENT_BINARY_DIR}/html/js
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
