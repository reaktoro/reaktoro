# Copy and configure the doxygen configuration files to binary dir
configure_file(doxyfile ${CMAKE_CURRENT_BINARY_DIR}/doxyfile @ONLY)

# Add target "docs" for building docs
add_custom_target(docs
  COMMAND cmake -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/html/cpp
  COMMAND doxygen doxyfile
  COMMAND cmake -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/css ${CMAKE_CURRENT_BINARY_DIR}/html/css
  COMMAND cmake -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/fonts ${CMAKE_CURRENT_BINARY_DIR}/html/fonts
  COMMAND cmake -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/js ${CMAKE_CURRENT_BINARY_DIR}/html/js
  COMMAND cmake -E copy_directory ${CMAKE_SOURCE_DIR}/docs/img ${CMAKE_CURRENT_BINARY_DIR}/html/img
  COMMAND cmake -E copy_directory ${CMAKE_SOURCE_DIR}/databases ${CMAKE_CURRENT_BINARY_DIR}/html/databases
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")

# Add target "docs-upload" for uploading the built documentation to reaktoro.org
add_custom_target(docs-upload
  COMMAND rsync -avz --delete ${CMAKE_CURRENT_BINARY_DIR}/html/ reaktoro:public_html
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
