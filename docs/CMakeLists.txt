# Copy and configure the doxygen configuration files to binary dir
configure_file(doxyfile-cpp      ${CMAKE_CURRENT_BINARY_DIR}/doxyfile-cpp @ONLY)
configure_file(doxyfile-markdown ${CMAKE_CURRENT_BINARY_DIR}/doxyfile-markdown @ONLY)

# Add target "docs" for building docs
add_custom_target(docs
  COMMAND cmake -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/html/cpp
  COMMAND doxygen doxyfile-cpp
  COMMAND doxygen doxyfile-markdown
  COMMAND cmake -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/css ${CMAKE_CURRENT_BINARY_DIR}/html/css
  COMMAND cmake -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/fonts ${CMAKE_CURRENT_BINARY_DIR}/html/fonts
  COMMAND cmake -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/js ${CMAKE_CURRENT_BINARY_DIR}/html/js
  COMMAND cmake -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/img ${CMAKE_CURRENT_BINARY_DIR}/html/img
  COMMAND cmake -E copy_directory ${CMAKE_SOURCE_DIR}/databases ${CMAKE_CURRENT_BINARY_DIR}/html/databases
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")

# Add target "docs-cpp" for building only C++ documentation files
add_custom_target(docs-cpp
  COMMAND cmake -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/html/cpp
  COMMAND doxygen doxyfile-cpp
  COMMAND cmake -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/css ${CMAKE_CURRENT_BINARY_DIR}/html/css
  COMMAND cmake -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/fonts ${CMAKE_CURRENT_BINARY_DIR}/html/fonts
  COMMAND cmake -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/js ${CMAKE_CURRENT_BINARY_DIR}/html/js
  COMMAND cmake -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/img ${CMAKE_CURRENT_BINARY_DIR}/html/img
  COMMAND cmake -E copy_directory ${CMAKE_SOURCE_DIR}/databases ${CMAKE_CURRENT_BINARY_DIR}/html/databases
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")

# Add target "docs-markdown" for building only markdown documentation files
add_custom_target(docs-markdown
  COMMAND doxygen doxyfile-markdown
  COMMAND cmake -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/css ${CMAKE_CURRENT_BINARY_DIR}/html/css
  COMMAND cmake -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/fonts ${CMAKE_CURRENT_BINARY_DIR}/html/fonts
  COMMAND cmake -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/js ${CMAKE_CURRENT_BINARY_DIR}/html/js
  COMMAND cmake -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/img ${CMAKE_CURRENT_BINARY_DIR}/html/img
  COMMAND cmake -E copy_directory ${CMAKE_SOURCE_DIR}/databases ${CMAKE_CURRENT_BINARY_DIR}/html/databases
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")

# Add target "docs-upload" for uploading the built documentation to reaktoro.org
add_custom_target(docs-upload
  COMMAND rsync -avz --delete ${CMAKE_CURRENT_BINARY_DIR}/html/ reaktoro:public_html
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")

