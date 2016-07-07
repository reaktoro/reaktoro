# Copy and configure the files to binary dir
configure_file(Doxyfile ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)

# # Set bootstrap variables
# set(BOOTSTRAP_ZIP ${CMAKE_SOURCE_DIR}/docs/dependencies/bootstrap-3.3.6-dist.zip)
# set(BOOTSTRAP_FROM ${CMAKE_CURRENT_BINARY_DIR}/bootstrap-3.3.6-dist)
# set(BOOTSTRAP_TO ${CMAKE_CURRENT_BINARY_DIR}/html/css/)

# # Set font-awesome variables
# set(FONT_AWESOME_ZIP ${CMAKE_SOURCE_DIR}/docs/dependencies/font-awesome-4.6.3.zip)
# set(FONT_AWESOME_FROM ${CMAKE_CURRENT_BINARY_DIR}/font-awesome-4.6.3)
# set(FONT_AWESOME_TO ${CMAKE_CURRENT_BINARY_DIR}/html/fonts/)

# # Add target "docs" for building docs
# add_custom_target(docs
#   COMMAND cmake -E tar x ${BOOTSTRAP_ZIP}
#   COMMAND cmake -E tar x ${FONT_AWESOME_ZIP}
#   COMMAND cmake -E tar x ${DOXYSTRAP_ZIP}
#   COMMAND cmake -E copy_directory ${BOOTSTRAP_FROM} ${BOOTSTRAP_TO}
#   COMMAND cmake -E copy_directory ${FONT_AWESOME_FROM} ${FONT_AWESOME_TO}
#   COMMAND cmake -E remove_directory ${DOXYSTRAP_TO}
#   COMMAND cmake -E rename ${DOXYSTRAP_FROM} ${DOXYSTRAP_TO}
#   COMMAND doxygen Doxyfile
#   WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")

# Add target "docs" for building docs
add_custom_target(docs
  COMMAND doxygen Doxyfile
  COMMAND cmake -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/css ${CMAKE_CURRENT_BINARY_DIR}/html/css
  COMMAND cmake -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/fonts ${CMAKE_CURRENT_BINARY_DIR}/html/fonts
  COMMAND cmake -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/js ${CMAKE_CURRENT_BINARY_DIR}/html/js
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
