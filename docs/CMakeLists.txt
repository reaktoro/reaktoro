#------------------------------------------------------------------------------
# Version Requirement
#------------------------------------------------------------------------------
cmake_minimum_required(VERSION 2.8)

#------------------------------------------------------------------------------
# Project Details
#------------------------------------------------------------------------------
project(docs)

# Copy and configure the Doxyfile
configure_file(Doxyfile ${CMAKE_BINARY_DIR}/docs/Doxyfile @ONLY)

# Add target "docs" for building docs
add_custom_target(docs
  COMMAND doxygen Doxyfile
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
