# Build doxystrap documentation
add_subdirectory(doxystrap)

# Build sphinx documentation (readthedocs)
add_subdirectory(sphinx)

# Create a target to build website
add_custom_target(website
  COMMAND ${CMAKE_COMMAND} --build . --target sphinx
  COMMAND ${CMAKE_COMMAND} --build . --target doxystrap
  COMMAND ${CMAKE_COMMAND} -E copy_directory docs/sphinx/html docs/website/
  COMMAND ${CMAKE_COMMAND} -E copy_directory docs/doxystrap/html docs/website/cpp
  WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
)

# Add target "docs-upload" for uploading the built documentation to reaktoro.org
add_custom_target(website-upload
  COMMAND rsync -avz --delete ${CMAKE_CURRENT_BINARY_DIR}/website/ reaktoro:public_html
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
