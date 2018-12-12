add_custom_target(sphinx
  COMMAND sphinx-build
    -b html                      # build html
    ${CMAKE_CURRENT_SOURCE_DIR}  # source directory
    ${CMAKE_CURRENT_BINARY_DIR}  # output directory
)
