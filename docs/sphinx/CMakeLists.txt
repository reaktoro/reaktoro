add_custom_target(sphinx
  COMMAND sphinx-build
    -W                           # warning as errors
    -b html                      # build html
    ${CMAKE_CURRENT_SOURCE_DIR}  # source directory
    ${CMAKE_CURRENT_BINARY_DIR}  # output directory
)
