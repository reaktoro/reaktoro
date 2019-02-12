# Create a sphinx-build target to build the sphinx documentation website
add_custom_target(sphinx
  COMMAND sphinx-build
    -b html                           # build html
    ${CMAKE_CURRENT_SOURCE_DIR}       # source directory
    ${CMAKE_CURRENT_BINARY_DIR}/html  # output directory
    )

    # Create a sphinx-autobuild target to automatically rebuild the sphinx documentation website
    add_custom_target(sphinx-serve
    COMMAND sphinx-autobuild
    -b html                           # build html
    ${CMAKE_CURRENT_SOURCE_DIR}       # source directory
    ${CMAKE_CURRENT_BINARY_DIR}/html  # output directory
)
