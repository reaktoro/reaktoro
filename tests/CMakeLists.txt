set(REAKTORO_LIBRARY_DIR ${CMAKE_BINARY_DIR}/REAKTORO)  # the directory where the libraries are compiled
set(REAKTORO_PYTHON_PACKAGE_DIR ${CMAKE_BINARY_DIR}/python/package/build/lib)  # the directory where the python package is built

file(TO_NATIVE_PATH ${REAKTORO_LIBRARY_DIR} REAKTORO_LIBRARY_DIR)  # correct slashes if in Windows
file(TO_NATIVE_PATH ${REAKTORO_PYTHON_PACKAGE_DIR} REAKTORO_PYTHON_PACKAGE_DIR)  # correct slashes if in Windows

# Add target tests that performs all C++ and Python tests
if(WIN32)
    add_custom_target(tests
        COMMENT "Running C++ tests..."
        COMMAND ${CMAKE_COMMAND} -E env
            "PATH=${REAKTORO_DEPS_PATH}\\public\\bin;$ENV{PATH}"
                $<TARGET_FILE:reaktoro-cpptests>
        # COMMENT "Ensure python package REAKTORO exists in ${REAKTORO_PYTHON_PACKAGE_DIR}"
        # COMMENT "Ensure REAKTORO C++ library exists in ${REAKTORO_LIBRARY_DIR}\\$<CONFIG>"
        # COMMENT "Starting parallel execution of pytest..."
        # COMMAND ${CMAKE_COMMAND} -E env
        #     "PATH=${REAKTORO_LIBRARY_DIR}\\$<CONFIG>;$ENV{PATH}"
        #     "PYTHONPATH=${REAKTORO_PYTHON_PACKAGE_DIR}"
        #         pytest ${CMAKE_CURRENT_SOURCE_DIR} -n auto -x
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
else() # UNIX
    add_custom_target(tests
        COMMENT "Running C++ tests..."
        COMMAND $<TARGET_FILE:reaktoro-cpptests>
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endif()
