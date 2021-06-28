# Set variable with the directory where the Reaktoro C++ library is located
set(REAKTORO_LIBRARY_DIR ${CMAKE_BINARY_DIR}/Reaktoro)

# Set variable with the directory where the python package reaktoro is built
set(REAKTORO_PYTHON_PACKAGE_DIR ${CMAKE_BINARY_DIR}/python/package/installed/$<CONFIG>/Lib/site-packages)

# Set variable with the directory where the python package autodiff is built
set(AUTODIFF_PYTHON_PACKAGE_DIR ${CMAKE_SOURCE_DIR}/deps/build/etc/src/AUTODIFF-build/python/package/autodiff/build/lib)

# Set variable with the directory where the python package optima is built
set(OPTIMA_PYTHON_PACKAGE_DIR ${CMAKE_SOURCE_DIR}/deps/build/etc/src/OPTIMA-build/python/package/build/lib)

# Correct the above paths if in Windows (e.g., slashes / to backslashes \)
file(TO_NATIVE_PATH ${REAKTORO_LIBRARY_DIR}        REAKTORO_LIBRARY_DIR)
file(TO_NATIVE_PATH ${REAKTORO_PYTHON_PACKAGE_DIR} REAKTORO_PYTHON_PACKAGE_DIR)
file(TO_NATIVE_PATH ${AUTODIFF_PYTHON_PACKAGE_DIR} AUTODIFF_PYTHON_PACKAGE_DIR)
file(TO_NATIVE_PATH ${OPTIMA_PYTHON_PACKAGE_DIR}   OPTIMA_PYTHON_PACKAGE_DIR)

# Add target tests that performs all C++ and Python tests
if(WIN32)

    set(CUSTOM_PATH "${REAKTORO_LIBRARY_DIR}\\$<CONFIG>\;${REAKTORO_DEPS_PATH}\\public\\bin\;$ENV{PATH}")
    set(CUSTOM_PYTHONPATH "${REAKTORO_PYTHON_PACKAGE_DIR};${AUTODIFF_PYTHON_PACKAGE_DIR};${OPTIMA_PYTHON_PACKAGE_DIR};$ENV{PYTHONPATH}")

    add_custom_target(tests-cpp
        COMMENT "Running C++ tests..."
        COMMAND ${CMAKE_COMMAND} -E env
            "PATH=${CUSTOM_PATH}"
                $<TARGET_FILE:reaktoro-cpptests>
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )

    add_custom_target(tests-py
        COMMENT "Running Python tests..."
        COMMENT "Ensure python package reaktoro exists in ${REAKTORO_PYTHON_PACKAGE_DIR}"
        COMMENT "Ensure python package autodiff exists in ${AUTODIFF_PYTHON_PACKAGE_DIR}"
        COMMENT "Ensure python package optima exists in ${OPTIMA_PYTHON_PACKAGE_DIR}"
        COMMENT "Starting parallel execution of pytest..."
        COMMAND ${CMAKE_COMMAND} -E echo "CUSTOM_PATH = ${CUSTOM_PATH}"
        COMMAND ${CMAKE_COMMAND} -E echo "CUSTOM_PYTHONPATH = ${CUSTOM_PYTHONPATH}"
        COMMAND ${CMAKE_COMMAND} -E env
            "PATH=${CUSTOM_PATH}"
            "PYTHONPATH=${CUSTOM_PYTHONPATH}"
                pytest ${CMAKE_SOURCE_DIR}/Reaktoro -n auto -x
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )

elseif(UNIX)

    add_custom_target(tests-cpp
        COMMENT "Running C++ tests..."
        COMMAND $<TARGET_FILE:reaktoro-cpptests>
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )

    set(CUSTOM_PYTHONPATH "${REAKTORO_PYTHON_PACKAGE_DIR}:${AUTODIFF_PYTHON_PACKAGE_DIR}:${OPTIMA_PYTHON_PACKAGE_DIR}")

    add_custom_target(tests-py
        COMMENT "Running Python tests..."
        COMMENT "Ensure python package reaktoro exists in ${REAKTORO_PYTHON_PACKAGE_DIR}"
        COMMENT "Ensure python package autodiff exists in ${AUTODIFF_PYTHON_PACKAGE_DIR}"
        COMMENT "Ensure python package optima exists in ${OPTIMA_PYTHON_PACKAGE_DIR}"
        COMMENT "Starting parallel execution of pytest..."
        COMMAND ${CMAKE_COMMAND} -E env
            "PYTHONPATH=${CUSTOM_PYTHONPATH}"
                pytest ${CMAKE_SOURCE_DIR}/Reaktoro -n auto -x
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
else()

    message(FATAL_ERROR "Could not create cmake targets tests-cpp and tests-py. Expecting Unix or Windows systems.")

endif()

# Combine both targets `tests-cpp` and `tests-py` to create target `tests`
add_custom_target(tests DEPENDS tests-py tests-cpp)
