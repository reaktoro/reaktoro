# Configure local pytest.ini file to permit pytest execution without path specification
configure_file(pytest.ini.in ${CMAKE_BINARY_DIR}/pytest.ini)

# Create cache variable with the number of parallel jobs to be used when executing pytest
set(REAKTORO_PYTEST_PARALLEL_JOBS 1 CACHE STRING "The number of parallel jobs when executing pytest.")

# Create auxiliary variable with parallel specification for pytest command
if(REAKTORO_PYTEST_PARALLEL_JOBS MATCHES 1)
    set(REAKTORO_PYTEST_PARALLEL_OPTION "")  # avoid passing "-n 1" to pytest to avoid overhead in creating a single parallel job
else()
    set(REAKTORO_PYTEST_PARALLEL_OPTION "-n ${REAKTORO_PYTEST_PARALLEL_JOBS}")
endif()

# Create target `tests-cpp` to execute C++ tests
add_custom_target(tests-cpp
    DEPENDS reaktoro-cpptests
    COMMENT "Running C++ tests..."
    COMMAND ${CMAKE_COMMAND} -E env
        "PATH=${REAKTORO_PATH}"
            $<TARGET_FILE:reaktoro-cpptests>
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR})

# Create target `tests-py` to execute Python tests
add_custom_target(tests-py
    DEPENDS reaktoro4py
    COMMENT "Running Python tests..."
    COMMENT "Starting parallel execution of pytest..."
    COMMAND ${CMAKE_COMMAND} -E env
        "PATH=${REAKTORO_PATH}"
        "PYTHONPATH=${REAKTORO_PYTHONPATH}"
            pytest -x ${REAKTORO_PYTEST_PARALLEL_OPTION}
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR})

# Combine both targets `tests-cpp` and `tests-py` to create target `tests`
add_custom_target(tests DEPENDS tests-py tests-cpp)
