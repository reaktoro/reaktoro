# Append the path to the installed dependency libraries to CMAKE_PREFIX_PATH
list(APPEND CMAKE_PREFIX_PATH
    ${REAKTORO_DEPS_PATH}/public
    ${REAKTORO_DEPS_PATH}/private)

if(REAKTORO_USE_OPENLIBM)
    find_package(openlibm REQUIRED)
endif()

find_package(autodiff REQUIRED)
find_package(Optima REQUIRED)
find_package(ThermoFun REQUIRED)
find_package(nlohmann_json 3.4.0 REQUIRED)
find_package(yaml-cpp 0.6.3 REQUIRED)

if(ThermoFun_FOUND)
    add_compile_definitions(REAKTORO_USING_THERMOFUN)
endif()

if(REAKTORO_BUILD_TESTS)
    # Find catch2, which is used as the testing framework for Reaktoro
    find_package(Catch2 REQUIRED)
    if(NOT Catch2_FOUND)
        message(WARNING "Could not find Catch2 - the C++ tests of Reaktoro will not be built!")
        message(WARNING "Setting REAKTORO_BUILD_TESTS to OFF.")
        set(REAKTORO_BUILD_TESTS OFF)
    else()
        message(STATUS "Found Catch2 v${Catch2_VERSION}")
    endif()
endif()

if(REAKTORO_BUILD_PYTHON)
    # Find pybind11 library (used to generate the python bindings for Reaktoro)
    find_package(pybind11 REQUIRED)
    if(NOT pybind11_FOUND)
        message(WARNING "Could not find pybind11 - the Python module reaktoro will not be built!")
        message(WARNING "Setting REAKTORO_BUILD_PYTHON to OFF.")
        set(REAKTORO_BUILD_PYTHON OFF)
    endif()
endif()
