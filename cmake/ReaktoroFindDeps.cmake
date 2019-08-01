find_package(Boost REQUIRED)
find_package(nlohmann_json 3.6.1 REQUIRED)

if(REAKTORO_BUILD_PYTHON)

    find_package(pybind11 2.3.0 REQUIRED)
    if(NOT pybind11_FOUND)
        message(WARNING "Could not find pybind11 - the Python module `reaktoro` will not be built.")
        set(REAKTORO_BUILD_PYTHON OFF)
    else()
        message(STATUS "Found pybind11 v${pybind11_VERSION}: ${pybind11_INCLUDE_DIRS}")
    endif()

endif()
