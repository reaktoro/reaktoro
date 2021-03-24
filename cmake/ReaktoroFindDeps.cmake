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

# Find pybind11 library (if needed)
if(REAKTORO_BUILD_PYTHON)
    find_package(pybind11 REQUIRED)
    if(NOT pybind11_FOUND)
        message(WARNING "Could not find pybind11 - the Python module `reaktoro` will not be built.")
        set(REAKTORO_BUILD_PYTHON OFF)
    else()
        message(STATUS "Found pybind11 v${pybind11_VERSION}: ${pybind11_INCLUDE_DIRS}")
    endif()
endif()
