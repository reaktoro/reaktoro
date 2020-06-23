find_package(Boost REQUIRED)
find_package(nlohmann_json 3.6.1 REQUIRED)

if(REAKTORO_USE_OPENLIBM)
    find_package(openlibm REQUIRED)
endif()

# Find ThermoFun library
find_package(ThermoFun)

# Define the pre-processor macro REAKTORO_USING_THERMOFUN if found.
if(ThermoFun_FOUND)
    add_compile_definitions(REAKTORO_USING_THERMOFUN)
endif()

# Find pybind11 library (if needed)
if(REAKTORO_BUILD_PYTHON)
    set(PYBIND11_CPP_STANDARD -std=c++17)  # Ensure pybind11 really uses c++17
    find_package(pybind11 2.3.0 REQUIRED)
    if(NOT pybind11_FOUND)
        message(WARNING "Could not find pybind11 - the Python module `reaktoro` will not be built.")
        set(REAKTORO_BUILD_PYTHON OFF)
    else()
        message(STATUS "Found pybind11 v${pybind11_VERSION}: ${pybind11_INCLUDE_DIRS}")
    endif()
endif()
