if(REAKTORO_ENABLE_OPENLIBM)
    find_package(openlibm REQUIRED)
endif()

set(REAKTORO_USE_autodiff        "" CACHE PATH "Specify this option in case a specific autodiff library should be used.")
set(REAKTORO_USE_Catch2          "" CACHE PATH "Specify this option in case a specific Catch2 library should be used.")
set(REAKTORO_USE_Eigen3          "" CACHE PATH "Specify this option in case a specific Eigen library should be used.")
set(REAKTORO_USE_nlohmann_json   "" CACHE PATH "Specify this option in case a specific nlohmann-json library should be used.")
set(REAKTORO_USE_Optima          "" CACHE PATH "Specify this option in case a specific Optima library should be used.")
set(REAKTORO_USE_phreeqc4rkt     "" CACHE PATH "Specify this option in case a specific Phreeqc library should be used.")
set(REAKTORO_USE_pybind11        "" CACHE PATH "Specify this option in case a specific pybind11 library should be used.")
set(REAKTORO_USE_ThermoFun       "" CACHE PATH "Specify this option in case a specific ThermoFun library should be used.")
set(REAKTORO_USE_yaml-cpp        "" CACHE PATH "Specify this option in case a specific yaml-cpp library should be used.")
set(REAKTORO_USE_tsl-ordered-map "" CACHE PATH "Specify this option in case a specific tsl-ordered-map library should be used.")

function(ReaktoroFindPackage package)
    # ARGV0: package name (required)
    # ARGV1: optional package version (optional)
    if(DEFINED REAKTORO_USE_${ARGV0} AND NOT REAKTORO_USE_${ARGV0} STREQUAL "")
        find_package(${ARGV0} ${ARGV1} REQUIRED EXACT QUIET PATHS ${REAKTORO_USE_${ARGV0}} NO_DEFAULT_PATH)
    else()
        find_package(${ARGV0} ${ARGV1} REQUIRED QUIET)
    endif()
    if(${ARGV0}_FOUND)
        message(STATUS "Found ${ARGV0}: ${${ARGV0}_DIR} (found version \"${${ARGV0}_VERSION}\")")
    endif()
    set(${ARGV0}_FOUND ${${ARGV0}_FOUND} PARENT_SCOPE)  # export to parent (e.g., Eigen3_FOUND)
    set(${ARGV0}_DIR ${${ARGV0}_DIR} PARENT_SCOPE)  # export to parent (e.g., Eigen3_DIR)
endfunction()

ReaktoroFindPackage(autodiff 0.6.12)
ReaktoroFindPackage(Eigen3 3.3.90)
ReaktoroFindPackage(nlohmann_json 3.6.1)
<<<<<<< HEAD
ReaktoroFindPackage(Optima 0.2.11)
=======
ReaktoroFindPackage(Optima 0.2.8)
>>>>>>> fa1dca8b4 (Depends on Optima v0.2.8 now)
ReaktoroFindPackage(phreeqc4rkt 3.6.2.1)
ReaktoroFindPackage(tabulate 1.4.0)
ReaktoroFindPackage(ThermoFun 0.3.8)
ReaktoroFindPackage(yaml-cpp 0.6.3)
ReaktoroFindPackage(tsl-ordered-map 1.0.0)

if(REAKTORO_BUILD_TESTS)
    ReaktoroFindPackage(Catch2 2.6.2)
    if(NOT Catch2_FOUND)
        message(WARNING "Could not find Catch2. The C++ tests of Reaktoro will not be built!")
        set(REAKTORO_BUILD_TESTS OFF)
    endif()
endif()

if(REAKTORO_BUILD_PYTHON)
    ReaktoroFindPackage(pybind11)
    find_program(PYBIND11_STUBGEN pybind11-stubgen)
    if(NOT pybind11_FOUND)
        message(WARNING "Could not find pybind11. The Python package reaktoro will not be built!")
        set(REAKTORO_BUILD_PYTHON OFF)
    endif()
    if(NOT PYBIND11_STUBGEN)
        message(WARNING "Could not find pybind11-stubgen (available via pip or conda). There will be no stubs generated for the python package reaktoro.")
    endif()
endif()
