# Recursively collect all hpp files from the current directory
file(GLOB_RECURSE HPP_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.hpp)

# Recursively collect all cpp files from the current directory
file(GLOB_RECURSE CPP_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

# Recursively collect all .test.cxx files from the current directory
file(GLOB_RECURSE CXX_FILES_TEST RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.test.cxx)

# Recursively collect all .py.cxx files from the current directory
file(GLOB_RECURSE CXX_FILES_PY RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.py.cxx)

# Enable automatic creation of a module definition (.def) file for a SHARED library on Windows.
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS TRUE)

# Create a library target for Reaktoro
add_library(Reaktoro ${HPP_FILES} ${CPP_FILES})

# Add an alias Reaktoro::Reaktoro to the target library Reaktoro
add_library(Reaktoro::Reaktoro ALIAS Reaktoro)

# Specify the external dependencies of Reaktoro to ensure proper build sequence
add_dependencies(Reaktoro Reaktoro::Embedded)

# Add the include paths to Reaktoro library target
target_include_directories(Reaktoro
    PUBLIC $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    PRIVATE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>)

# Link Reaktoro library against external dependencies
target_link_libraries(Reaktoro
    PRIVATE Reaktoro::Embedded
    PRIVATE tabulate::tabulate
    PUBLIC Eigen3::Eigen
    PUBLIC autodiff::autodiff
    PUBLIC Eigen3::Eigen
    PUBLIC nlohmann_json::nlohmann_json
    PUBLIC Optima::Optima
    PUBLIC phreeqc4rkt::phreeqc4rkt
    PUBLIC ThermoFun::ThermoFun
    PUBLIC tsl::ordered_map
    PUBLIC yaml-cpp
)

# Enable implicit conversion of autodiff::real to double
target_compile_definitions(Reaktoro PUBLIC AUTODIFF_ENABLE_IMPLICIT_CONVERSION_REAL=1)

if(REAKTORO_ENABLE_OPENLIBM)
    configure_target_to_use_openlibm(Reaktoro)
    target_compile_definitions(Reaktoro PUBLIC REAKTORO_ENABLE_OPENLIBM=1)
endif()

# Set compilation features to be propagated to dependent codes.
target_compile_features(Reaktoro PUBLIC cxx_std_17)

# Set some target properties
set_target_properties(Reaktoro PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    CXX_EXTENSIONS OFF)

#==============================================================================
#
# COMPILATION CONFIGURATION FOR C++ TESTS
#
#==============================================================================

if(REAKTORO_BUILD_TESTS)

    # Find catch2, which is used as the testing framework for Reaktoro
    find_package(Catch2 REQUIRED)

    # The global path to the databases directory (used below as a compile definition for the C++ tests)
    set(REAKTORO_EMBEDDED_DIR  ${PROJECT_SOURCE_DIR}/embedded)
    set(REAKTORO_DATABASES_DIR ${PROJECT_SOURCE_DIR}/embedded/databases)
    set(REAKTORO_PARAMS_DIR    ${PROJECT_SOURCE_DIR}/embedded/params)

    # Create a test executable target for Reaktoro
    add_executable(reaktoro-cpptests ${CXX_FILES_TEST})
    target_link_libraries(reaktoro-cpptests Reaktoro Catch2::Catch2)
    target_include_directories(reaktoro-cpptests PUBLIC ${PROJECT_SOURCE_DIR})
    target_compile_definitions(reaktoro-cpptests
        PRIVATE $<$<CXX_COMPILER_ID:MSVC>:EIGEN_STRONG_INLINE=inline> # This is needed (`EIGEN_STRONG_INLINE=inline` instead of `EIGEN_STRONG_INLINE=__forceinline`) with MSVC, otherwise some tests (ChemicalProps.test.cxx) takes too long to compile.
        PRIVATE REAKTORO_EMBEDDED_DIR="${REAKTORO_EMBEDDED_DIR}"      # This permits the C++ tests to easily load embedded resources via global addresses so that the tests can be executed from anywhere without errors.
        PRIVATE REAKTORO_DATABASES_DIR="${REAKTORO_DATABASES_DIR}"    # This permits the C++ tests to easily load embedded databases via global addresses so that the tests can be executed from anywhere without errors.
        PRIVATE REAKTORO_PARAMS_DIR="${REAKTORO_PARAMS_DIR}"          # This permits the C++ tests to easily load embedded model parameters via global addresses so that the tests can be executed from anywhere without errors.
    )

endif()

#==============================================================================
#
# COMPILATION CONFIGURATION FOR PYTHON BINDINGS
#
#==============================================================================

if(REAKTORO_BUILD_PYTHON)

    # Find pybind11, needed to compile the python bindings
    find_package(pybind11 REQUIRED)

    # Compile Reaktoro *.py.cxx files into a Python module
    pybind11_add_module(reaktoro4py ${CXX_FILES_PY})

    # Add the root directory of the project to the include list
    target_include_directories(reaktoro4py PRIVATE ${PROJECT_SOURCE_DIR})

    # Link reaktoro4py against Reaktoro C++ library
    target_link_libraries(reaktoro4py LINK_PRIVATE Reaktoro::Reaktoro)

    # This is needed to avoid compilation issues when using clang-9 and c++17
    if(${CMAKE_CXX_COMPILER_ID} STREQUAL Clang AND ${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS 10.0.0)
        target_compile_options(reaktoro4py PRIVATE -fsized-deallocation)
    endif()

endif()

#==============================================================================
#
# INSTALL CONFIGURATION
#
#==============================================================================

# Install Reaktoro C++ library
install(TARGETS Reaktoro
    EXPORT ReaktoroTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT libraries
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT libraries
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT libraries)

# Install Reaktoro debug symbols for MSVC builds
if(MSVC)
    install(FILES $<TARGET_PDB_FILE:Reaktoro>
        DESTINATION ${CMAKE_INSTALL_BINDIR}
        COMPONENT libraries
        OPTIONAL)
endif()

# Install Reaktoro header files
install(DIRECTORY ${PROJECT_SOURCE_DIR}/Reaktoro
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} COMPONENT headers
    FILES_MATCHING
        PATTERN "*.hpp"
        PATTERN "__pycache__" EXCLUDE)
