# Recursively collect all header files from the current directory
file(GLOB_RECURSE HPP_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.hpp *.hxx)

# Recursively collect all source files from the current directory
file(GLOB_RECURSE CPP_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

# Recursively collect all test files from the current directory
file(GLOB_RECURSE TEST_CXX_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.test.cxx)

# Enable automatic creation of a module definition (.def) file for a SHARED library on Windows.
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS TRUE)

# Add the library path of the external dependencies to link directories
# Note: cmake v3.13 is needed to use target_link_directories command!
link_directories(${REAKTORO_DEPS_LIBRARY_PATH})

# Create a library target for Reaktoro
add_library(Reaktoro ${HPP_FILES} ${CPP_FILES})

# Add an alias Reaktoro::Reaktoro to the target library Reaktoro
add_library(Reaktoro::Reaktoro ALIAS Reaktoro)

# Specify the external dependencies of Reaktoro to ensure proper build sequence
add_dependencies(Reaktoro Reaktoro::Databases)

# Add the include paths to Reaktoro library target
target_include_directories(Reaktoro
    PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/deps/eigen3>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/Reaktoro/deps/eigen3>
    PRIVATE ${REAKTORO_INCLUDE_PATH_DEPS})

# Link Reaktoro library against external dependencies
target_link_libraries(Reaktoro
    PRIVATE pugixml
    PRIVATE sundials_cvode
    PRIVATE sundials_nvecserial
    PRIVATE phreeqc
    # PRIVATE gems
    PRIVATE nlohmann_json::nlohmann_json
    PRIVATE Reaktoro::Databases
    PUBLIC autodiff::autodiff
    PUBLIC Optima::Optima
    PUBLIC yaml-cpp
)

# Enable implicit conversion of autodiff::real to double
target_compile_definitions(Reaktoro PUBLIC AUTODIFF_ENABLE_IMPLICIT_CONVERSION_REAL=1)

if(REAKTORO_USE_OPENLIBM)
    configure_target_to_use_openlibm(Reaktoro)
    target_compile_definitions(Reaktoro PUBLIC REAKTORO_USE_OPENLIBM=1)
endif()

# Link Reaktoro library against ThermoFun if found
if(ThermoFun_FOUND)
    target_link_libraries(Reaktoro PUBLIC ThermoFun::ThermoFun)
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

    # Create a test executable target for Reaktoro
    add_executable(reaktoro-cpptests ${TEST_CXX_FILES})
    target_link_libraries(reaktoro-cpptests Reaktoro Catch2::Catch2)
    target_include_directories(reaktoro-cpptests PUBLIC ${PROJECT_SOURCE_DIR})
    target_compile_definitions(reaktoro-cpptests PRIVATE $<$<CXX_COMPILER_ID:MSVC>:EIGEN_STRONG_INLINE=inline>)  # use EIGEN_STRONG_INLINE=inline instead of EIGEN_STRONG_INLINE=__forceinline with MSVC, otherwise some tests (ChemicalProps.test.cxx) takes too long to compile

    # Add target tests that performs all C++ and Python tests
    add_custom_target(tests
        COMMENT "Running C++ tests..."
        COMMAND $<TARGET_FILE:reaktoro-cpptests>
        # COMMENT "Running Python tests..."  # TODO: Implement Python tests.
        # COMMAND $<TARGET_FILE:reaktoro-pytests>
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

endif()

#==============================================================================
#
# COMPILATION CONFIGURATION FOR PYTHON BINDINGS
#
#==============================================================================

if(REAKTORO_BUILD_PYTHON)

    # Include the Eigen include path for pybind11
    include_directories(${PROJECT_SOURCE_DIR}/Reaktoro/deps/eigen3)

    # Find pybind11, needed to compile the python bindings
    find_package(pybind11 REQUIRED)

    # Recursively collect all python binding files from the current directory
    file(GLOB_RECURSE PY_CXX_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.py.cxx)

    # Compile Reaktoro *.py.cxx files into a Python module
    pybind11_add_module(reaktoro4py ${PY_CXX_FILES})

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

# Install debug symbols
if(MSVC)
    install(
        FILES $<TARGET_PDB_FILE:Reaktoro>
        DESTINATION ${CMAKE_INSTALL_BINDIR}
        COMPONENT libraries
        OPTIONAL
    )
endif()

# Install Reaktoro header files
install(DIRECTORY ${PROJECT_SOURCE_DIR}/Reaktoro
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} COMPONENT headers
    FILES_MATCHING PATTERN "*.hpp" PATTERN "*.hxx")

# Install Reaktoro/deps header files
install(DIRECTORY ${PROJECT_SOURCE_DIR}/Reaktoro/deps
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Reaktoro COMPONENT headers)

# Install public external deps libraries
install(DIRECTORY ${PROJECT_SOURCE_DIR}/deps/build/install/public/
    DESTINATION ${CMAKE_INSTALL_PREFIX} COMPONENT libraries)
