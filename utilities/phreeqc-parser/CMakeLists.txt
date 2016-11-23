# Require a certain version of cmake
cmake_minimum_required(VERSION 2.8)

file(GLOB CPPFILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

foreach(CPPFILE ${CPPFILES})
    get_filename_component(CPPNAME ${CPPFILE} NAME_WE)
    add_executable(${CPPNAME} ${CPPFILE})
    target_link_libraries(${CPPNAME} ReaktoroShared)
endforeach()
