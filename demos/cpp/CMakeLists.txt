file(GLOB_RECURSE CPPFILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)
find_package(Boost COMPONENTS system filesystem REQUIRED)

foreach(CPPFILE ${CPPFILES})
    get_filename_component(CPPNAME ${CPPFILE} NAME_WE)
    add_executable(${CPPNAME} ${CPPFILE})
	target_link_libraries(${CPPNAME} Reaktoro::Reaktoro)
endforeach()
