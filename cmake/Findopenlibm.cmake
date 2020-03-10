# Finds the openlibm library
#
# The import target 'openlibm::openlibm' is created.
# Note that only the static library is supported for now.

find_path(OPENLIBM_INCLUDE_DIR openlibm)
if(WIN32)
    find_library(OPENLIBM_LIBRARIES libopenlibm_static)
else()
    find_library(OPENLIBM_LIBRARIES libopenlibm.a)
endif()

if(NOT TARGET openlibm::openlibm)
    add_library(openlibm::openlibm INTERFACE IMPORTED)
    if(WIN32)
        set_target_properties(openlibm::openlibm PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${OPENLIBM_INCLUDE_DIR}"
            INTERFACE_LINK_LIBRARIES "${OPENLIBM_LIBRARIES}"
        )
    elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        set_target_properties(openlibm::openlibm PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${OPENLIBM_INCLUDE_DIR}"
            INTERFACE_LINK_LIBRARIES "-Wl,--whole-archive ${OPENLIBM_LIBRARIES} -Wl,--no-whole-archive"
        )
    else()
        set_target_properties(openlibm::openlibm PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${OPENLIBM_INCLUDE_DIR}"
            INTERFACE_LINK_LIBRARIES "-Wl,-force_load,${OPENLIBM_LIBRARIES}"
        )
    endif()
endif()

function(configure_target_to_use_openlibm target_name)
    message(STATUS "Configuring to use openlibm: ${target_name}")
    target_link_libraries(
        ${target_name}
        PRIVATE
            openlibm::openlibm
    )
    if(WIN32)
        target_link_options(
            ${target_name}
            PRIVATE
                /FORCE:MULTIPLE
        )
    elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        target_compile_options(
            ${target_name}
            PRIVATE
                -fno-builtin
        )
        target_link_options(
            ${target_name}
            PRIVATE
                -Wl,-Bsymbolic-functions
        )
    else()
        target_compile_options(
            ${target_name}
            PRIVATE
                -fno-builtin
        )
    endif()

    # Remove transiency of openlibm::openlibm
    get_property(_target_interface_link_libraries TARGET ${target_name} PROPERTY INTERFACE_LINK_LIBRARIES)
    list(REMOVE_ITEM _target_interface_link_libraries "$<LINK_ONLY:openlibm::openlibm>")
    set_property(TARGET ${target_name} PROPERTY INTERFACE_LINK_LIBRARIES ${_target_interface_link_libraries})
endfunction()
