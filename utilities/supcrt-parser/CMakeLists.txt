add_custom_target(update-databases-supcrt
    COMMENT "Updating supcrt json and yaml databases..."
    COMMAND ${CMAKE_COMMAND} -E echo "This operation expects files linearized.supcrt**.dat have been updated!"
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/generate-supcrt-databases
        --filename ${CMAKE_CURRENT_SOURCE_DIR}/linearized.supcrt98.dat
        --outdir ${CMAKE_SOURCE_DIR}/databases/reaktoro
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/generate-supcrt-databases
        --filename ${CMAKE_CURRENT_SOURCE_DIR}/linearized.supcrt07.dat
        --outdir ${CMAKE_SOURCE_DIR}/databases/reaktoro
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/generate-supcrt-databases
        --filename ${CMAKE_CURRENT_SOURCE_DIR}/linearized.supcrt16.dat
        --outdir ${CMAKE_SOURCE_DIR}/databases/reaktoro
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
