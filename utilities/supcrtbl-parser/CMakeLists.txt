add_custom_target(update-supcrtbl-database
    COMMENT "Updating supcrtbl json and yaml databases..."
    COMMAND ${CMAKE_COMMAND} -E echo "This operation expects databases/supcrtbl/spronsbl.dat has been updated!"
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/parse-supcrtbl-database
        -i     ${CMAKE_SOURCE_DIR}/databases/supcrtbl/spronsbl.dat
        -ojson ${CMAKE_SOURCE_DIR}/databases/supcrtbl/supcrtbl.json
        -oyaml ${CMAKE_SOURCE_DIR}/databases/supcrtbl/supcrtbl.yaml
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
