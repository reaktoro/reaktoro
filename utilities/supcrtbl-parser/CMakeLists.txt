add_custom_target(update-databases-supcrtbl
    COMMENT "Updating supcrtbl json and yaml databases..."
    COMMAND ${CMAKE_COMMAND} -E echo "This operation expects spronsbl.dat has been updated!"
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/generate-supcrtbl-databases
        -i     ${CMAKE_CURRENT_SOURCE_DIR}/spronsbl.dat
        -ojson ${CMAKE_SOURCE_DIR}/databases/reaktoro/supcrtbl.json
        -oyaml ${CMAKE_SOURCE_DIR}/databases/reaktoro/supcrtbl.yaml
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
