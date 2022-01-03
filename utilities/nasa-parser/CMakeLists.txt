add_custom_target(update-databases-nasa
    COMMENT "Updating nasa json and yaml databases..."
    COMMAND ${CMAKE_COMMAND} -E echo "This operation expects nasa-cea.dat has been updated!"
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/generate-nasa-databases
        -i     ${CMAKE_CURRENT_SOURCE_DIR}/nasa-cea.dat
        -ojson ${CMAKE_SOURCE_DIR}/databases/reaktoro/nasa-cea.json
        -oyaml ${CMAKE_SOURCE_DIR}/databases/reaktoro/nasa-cea.yaml
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
