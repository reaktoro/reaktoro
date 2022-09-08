# Recursively collect all parameter database files from the current directory
file(GLOB_RECURSE FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *)

# Create a resource library containing embedded parameter databases
cmrc_add_resource_library(ReaktoroParams
    ALIAS Reaktoro::Params
    PREFIX params
    ${FILES}
)
