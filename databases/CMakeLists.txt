# Recursively collect all database files from the current directory
file(GLOB_RECURSE FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *)

# Create a resource library containing embedded databases
cmrc_add_resource_library(ReaktoroDatabases
    ALIAS Reaktoro::Databases
    PREFIX databases
    ${FILES}
)

# Set some target properties
set_target_properties(ReaktoroDatabases PROPERTIES
    POSITION_INDEPENDENT_CODE ON)
