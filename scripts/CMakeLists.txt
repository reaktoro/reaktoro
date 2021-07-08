# Copy the bash/batch scripts to build directory for easy setup of environment variables.
# Create target "envs" for Unix and Windows to export the environment variables

if(UNIX)
    configure_file(envs.in ${PROJECT_BINARY_DIR}/envs)
    add_custom_target(envs
        COMMAND bash -c "source envs"
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR})
endif()

if(WIN32)
    file(COPY envs4debug.bat DESTINATION ${PROJECT_BINARY_DIR})
    file(COPY envs4release.bat DESTINATION ${PROJECT_BINARY_DIR})
    add_custom_target(envs
        COMMAND cmd /c "$<$<CONFIG:Release>:envs4release.bat>$<$<CONFIG:Debug>:envs4debug.bat>"
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR})
endif()
