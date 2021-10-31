# Copy the bash/batch scripts to build directory for easy setup of environment variables.
# Create target "envs" for Unix and Windows to export the environment variables

if(UNIX)
    configure_file(envs.in ${PROJECT_BINARY_DIR}/envs)
    add_custom_target(envs
        COMMENT "ATTENTION: Execute `source envs` instead! CMake will not propagate the changes below to its parent process."
        COMMAND bash -c "source envs"
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR})
endif()

if(WIN32)
    # Create a native path for binary dir to be used when configuring envs4debug.bat.in and envs4release.bat.in
    file(TO_NATIVE_PATH ${PROJECT_BINARY_DIR} PROJECT_BINARY_DIR_NATIVE)

    configure_file(envs4debug.bat.in ${PROJECT_BINARY_DIR}/envs4debug.bat)
    configure_file(envs4release.bat.in ${PROJECT_BINARY_DIR}/envs4release.bat)

    add_custom_target(envs
        COMMENT "ATTENTION: Execute `envs4debug.bar` or `envs4release.bar` instead! CMake will not propagate the changes below to its parent process."
        COMMAND cmd /c "$<$<CONFIG:Release>:envs4release.bat>$<$<CONFIG:Debug>:envs4debug.bat>"
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR})
endif()
