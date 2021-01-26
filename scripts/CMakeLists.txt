# Copy the bash/batch scripts to build directory for easy setup of environment variables.

if(UNIX)
    file(COPY envs DESTINATION ${CMAKE_BINARY_DIR})
endif()

if(WIN32)
    file(COPY envs4debug.bat DESTINATION ${CMAKE_BINARY_DIR})
    file(COPY envs4release.bat DESTINATION ${CMAKE_BINARY_DIR})
endif()
