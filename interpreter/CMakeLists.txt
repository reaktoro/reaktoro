# Build the Reaktoro interpreter
add_executable(ireaktoro ireaktoro.cpp)

# Link the Reaktoro interpreter against Reaktoro C++ library
target_link_libraries(ireaktoro Reaktoro)

# Install Reaktoro interpreter
install(TARGETS ireaktoro
    EXPORT ReaktoroTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT interpreter
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT interpreter
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT interpreter)
