# Specify the minimum cmake version
cmake_minimum_required(VERSION 3.9)

# Specify the name of your project
project(reaktoro-app)

# Find the Reaktoro dependency.
find_package(Reaktoro REQUIRED)

# Create an executable named app with source code main.cpp
add_executable(app main.cpp)

# Link your application against Reaktoro C++ library.
target_link_libraries(app Reaktoro::Reaktoro)
