cmake_minimum_required(VERSION 3.10)
project(Physica)

# Set C++ Standard
set(CMAKE_CXX_STANDARD 11)

# Add include directory
include_directories(include)

# Add source files for the library
add_library(physica
    src/Vector2D.cpp
    src/Physics.cpp
)

# Create an executable for the example program
add_executable(example examples/main.cpp)

# Link the example with the Physica library
target_link_libraries(example physica)
