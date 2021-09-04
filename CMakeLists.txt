cmake_minimum_required(VERSION 3.11)
project(delaunator)

# C++ version
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(delaunator INTERFACE)
target_include_directories(delaunator INTERFACE .)

# add_executable(example example/example.cpp)
# target_link_libraries(example delaunator)