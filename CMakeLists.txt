cmake_minimum_required(VERSION 3.7)

project(optimizations)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "-O3 -march=native")

#find_package(OpenMP)
#target_link_libraries(test OpenMP::OpenMP_CXX)

add_executable(test test.cpp)
target_link_libraries(test)
