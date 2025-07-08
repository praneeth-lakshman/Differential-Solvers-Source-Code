cmake_minimum_required(VERSION 3.13.0)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
project(ODESolvers VERSION 0.1.0 LANGUAGES C CXX)
add_subdirectory(pybind11)

pybind11_add_module(ODESolvers main.cpp)

