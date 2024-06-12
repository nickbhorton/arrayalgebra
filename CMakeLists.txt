cmake_minimum_required(VERSION 3.29)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED On)

project(
    arrayalgebra
    LANGUAGES CXX
)

add_compile_options(
    -g
    -O3
    -Wall
    -Werror
    -Wextra
)
