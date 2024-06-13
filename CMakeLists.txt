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

add_executable(
    mult
    benchmark/multiplication.cc
)
target_include_directories(
    mult
    PUBLIC
    ${CMAKE_SOURCE_DIR}
)
