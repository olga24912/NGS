cmake_minimum_required(VERSION 3.9)
project(hw2_De_Bruijn_graph)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES main.cpp Hash.cpp Hash.h)
add_executable(hw2_De_Bruijn_graph ${SOURCE_FILES} DeBruijnGraph.cpp DeBruijnGraph.h GraphConstructor.cpp GraphConstructor.h)