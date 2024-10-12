cmake_minimum_required(VERSION 3.28)
project(MES_proj)

set(CMAKE_CXX_STANDARD 17)

add_executable(MES_proj main.cpp
        structures/GlobalData.cpp
        structures/GlobalData.h
        structures/structs.cpp
        structures/structs.h
        structures/grid.cpp
        structures/grid.h
)
