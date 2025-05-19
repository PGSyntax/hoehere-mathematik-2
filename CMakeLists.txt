cmake_minimum_required(VERSION 3.30)
project(hm2_3)

set(CMAKE_CXX_STANDARD 20)

add_executable(hm2_3 main.cpp CMyMatrix.cpp CMyVektor.cpp
        C_DGLSolver.cpp)
