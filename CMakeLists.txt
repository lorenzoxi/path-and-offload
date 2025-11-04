cmake_minimum_required(VERSION 3.16)
project(2025-journal-drone CXX)

# Use the Release MSVC CRT even in Debug (Gurobi ships release libs)
set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreadedDLL")  # /MD for all configs
SET(CMAKE_CXX_STANDARD 17)

set(GUROBI_ROOT "C:/gurobi1203/win64" CACHE PATH "C:/gurobi1203/win64")
include_directories("${GUROBI_ROOT}/include")

set(GUROBI_LIB_DIR "${GUROBI_ROOT}/lib")
# Pick the C++ wrapper that exists in your lib folder; try 2022/2019/2017 in order
find_library(GUROBI_CPP_LIB NAMES gurobi_c++md2022 gurobi_c++md2019 gurobi_c++md2017
        PATHS "${GUROBI_LIB_DIR}" NO_DEFAULT_PATH REQUIRED)
find_library(GUROBI_CORE_LIB NAMES gurobi120
        PATHS "${GUROBI_LIB_DIR}" NO_DEFAULT_PATH REQUIRED)

add_executable(2025-journal-drone
        main.cpp
        energy.h)
target_link_libraries(2025-journal-drone "${GUROBI_CPP_LIB}" "${GUROBI_CORE_LIB}")
