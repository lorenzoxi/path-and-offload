cmake_minimum_required(VERSION 3.16)
project(path-and-offload CXX)

SET(CMAKE_CXX_STANDARD 17)

if (APPLE)
    set(CPLEX_ROOT_DIR "/Applications/CPLEX_Studio2211")
    set(CPLEX_INCLUDE_DIR "${CPLEX_ROOT_DIR}/cplex/include")
    set(CONCERT_INCLUDE_DIR "${CPLEX_ROOT_DIR}/concert/include")
    set(CPLEX_LIB_DIR "${CPLEX_ROOT_DIR}/cplex/lib/arm64_osx/static_pic")
    set(CONCERT_LIB_DIR "${CPLEX_ROOT_DIR}/concert/lib/arm64_osx/static_pic")
    set(Boost_INCLUDE_DIR "/usr/local/Cellar/boost/1.82.0_1/include")
    set(Boost_LIBRARY_DIR "/usr/local/Cellar/boost/1.82.0_1/lib")
    set(GUROBI_HOME "/Library/gurobi1102/macos_universal2")
    set(CGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE TRUE)

    include_directories(${GUROBI_HOME}/include)
    include_directories(${CPLEX_INCLUDE_DIR})
    include_directories(${CONCERT_INCLUDE_DIR})
    link_directories(${CPLEX_LIB_DIR})
    link_directories(${CONCERT_LIB_DIR})
    link_directories(${GUROBI_HOME}/lib)
    find_package(CGAL REQUIRED)


    add_executable(
            path-and-offload
            main.cpp
    )
    target_link_libraries(path-and-offload ilocplex concert cplex pthread gurobi_c++ gurobi110 CGAL::CGAL)

elseif (WIN32)
    #Windows Setup
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
else ()
    message(FATAL_ERROR "Unsupported platform")
endif ()