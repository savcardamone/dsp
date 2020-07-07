cmake_minimum_required(VERSION 3.16)
project(dsp)

list(APPEND CMAKE_MODULE_PATH /home/sav/eigen/cmake)

set(BOOST_INCLUDEDIR /home/sav/boost)
set(BOOST_LIBRARYDIR /home/sav/boost/stage/lib)
set(EIGEN3_INCLUDE_DIR /home/sav/eigen)

find_package(Eigen3 3.3 REQUIRED)
if(NOT EIGEN3_FOUND)
    message(STATUS "Eigen3 Not Found.")
endif()
find_package(Boost 1.73 REQUIRED)
if(NOT Boost_FOUND)
    message(STATUS "Boost Not Found.")
endif()

add_library(${PROJECT_NAME} INTERFACE)

target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_20)

enable_testing()
add_subdirectory(test)

