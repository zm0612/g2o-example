cmake_minimum_required(VERSION 3.1)
project(untitled2)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_BUILD_TYPE Release)

set(ALL_TARGET_LIBRARIES "")

include(cmake/FindEigen3.cmake)
include(cmake/FindG2O.cmake)

add_executable(fit_curve fit_curve.cpp)
target_link_libraries(fit_curve ${ALL_TARGET_LIBRARIES})
