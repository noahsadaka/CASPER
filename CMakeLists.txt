cmake_minimum_required(VERSION 3.26)
project(spicetest LANGUAGES CXX)

find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
find_package(Boost REQUIRED)
add_library(spicetest SHARED spicetest.cpp)
target_include_directories(spicetest PUBLIC "pybind11" "${CMAKE_SOURCE_DIR}/cspice/include" "${CMAKE_SOURCE_DIR}cspice/lib/cspice.a")
target_link_libraries(spicetest PRIVATE "${CMAKE_SOURCE_DIR}/cspice/include/SpiceUsr.h" "${CMAKE_SOURCE_DIR}/cspice/lib/cspice.a")
set_property(TARGET spicetest PROPERTY PREFIX "")
