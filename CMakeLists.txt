cmake_minimum_required(VERSION 3.26)
project(casper LANGUAGES CXX)

find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
find_package(Boost REQUIRED)
add_library(casper SHARED nbody_prop.cpp)
target_include_directories(casper PUBLIC "pybind11" "${CMAKE_SOURCE_DIR}/cspice/include" "${CMAKE_SOURCE_DIR}cspice/lib/cspice.a")
target_compile_definitions(casper PRIVATE PYTHON_COMPILE)
target_link_libraries(casper PRIVATE "${CMAKE_SOURCE_DIR}/cspice/include/SpiceUsr.h" "${CMAKE_SOURCE_DIR}/cspice/lib/cspice.a")

set_property(TARGET casper PROPERTY PREFIX "")


add_executable(testing nbody_prop.cpp)
target_include_directories(testing PUBLIC "${CMAKE_SOURCE_DIR}/cspice/include" "${CMAKE_SOURCE_DIR}cspice/lib/cspice.a")
target_link_libraries(testing PRIVATE "${CMAKE_SOURCE_DIR}/cspice/include/SpiceUsr.h" "${CMAKE_SOURCE_DIR}/cspice/lib/cspice.a")

