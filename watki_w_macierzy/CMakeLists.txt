cmake_minimum_required (VERSION 2.6)
project (Studnia)

set(CMAKE_BUILD_TYPE Debug)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")

find_package (Threads)
find_package(GSL REQUIRED)
add_executable (prog studnia.c)
target_link_libraries (prog ${CMAKE_THREAD_LIBS_INIT} GSL::gsl)


