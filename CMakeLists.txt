    

cmake_minimum_required (VERSION 2.8)

project (ENU_K_BUNDLE)

aux_source_directory(. DIR_SRCS)

set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -Wall -g -ggdb -pg")
set(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 -Wall")

add_executable(enplex ${DIR_SRCS})