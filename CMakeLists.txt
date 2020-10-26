cmake_minimum_required (VERSION 3.12)
project (banquet C CXX)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING "Choose the type of build." FORCE)
endif()

include(CheckCCompilerFlag)

SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
FIND_PACKAGE(NTL REQUIRED)
FIND_PACKAGE(GMP REQUIRED)

add_subdirectory(keccak)

set(BANQUET_SRCS
  aes.c
  banquet.c
  banquet_instances.c
  tree.c
  utils.c
  )

add_library(banquet_static STATIC ${BANQUET_SRCS})
target_link_libraries(banquet_static PUBLIC ${NTL_LIB} keccak)

add_library(bench_utils STATIC tools/bench_utils.c tools/bench_timing.c)
add_executable(bench tools/bench.c)
target_link_libraries(bench picnic bench_utils)