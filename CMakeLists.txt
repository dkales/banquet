cmake_minimum_required (VERSION 3.12)
project (banquet C CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING "Choose the type of build." FORCE)
endif()

include(CheckCCompilerFlag)

set(BUILD_TESTS ON CACHE BOOL "Build unit tests.")

SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
FIND_PACKAGE(NTL REQUIRED)
FIND_PACKAGE(GMP REQUIRED)

add_subdirectory(keccak)

set(BANQUET_SRCS
  aes.cpp
  banquet.cpp
  banquet_instances.cpp
  tree.cpp
  tape.cpp
  utils.cpp
  )

add_library(banquet_static STATIC ${BANQUET_SRCS})
target_link_libraries(banquet_static PUBLIC ${NTL_LIB} keccak)

#add_library(bench_utils STATIC tools/bench_utils.c tools/bench_timing.c)
#add_executable(bench tools/bench.c)
#target_link_libraries(bench picnic bench_utils)

if(BUILD_TESTS)
find_package(Catch2 REQUIRED)
add_executable(aes_test tests/aes_test.cpp)
add_executable(tree_test tests/tree_test.cpp)
add_executable(banquet_test tests/banquet_test.cpp)
target_link_libraries(aes_test Catch2::Catch2 banquet_static)
target_link_libraries(tree_test Catch2::Catch2 banquet_static)
target_link_libraries(banquet_test Catch2::Catch2 banquet_static)

include(CTest)
include(Catch)
catch_discover_tests(tree_test)
catch_discover_tests(aes_test)
catch_discover_tests(banquet_test)
endif()