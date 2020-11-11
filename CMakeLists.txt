cmake_minimum_required (VERSION 3.12)
project (banquet C CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 20) # for span

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING "Choose the type of build." FORCE)
endif()

include(CheckCCompilerFlag)

set(BUILD_TESTS ON CACHE BOOL "Build unit tests.")

SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
FIND_PACKAGE(NTL REQUIRED)
FIND_PACKAGE(GMP REQUIRED)
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

add_subdirectory(keccak)

set(BANQUET_SRCS
  aes.cpp
  banquet.cpp
  banquet_instances.cpp
  tree.cpp
  tape.cpp
  utils.cpp
  randomness.c
  )

add_library(banquet_static STATIC ${BANQUET_SRCS})
target_link_libraries(banquet_static PUBLIC keccak ${NTL_LIBRARY}  ${GMP_LIBRARIES} Threads::Threads)
target_include_directories(banquet_static PUBLIC ${NTL_INCLUDE_DIR})

add_library(bench_utils STATIC tools/bench_utils.cpp tools/bench_timing.cpp)
add_executable(bench tools/bench.cpp)
target_link_libraries(bench banquet_static bench_utils)

if(BUILD_TESTS)
Include(FetchContent)

FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v2.13.1)

FetchContent_MakeAvailable(Catch2)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${Catch2_SOURCE_DIR}/contrib/")

add_executable(aes_test tests/aes_test.cpp)
add_executable(tree_test tests/tree_test.cpp)
add_executable(banquet_test tests/banquet_test.cpp)
add_executable(util_test tests/util_test.cpp)
target_link_libraries(aes_test Catch2::Catch2 banquet_static)
target_link_libraries(tree_test Catch2::Catch2 banquet_static)
target_link_libraries(banquet_test Catch2::Catch2 banquet_static)
target_link_libraries(util_test Catch2::Catch2 banquet_static)

include(CTest)
include(Catch)
catch_discover_tests(tree_test)
catch_discover_tests(aes_test)
catch_discover_tests(banquet_test)
catch_discover_tests(util_test)
endif()