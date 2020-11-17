cmake_minimum_required (VERSION 3.12)
project (banquet C CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17) 

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
endif()

include(CheckCXXCompilerFlag)

set(BUILD_TESTS OFF CACHE BOOL "Build unit tests.")

SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

add_compile_options("-Wall")

add_subdirectory(keccak)

set(BANQUET_SRCS
  aes.cpp
  banquet.cpp
  banquet_instances.cpp
  field.cpp
  tree.cpp
  tape.cpp
  randomness.c
  )

set(PCLMUL_FILES field.cpp)
CHECK_CXX_COMPILER_FLAG(-mpclmul COMPILER_HAS_M_PCLMUL)
CHECK_CXX_COMPILER_FLAG(-msse2 COMPILER_HAS_M_SSE2)
CHECK_CXX_COMPILER_FLAG(-msse4 COMPILER_HAS_M_SSE4)
if (COMPILER_HAS_M_PCLMUL AND COMPILER_HAS_M_SSE2 AND COMPILER_HAS_M_SSE4)
  add_compile_options(-mpclmul -msse2 -msse4)
else()
  message(
    ERROR
    "compiler does not have at least one of flag (pclmul, sse2, sse4) which are needed"  )
endif()

add_library(banquet_static STATIC ${BANQUET_SRCS})
target_link_libraries(banquet_static PUBLIC keccak) 

add_library(bench_utils STATIC tools/bench_utils.cpp tools/bench_timing.cpp)
add_executable(bench tools/bench.cpp)
add_executable(bench_free tools/bench_free.cpp)
target_link_libraries(bench banquet_static bench_utils)
target_link_libraries(bench_free banquet_static bench_utils)

if(BUILD_TESTS)
FIND_PACKAGE(NTL REQUIRED)
FIND_PACKAGE(GMP REQUIRED)
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
Include(FetchContent)

FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v2.13.1)

FetchContent_MakeAvailable(Catch2)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${Catch2_SOURCE_DIR}/contrib/")
add_library(test_utils STATIC tests/utils.cpp)
target_link_libraries(test_utils PUBLIC ${NTL_LIBRARY} ${GMP_LIBRARIES} Threads::Threads banquet_static)

add_executable(aes_test tests/aes_test.cpp)
add_executable(tree_test tests/tree_test.cpp)
add_executable(banquet_test tests/banquet_test.cpp)
add_executable(util_test tests/util_test.cpp)
add_executable(field_test tests/field_test.cpp)
target_link_libraries(aes_test Catch2::Catch2 test_utils)
target_link_libraries(tree_test Catch2::Catch2 test_utils)
target_link_libraries(banquet_test Catch2::Catch2 test_utils)
target_link_libraries(util_test Catch2::Catch2 test_utils)
target_link_libraries(field_test Catch2::Catch2  test_utils)

include(CTest)
include(Catch)
catch_discover_tests(tree_test)
catch_discover_tests(aes_test)
catch_discover_tests(banquet_test)
catch_discover_tests(field_test)
catch_discover_tests(util_test)
endif()