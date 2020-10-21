cmake_minimum_required (VERSION 3.12)
project (banquet C CXX)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING "Choose the type of build." FORCE)
endif()

include(CheckCCompilerFlag)

SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
FIND_PACKAGE(NTL REQUIRED)
FIND_PACKAGE(GMP REQUIRED)


set(SHA3_IMPL "opt64" CACHE STRING "Select SHA3 implementation.")

include_directories(${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR} sha3 "sha3/${SHA3_IMPL}")

list(APPEND SHA3_SOURCES
     sha3/KeccakHash.c
     sha3/KeccakSpongeWidth1600.c)

if(${SHA3_IMPL} STREQUAL "avx2")
  list(APPEND SHA3_SOURCES sha3/avx2/KeccakP-1600-AVX2.s)
  set_property(SOURCE sha3/avx2/KeccakP-1600-AVX2.s PROPERTY LANGUAGE C)
  # the times4 variant of Keccack using avx2
  list(APPEND SHA3_SOURCES
       sha3/avx2/KeccakP-1600-times4-SIMD256.c
       sha3/KeccakSpongeWidth1600times4.c
       sha3/KeccakHashtimes4.c)
elseif(${SHA3_IMPL} STREQUAL "opt64")
  list(APPEND SHA3_SOURCES sha3/opt64/KeccakP-1600-opt64.c)
elseif(${SHA3_IMPL} STREQUAL "armv8a-neon")
  list(APPEND SHA3_SOURCES sha3/armv8a-neon/KeccakP-1600-armv8a-neon.s)
  set_property(SOURCE sha3/armv8a-neon/KeccakP-1600-armv8a-neon.s PROPERTY LANGUAGE C)
else()
  message(FATAL_ERROR "Unknown SHA3 implementation")
endif()

FILE(GLOB SOURCES "*.c" "*.h")

add_library(banquet_static STATIC ${SOURCES} ${SHA3_SOURCES})
target_link_libraries(banquet_static PUBLIC ${NTL_LIB})
add_library(bench_utils STATIC tools/bench_utils.c tools/bench_timing.c)
add_executable(bench tools/bench.c)
target_link_libraries(bench picnic bench_utils)