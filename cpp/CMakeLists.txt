cmake_minimum_required(VERSION 3.14)
project(anf_simplify VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Release)

# ---- Library ----
add_library(anf_core
    src/sparse_anf.cpp
    src/gf2_linalg.cpp
    src/simplify.cpp
)

target_include_directories(anf_core PUBLIC include)

# ---- Main executable ----
add_executable(anf_simplify
    src/main.cpp
    src/sparse_anf.cpp
    src/gf2_linalg.cpp
    src/simplify.cpp
)

target_include_directories(anf_simplify PRIVATE include)

# ---- Benchmarks ----
add_executable(anf_bench bench/bench_main.cpp)
target_include_directories(anf_bench PRIVATE include)
target_link_libraries(anf_bench anf_core)

# ---- Tests ----
add_executable(anf_test test/test_main.cpp)
target_include_directories(anf_test PRIVATE include)
target_link_libraries(anf_test anf_core)

# ---- Build options ----
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(anf_core PRIVATE -O3 -march=native -mtune=native)
    target_compile_options(anf_simplify PRIVATE -O3 -march=native -mtune=native)
endif()
