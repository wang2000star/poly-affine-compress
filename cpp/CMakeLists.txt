cmake_minimum_required(VERSION 3.14)
project(optimize_anf VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Release)

# Common flags
set(COMMON_FLAGS "-O3 -march=native -mtune=native -pthread")
set(CMAKE_CXX_FLAGS_RELEASE "${COMMON_FLAGS}")

# Core sources
set(CORE_SRC
    src/circuit.cpp
    src/truth_table.cpp
    src/moebius.cpp
    src/walsh.cpp
    src/gf2.cpp
    src/anf.cpp
    src/affine.cpp
    src/search.cpp
    src/io.cpp
    src/verify.cpp
)

# Main ANF optimizer
add_executable(optimize_anf
    src/main.cpp
    ${CORE_SRC}
)
target_include_directories(optimize_anf PRIVATE include)

# Legacy: compute_raw_anf (standalone, uses old code)
add_executable(compute_raw_anf src/compute_raw_anf.cpp)
target_compile_options(compute_raw_anf PRIVATE -O3 -march=native -mtune=native -pthread)

# Standalone verifier: read circuit + trans + expr, verify f(x) == g(Mx+b)
add_executable(verify_anf
    src/verify_anf.cpp
    src/circuit.cpp
)
target_include_directories(verify_anf PRIVATE include)
target_compile_options(verify_anf PRIVATE -O3 -march=native -mtune=native -pthread)
