cmake_minimum_required(VERSION 3.14)
project(optimize_anf VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Release)

# Common flags
set(COMMON_FLAGS "-O3 -march=native -mtune=native -pthread")
set(CMAKE_CXX_FLAGS_RELEASE "${COMMON_FLAGS}")

# Core sources (truth-table pipeline)
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

# ====================================================================
# Direction 1a: Truth-table-based search
#   opt  = single-output (k=1, shared=per-output are equivalent)
#   opt1 = multi-output, shared transform
#   opt2 = multi-output, per-output transforms (large-n z-space eval)
# ====================================================================

add_executable(optimize_anf_d1a_opt
    src/main_opt1.cpp
    ${CORE_SRC}
)
target_include_directories(optimize_anf_d1a_opt PRIVATE include)
target_compile_definitions(optimize_anf_d1a_opt PRIVATE DIRECTION_TAG="d1a" STRATEGY_TAG="opt")

add_executable(optimize_anf_d1a_opt1
    src/main_opt1.cpp
    ${CORE_SRC}
)
target_include_directories(optimize_anf_d1a_opt1 PRIVATE include)
target_compile_definitions(optimize_anf_d1a_opt1 PRIVATE DIRECTION_TAG="d1a" STRATEGY_TAG="opt1")

add_executable(optimize_anf_d1a_opt2
    src/main_large.cpp
    ${CORE_SRC}
)
target_include_directories(optimize_anf_d1a_opt2 PRIVATE include)
target_compile_definitions(optimize_anf_d1a_opt2 PRIVATE STRATEGY_TAG="opt2")

# ====================================================================
# Direction 1b: Gate-builder-based construction
#   opt  = single-output (opt=opt2 when k=1)
#   opt2 = per-output transforms
#   Note: opt1 (shared) not applicable — gate builder is per-output by design
# ====================================================================

add_executable(optimize_anf_d1b_opt
    src/main_gate_builder.cpp
    src/sparse_anf.cpp
    src/gate_builder.cpp
    src/circuit.cpp
)
target_include_directories(optimize_anf_d1b_opt PRIVATE include)
target_compile_definitions(optimize_anf_d1b_opt PRIVATE STRATEGY_TAG="opt")
target_compile_options(optimize_anf_d1b_opt PRIVATE -O3 -march=native -mtune=native -pthread -std=c++17)

add_executable(optimize_anf_d1b_opt2
    src/main_gate_builder.cpp
    src/sparse_anf.cpp
    src/gate_builder.cpp
    src/circuit.cpp
)
target_include_directories(optimize_anf_d1b_opt2 PRIVATE include)
target_compile_definitions(optimize_anf_d1b_opt2 PRIVATE STRATEGY_TAG="opt2")
target_compile_options(optimize_anf_d1b_opt2 PRIVATE -O3 -march=native -mtune=native -pthread -std=c++17)

# ====================================================================
# Direction 2: Sparse-g enumeration (pool search + c-correction)
#   opt  = single-output (opt=opt2 when k=1)
#   opt2 = per-output transforms
#   Note: opt1 (shared) not applicable — sparse g is per-output by design
# ====================================================================

add_executable(optimize_anf_d2_opt
    src/main_sparse_g.cpp
    ${CORE_SRC}
)
target_include_directories(optimize_anf_d2_opt PRIVATE include)
target_compile_definitions(optimize_anf_d2_opt PRIVATE STRATEGY_TAG="opt")

add_executable(optimize_anf_d2_opt2
    src/main_sparse_g.cpp
    ${CORE_SRC}
)
target_include_directories(optimize_anf_d2_opt2 PRIVATE include)
target_compile_definitions(optimize_anf_d2_opt2 PRIVATE STRATEGY_TAG="opt2")

# ====================================================================
# Direction 3: Per-output Walsh + random + hill-climb
#   opt  = single-output (opt=opt2 when k=1)
#   opt1 = shared transform (uses same shared-search as d1a_opt1)
#   opt2 = per-output transforms
# ====================================================================

add_executable(optimize_anf_d3_opt
    src/main_opt2.cpp
    ${CORE_SRC}
)
target_include_directories(optimize_anf_d3_opt PRIVATE include)
target_compile_definitions(optimize_anf_d3_opt PRIVATE STRATEGY_TAG="opt")

add_executable(optimize_anf_d3_opt1
    src/main_opt1.cpp
    ${CORE_SRC}
)
target_include_directories(optimize_anf_d3_opt1 PRIVATE include)
target_compile_definitions(optimize_anf_d3_opt1 PRIVATE DIRECTION_TAG="d3" STRATEGY_TAG="opt1")

add_executable(optimize_anf_d3_opt2
    src/main_opt2.cpp
    ${CORE_SRC}
)
target_include_directories(optimize_anf_d3_opt2 PRIVATE include)
target_compile_definitions(optimize_anf_d3_opt2 PRIVATE STRATEGY_TAG="opt2")

# ====================================================================
# Direction 1c: Complement-only search (z = x ⊕ b, M = I)
#   Phase 1b in main_opt2.cpp, only identity + complement b search
#   opt  = single-output
#   opt2 = per-output transforms
#   opt1 (shared) not applicable — complement is per-output by nature
# ====================================================================

add_executable(optimize_anf_d1c_opt
    src/main_opt2.cpp
    ${CORE_SRC}
)
target_include_directories(optimize_anf_d1c_opt PRIVATE include)
target_compile_definitions(optimize_anf_d1c_opt PRIVATE DIRECTION_TAG="d1c" STRATEGY_TAG="opt")

add_executable(optimize_anf_d1c_opt2
    src/main_opt2.cpp
    ${CORE_SRC}
)
target_include_directories(optimize_anf_d1c_opt2 PRIVATE include)
target_compile_definitions(optimize_anf_d1c_opt2 PRIVATE DIRECTION_TAG="d1c" STRATEGY_TAG="opt2")

# ====================================================================
# Auxiliary tools
# ====================================================================

# Standalone raw ANF computation
add_executable(compute_raw_anf src/compute_raw_anf.cpp)
target_compile_options(compute_raw_anf PRIVATE -O3 -march=native -mtune=native -pthread)

# Standalone verifier: read circuit + trans + expr, verify f(x) == g(Mx+b)
add_executable(verify_anf
    src/verify_anf.cpp
    src/circuit.cpp
)
target_include_directories(verify_anf PRIVATE include)
target_compile_options(verify_anf PRIVATE -O3 -march=native -mtune=native -pthread)

# Preprocess: one-time variable remapping and output classification
add_executable(preprocess
    src/preprocess.cpp
    src/circuit.cpp
)
target_include_directories(preprocess PRIVATE include)
target_compile_options(preprocess PRIVATE -O3 -march=native -mtune=native -pthread -std=c++17)
