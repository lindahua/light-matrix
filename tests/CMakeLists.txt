# The CMakeLists file for Light-Matrix (Testing)

cmake_minimum_required(VERSION 2.8)
project(LightMatrix)

enable_testing()

include_directories(..)
include_directories($ENV{LIGHT_TEST_HOME})

# Compiler configuration

set(CMAKE_BUILD_TYPE "Release")

if (MSVC)
	set(LANG_FLAGS "/arch:AVX /EHsc")
	set(WARNING_FLAGS "/W4")
else (MSVC)
	set(LANG_FLAGS "-std=c++0x -pedantic -mavx")
	set(WARNING_FLAGS "-Wall -Wextra -Wconversion -Wformat -Wno-unused-parameter ")
endif (MSVC)

if (${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
	set(LANG_FLAGS "${LANG_FLAGS} -stdlib=libc++ -Qunused-arguments")
	set(CMAKE_CXX_COMPILER "clang++")
endif (${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")

set(TEST_CXX_FLAGS "${LANG_FLAGS} ${WARNING_FLAGS}")

set(CMAKE_CXX_FLAGS "${TEST_CXX_FLAGS}")


# Find third-party libraries

# set(CMAKE_MODULE_PATH "../cmake_modules" ${CMAKE_MODULE_PATH})
# find_package(MKL)


# Header file groups (to be used as dependencies)

set(INC ../light_mat)
set(EXECUTABLE_OUTPUT_PATH bin)

set(COMMON_HS 
    ${INC}/config/user_config.h 
    ${INC}/config/platform_config.h
    ${INC}/config/config.h
    ${INC}/common/prim_types.h
    ${INC}/common/meta_base.h
    ${INC}/common/arg_check.h
    ${INC}/common/args_alg.h
    ${INC}/common/range.h
    ${INC}/common/basic_defs.h    
    ${INC}/common/memory.h
    ${INC}/common/memalloc.h
    ${INC}/common/block.h
    ${INC}/math/math_base.h
    ${INC}/math/math_constants.h)
    
set(SIMD_BASE_HS
    ${INC}/math/simd_arch.h
    ${INC}/math/simd_base.h
    ${INC}/math/simd_debug.h
    ${INC}/math/internal/numrepr_format.h)
    
set(SSE_HS
    ${INC}/math/internal/sse_helpers.h
    ${INC}/math/internal/sse_testz_impl.h
    ${INC}/math/internal/sse2_round_impl.h
    ${INC}/math/sse_packs.h
    ${INC}/math/sse_bpacks.h
    ${INC}/math/sse_ops.h
    ${INC}/math/sse_reduce.h)
    
set(AVX_HS
    ${INC}/math/internal/avx_helpers.h
    ${INC}/math/avx_packs.h
    ${INC}/math/avx_bpacks.h
    ${INC}/math/avx_ops.h)
  
set(MATRIX_BASE_HS
    ${INC}/matrix/matrix_shape.h
    ${INC}/matrix/matrix_layout.h
    ${INC}/matrix/matrix_fwd.h 
    ${INC}/matrix/matrix_meta.h
    ${INC}/matrix/matrix_shape.h
    ${INC}/matrix/matrix_concepts.h
    ${INC}/matrix/matrix_properties.h
    ${INC}/matrix/internal/matrix_fill_internal.h
    ${INC}/matrix/internal/matrix_copy_internal.h
    ${INC}/matrix/internal/matrix_colviews_internal.h
    ${INC}/matrix/internal/matrix_rowviews_internal.h
    ${INC}/matrix/internal/matrix_matviews_internal.h
    ${INC}/matrix/matrix_fill.h
    ${INC}/matrix/matrix_copy.h
    ${INC}/matrix/matrix_print.h
    ${INC}/matrix/matrix_subviews.h
    ${INC}/matrix/matrix_base.h)
    
set(MATRIX_CLASS_HS
    ${INC}/matrix/regular_mat_base.h
    ${INC}/matrix/dense_matrix.h
    ${INC}/matrix/ref_matrix.h
    ${INC}/matrix/ref_block.h
    ${INC}/matrix/ref_grid.h
    ${INC}/matrix/step_vecs.h
    ${INC}/matrix/dense_mutable_view.h
    ${INC}/matrix/matrix_classes.h)
              
set(MATRIX_HS
    ${COMMON_HS}
    ${MATRIX_BASE_HS}
    ${MATRIX_CLASS_HS})    
    
# Executables

add_library(test_main STATIC test_main.cpp)

add_executable(test_meta ${COMMON_HS} test_meta.cpp)
add_executable(test_memory ${COMMON_HS} test_memory.cpp)
add_executable(test_blocks ${COMMON_HS} test_blocks.cpp)

set(SSE_TEST_HS
    ${COMMON_HS}
    ${SIMD_BASE_HS}
    ${SSE_HS})

add_executable(test_sse_packs ${SSE_TEST_HS} test_sse_packs.cpp)
add_executable(test_sse_bpacks ${SSE_TEST_HS} test_sse_bpacks.cpp)
add_executable(test_sse_ops ${SSE_TEST_HS} test_sse_ops.cpp)
add_executable(test_sse_round ${SSE_TEST_HS} test_sse_round.cpp)
add_executable(test_sse_reduce ${SSE_TEST_HS} test_sse_reduce.cpp)

set(AVX_TEST_HS
    ${COMMON_HS}
    ${SIMD_BASE_HS}
    ${AVX_HS})

add_executable(test_avx_packs ${AVX_TEST_HS} test_avx_packs.cpp)
add_executable(test_avx_bpacks ${AVX_TEST_HS} test_avx_bpacks.cpp)
add_executable(test_avx_ops ${AVX_TEST_HS} test_avx_ops.cpp)
add_executable(test_avx_round ${AVX_TEST_HS} test_avx_round.cpp)

set(DENSE_MAT_HS
    ${COMMON_HS}
    ${MATRIX_BASE_HS}
    ${INC}/matrix/dense_matrix.h)

add_executable(test_dense_mat ${DENSE_MAT_HS} test_dense_mat.cpp)
add_executable(test_dense_vec ${DENSE_MAT_HS} test_dense_vec.cpp) 

set(REF_MAT_HS
    ${COMMON_HS}
    ${MATRIX_BASE_HS}
    ${INC}/matrix/regular_mat_base.h
    ${INC}/matrix/dense_matrix.h
    ${INC}/matrix/ref_matrix.h
    ${INC}/matrix/ref_block.h
    ${INC}/matrix/ref_grid.h
    ${INC}/matrix/step_vecs.h)
    
add_executable(test_ref_mat ${REF_MAT_HS} test_ref_mat.cpp)
add_executable(test_ref_vec ${REF_MAT_HS} test_ref_vec.cpp)
add_executable(test_ref_block ${REF_MAT_HS} test_ref_block.cpp)
add_executable(test_ref_grid ${REF_MAT_HS} test_ref_grid.cpp)
add_executable(test_step_vec ${REF_MAT_HS} test_step_vec.cpp)

add_executable(test_mat_props ${REF_MAT_HS} test_mat_props.cpp)
add_executable(test_mat_fill ${REF_MAT_HS} test_mat_fill.cpp)
add_executable(test_mat_copy ${REF_MAT_HS} test_mat_copy.cpp)

add_executable(test_dense_eval ${MATRIX_HS} test_dense_eval.cpp)

add_executable(test_mat_vecviews ${MATRIX_HS} test_mat_vecviews.cpp)
add_executable(test_mat_matviews ${MATRIX_HS} test_mat_matviews.cpp)


# suppress bogus free-nonheap-object complaint under GCC

if (${CMAKE_CXX_COMPILER_ID} MATCHES "GNU")
set_target_properties(test_dense_eval PROPERTIES COMPILE_FLAGS "-Wno-free-nonheap-object")
endif (${CMAKE_CXX_COMPILER_ID} MATCHES "GNU")

# Link to test_main

set(LMAT_TEST_NAMES
    test_meta
    test_memory
	test_blocks
	test_sse_packs
	test_sse_bpacks
	test_sse_ops
	test_sse_round
	test_sse_reduce
	test_avx_packs
	test_avx_bpacks
	test_avx_ops
	test_avx_round
	test_dense_mat
	test_dense_vec
	test_ref_mat
	test_ref_vec
	test_ref_block
	test_ref_grid
	test_step_vec
	test_mat_props
	test_mat_fill
	test_mat_copy
	test_dense_eval
	test_mat_vecviews
	test_mat_matviews
)
	
foreach(tname ${LMAT_TEST_NAMES})
	target_link_libraries(${tname} test_main)
endforeach(tname)	


# Add tests

foreach(tname ${LMAT_TEST_NAMES})
	add_test(NAME ${tname} COMMAND ${tname})
endforeach(tname)





