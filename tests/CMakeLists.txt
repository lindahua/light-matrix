# The CMakeLists file for Light-Matrix (Testing)

cmake_minimum_required(VERSION 2.8)
project(LightMatrix)

enable_testing()

include_directories(..)
include_directories($ENV{LIGHT_TEST_HOME})

# Compiler configuration

set(HAS_AVX "YES")

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

set(CMAKE_MODULE_PATH "../cmake_modules" ${CMAKE_MODULE_PATH})

# Intel SVML 

set(SVML_PATH "$ENV{ICC_LIBPATH}")

if (SVML_PATH)

    if (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    	find_library(SVML svml_disp
    		HINTS ${SVML_PATH})
    	find_library(LIBIRC libirc 
    		HINTS ${SVML_PATH})
    	
    else (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    	find_library(SVML svml 
    		HINTS ${SVML_PATH})
    
    endif (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    
    if (SVML_NOT_FOUND)
        message(SEND_ERROR "Intel SVML (Short Vector Math Library) is NOT found.")
    endif (SVML_NOT_FOUND)

endif (SVML_PATH)

# AMD LibM

set(LIBM_HOME "$ENV{AMDLIBM_HOME}")

if ( LIBM_HOME )

    set(AMDLIBM_PATH "$ENV{AMDLIBM_HOME}/lib/dynamic")

    find_library(AMDLIBM amdlibm
        HINTS ${AMDLIBM_PATH})
    
    if (AMDLIBM_NOT_FOUND)
        message(SEND_ERROR "AMD LibM is NOT found.")
    endif (AMDLIBM_NOT_FOUND)

endif ( LIBM_HOME )


# Intel MKL

find_package(MKL)


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
    ${INC}/math/math_constants.h
    ${INC}/math/fun_tags.h
    ${INC}/math/math_functors.h)
    
set(SIMD_BASE_HS
    ${INC}/math/simd_arch.h
    ${INC}/math/simd_base.h
    ${INC}/math/simd_debug.h
    ${INC}/math/internal/numrepr_format.h)
    
set(SSE_HS
    ${INC}/math/internal/sse_helpers.h
    ${INC}/math/internal/sse_testz_impl.h
    ${INC}/math/internal/sse2_round_impl.h
    ${INC}/math/internal/sse_math_emulate.h
    ${INC}/math/sse_packs.h
    ${INC}/math/sse_bpacks.h
    ${INC}/math/sse_ops.h
    ${INC}/math/sse_math.h
    ${INC}/math/sse_reduce.h
    ${INC}/math/sse.h)
    
set(AVX_HS
    ${INC}/math/internal/avx_helpers.h
    ${INC}/math/internal/avx_math_emulate.h
    ${INC}/math/avx_packs.h
    ${INC}/math/avx_bpacks.h
    ${INC}/math/avx_ops.h
    ${INC}/math/avx_math.h
    ${INC}/math/avx_reduce.h
    ${INC}/math/avx.h)
  
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
    
set(MATRIX_EVAL_HS
    ${INC}/mateval/internal/ewise_eval_internal.h
    ${INC}/mateval/mateval_fwd.h
    ${INC}/mateval/fun_maps.h
    ${INC}/mateval/common_kernels.h
    ${INC}/mateval/vec_accessors.h
    ${INC}/mateval/multicol_accessors.h
    ${INC}/mateval/ewise_eval.h)
    
set(MAP_EXPR_HS
    ${INC}/mateval/internal/map_expr_internal.h
    ${INC}/mateval/map_accessors.h
    ${INC}/mateval/map_expr.h
    ${INC}/mateval/mat_arith.h
    ${INC}/mateval/mat_emath.h)
    
set(MATRIX_REDUC_HS
    ${INC}/mateval/internal/mat_fold_internal.h
    ${INC}/mateval/internal/mat_reduce_internal.h
    ${INC}/mateval/mat_fold.h
    ${INC}/mateval/mat_reduce.h)
                  
set(MATRIX_HS
    ${COMMON_HS}
    ${MATRIX_BASE_HS}
    ${MATRIX_CLASS_HS})    
    
set(BLAS_HS
    ${INC}/linalg/linalg_fwd.h
    ${INC}/linalg/internal/linalg_aux.h
    ${INC}/linalg/blas_l1.h
    ${INC}/linalg/blas_l2.h
    ${INC}/linalg/blas_l3.h
    ${INC}/linalg/blas.h)    
    
set(LAPACK_HS
    ${INC}/linalg/lapack_fwd.h
    ${INC}/linalg/lapack_lu.h
    ${INC}/linalg/lapack_chol.h
    ${INC}/linalg/lapack_qr.h
    ${INC}/linalg/lapack_syev.h
    ${INC}/linalg/lapack_svd.h)
    
# Executables

add_library(test_main STATIC test_main.cpp)

add_executable(test_memory ${COMMON_HS} common/test_memory.cpp)
add_executable(test_blocks ${COMMON_HS} common/test_blocks.cpp)

set(SSE_TEST_HS
    ${COMMON_HS}
    ${SIMD_BASE_HS}
    ${SSE_HS})

add_executable(test_sse_packs  ${SSE_TEST_HS} math/test_sse_packs.cpp)
add_executable(test_sse_bpacks ${SSE_TEST_HS} math/test_sse_bpacks.cpp)
add_executable(test_sse_ops    ${SSE_TEST_HS} math/test_sse_ops.cpp)
add_executable(test_sse_round  ${SSE_TEST_HS} math/test_sse_round.cpp)
add_executable(test_sse_reduce ${SSE_TEST_HS} math/test_sse_reduce.cpp)

set(AVX_TEST_HS
    ${COMMON_HS}
    ${SIMD_BASE_HS}
    ${AVX_HS})

if (HAS_AVX)
add_executable(test_avx_packs  ${AVX_TEST_HS} math/test_avx_packs.cpp)
add_executable(test_avx_bpacks ${AVX_TEST_HS} math/test_avx_bpacks.cpp)
add_executable(test_avx_ops    ${AVX_TEST_HS} math/test_avx_ops.cpp)
add_executable(test_avx_round  ${AVX_TEST_HS} math/test_avx_round.cpp)
add_executable(test_avx_reduce ${AVX_TEST_HS} math/test_avx_reduce.cpp)
endif (HAS_AVX)

if (SVML_PATH)

add_executable(test_sse_math_svml ${SSE_TEST_HS} math/test_sse_math.cpp)
if (HAS_AVX)
add_executable(test_avx_math_svml ${AVX_TEST_HS} math/test_avx_math.cpp)
endif (HAS_AVX)

set_target_properties(test_sse_math_svml
	PROPERTIES
	COMPILE_FLAGS "-DLMAT_USE_INTEL_SVML"
)

if (HAS_AVX)
set_target_properties(test_avx_math_svml
	PROPERTIES
	COMPILE_FLAGS "-DLMAT_USE_INTEL_SVML"
)
endif (HAS_AVX)

if (MSVC)
target_link_libraries(test_sse_math_svml test_main ${SVML} ${LIBIRC})
if (HAS_AVX)
target_link_libraries(test_avx_math_svml test_main ${SVML} ${LIBIRC})
endif (HAS_AVX)
else (MSVC)
target_link_libraries(test_sse_math_svml test_main ${SVML})
if (HAS_AVX)
target_link_libraries(test_avx_math_svml test_main ${SVML})
endif (HAS_AVX)
endif (MSVC)

endif (SVML_PATH)


set(DENSE_MAT_HS
    ${COMMON_HS}
    ${MATRIX_BASE_HS}
    ${INC}/matrix/dense_matrix.h)

add_executable(test_dense_mat ${DENSE_MAT_HS} matrix/test_dense_mat.cpp)
add_executable(test_dense_vec ${DENSE_MAT_HS} matrix/test_dense_vec.cpp) 

set(REF_MAT_HS
    ${COMMON_HS}
    ${MATRIX_BASE_HS}
    ${INC}/matrix/regular_mat_base.h
    ${INC}/matrix/dense_matrix.h
    ${INC}/matrix/ref_matrix.h
    ${INC}/matrix/ref_block.h
    ${INC}/matrix/ref_grid.h
    ${INC}/matrix/step_vecs.h)
    
add_executable(test_ref_mat   ${REF_MAT_HS} matrix/test_ref_mat.cpp)
add_executable(test_ref_vec   ${REF_MAT_HS} matrix/test_ref_vec.cpp)
add_executable(test_ref_block ${REF_MAT_HS} matrix/test_ref_block.cpp)
add_executable(test_ref_grid  ${REF_MAT_HS} matrix/test_ref_grid.cpp)
add_executable(test_step_vec  ${REF_MAT_HS} matrix/test_step_vec.cpp)

add_executable(test_mat_props  ${REF_MAT_HS} matrix/test_mat_props.cpp)
add_executable(test_mat_fill   ${REF_MAT_HS} matrix/test_mat_fill.cpp)
add_executable(test_mat_copy   ${REF_MAT_HS} matrix/test_mat_copy.cpp)
add_executable(test_dense_eval ${MATRIX_HS}  matrix/test_dense_eval.cpp)

# suppress bogus free-nonheap-object complaint under GCC
if (${CMAKE_CXX_COMPILER_ID} MATCHES "GNU")
set_target_properties(test_dense_eval PROPERTIES COMPILE_FLAGS "-Wno-free-nonheap-object")
endif (${CMAKE_CXX_COMPILER_ID} MATCHES "GNU")

add_executable(test_mat_vecviews ${MATRIX_HS} matrix/test_mat_vecviews.cpp)
add_executable(test_mat_matviews ${MATRIX_HS} matrix/test_mat_matviews.cpp)

set(EWISE_EVAL_HS
    ${MATRIX_HS}
    ${MATRIX_EVAL_HS})
    
add_executable(test_linear_ewise ${EWISE_EVAL_HS} mateval/test_linear_ewise.cpp)    
add_executable(test_percol_ewise ${EWISE_EVAL_HS} mateval/test_percol_ewise.cpp)
add_executable(test_ewise_accum ${EWISE_EVAL_HS}  mateval/test_ewise_accum.cpp)

set(MAP_EXPR_TEST_HS
    ${MATRIX_HS}
    ${MATRIX_EVAL_HS}
    ${MAP_EXPR_HS})
    
add_executable(test_map_expr ${MAP_EXPR_TEST_HS} mateval/test_map_expr.cpp)
add_executable(test_mat_arith ${MAP_EXPR_TEST_HS} mateval/test_mat_arith.cpp)
add_executable(test_mat_emath ${MAP_EXPR_TEST_HS} mateval/test_mat_emath.cpp)
add_executable(test_mat_cast ${MAP_EXPR_TEST_HS} mateval/test_mat_cast.cpp)

if (SVML_PATH)
set_target_properties(test_mat_emath
	PROPERTIES
	COMPILE_FLAGS "-DLMAT_USE_INTEL_SVML"
)
if (MSVC)
target_link_libraries(test_mat_emath test_main ${SVML} ${LIBIRC})
else (MSVC)
target_link_libraries(test_mat_emath test_main ${SVML})
endif (MSVC)
endif (SVML_PATH)


set(REDUC_HS
    ${MATRIX_HS}
    ${MATRIX_EVAL_HS}
    ${MATRIX_REDUC_HS})

add_executable(test_full_reduce    ${REDUC_HS} mateval/test_full_reduce.cpp)
add_executable(test_colwise_reduce ${REDUC_HS} mateval/test_colwise_reduce.cpp)
add_executable(test_rowwise_reduce ${REDUC_HS} mateval/test_rowwise_reduce.cpp)

set(BLAS_TEST_HS
    ${MATRIX_HS}
    ${BLAS_HS})

add_executable(test_blas_l1 ${BLAS_TEST_HS} linalg/test_blas_l1.cpp)
add_executable(test_blas_l2 ${BLAS_TEST_HS} linalg/test_blas_l2.cpp)
add_executable(test_blas_l3 ${BLAS_TEST_HS} linalg/test_blas_l3.cpp)

set(LAPACK_TEST_HS
    ${MATRIX_HS}
    ${BLAS_HS}
    ${LAPACK_HS})
    
add_executable(test_lapack_lu ${LAPACK_TEST_HS} linalg/test_lapack_lu.cpp)
add_executable(test_lapack_chol ${LAPACK_TEST_HS} linalg/test_lapack_chol.cpp)
add_executable(test_lapack_qr ${LAPACK_TEST_HS} linalg/test_lapack_qr.cpp)
add_executable(test_lapack_syev ${LAPACK_TEST_HS} linalg/test_lapack_syev.cpp)
add_executable(test_lapack_svd ${LAPACK_TEST_HS} linalg/test_lapack_svd.cpp)

# Groups of test units

# common

set(LMAT_COMMON_TESTS
    test_memory
    test_blocks)
    
# math    
    
set(LMAT_SSE_TESTS
    test_sse_packs
    test_sse_bpacks
    test_sse_ops
    test_sse_round
    test_sse_reduce)

if (HAS_AVX)
set(LMAT_AVX_TESTS
    test_avx_packs
    test_avx_bpacks
    test_avx_ops
    test_avx_round
    test_avx_reduce)
endif (HAS_AVX)

if (HAS_AVX)
    set(LMAT_MATH_TESTS
        ${LMAT_SSE_TESTS}
        ${LMAT_AVX_TESTS})
else (HAS_AVX)
    set(LMAT_MATH_TESTS
        ${LMAT_SSE_TESTS})
endif (HAS_AVX)        

# matrix

set(LMAT_MATRIX_TESTS
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
	test_mat_matviews)

# mateval
	
set(LMAT_MATEVAL_TESTS
    test_linear_ewise
	test_percol_ewise
	test_ewise_accum
	test_map_expr
	test_mat_arith
	test_mat_emath
	test_mat_cast
	test_full_reduce
	test_colwise_reduce
	test_rowwise_reduce)

# linalg

set(LMAT_LINALG_TESTS
    test_blas_l1
    test_blas_l2
    test_blas_l3
    test_lapack_lu
    test_lapack_chol
    test_lapack_qr
    test_lapack_syev
    test_lapack_svd)
    
# Link to MKL

foreach(tname ${LMAT_LINALG_TESTS})
    target_link_libraries(${tname} ${MKL_LIBRARY})
endforeach(tname)

# all

set(LMAT_TEST_NAMES
    ${LMAT_COMMON_TESTS}
    ${LMAT_MATH_TESTS}
    ${LMAT_MATRIX_TESTS}
    ${LMAT_MATEVAL_TESTS}
    ${LMAT_LINALG_TESTS}
)

# Link to test_main
	
foreach(tname ${LMAT_TEST_NAMES})
	target_link_libraries(${tname} test_main)
endforeach(tname)	


# Add tests

foreach(tname ${LMAT_TEST_NAMES})
	add_test(NAME ${tname} COMMAND ${tname})
endforeach(tname)

if (SVML_PATH)
add_test(NAME test_sse_math_svml COMMAND test_sse_math_svml)
if (HAS_AVX)
add_test(NAME test_avx_math_svml COMMAND test_avx_math_svml)
endif (HAS_AVX)
endif (SVML_PATH)




