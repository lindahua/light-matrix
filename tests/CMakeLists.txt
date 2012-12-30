# The CMakeLists file for Light-Matrix (Testing)

cmake_minimum_required(VERSION 2.8)
project(LightMatrix)

enable_testing()

include_directories(..)
include_directories($ENV{LIGHT_TEST_HOME})

#==========================================================
#
#    compiler configurations
#
#==========================================================

include("../cmake_modules/DetectISA.cmake")

set(CMAKE_BUILD_TYPE "Release")

if (MSVC)
	set(LANG_FLAGS "${ARCH_FLAG} /EHsc")
	set(WARNING_FLAGS "/W4")
else (MSVC)
	set(LANG_FLAGS "-std=c++0x -pedantic -m64 ${ARCH_FLAG}")
	set(WARNING_FLAGS "-Wall -Wextra -Wconversion -Wformat -Wno-unused-parameter ")
endif (MSVC)

if (${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
	set(LANG_FLAGS "${LANG_FLAGS} -stdlib=libc++ -Qunused-arguments")
	set(CMAKE_CXX_COMPILER "clang++")
endif (${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")

set(TEST_CXX_FLAGS "${LANG_FLAGS} ${WARNING_FLAGS}")

set(CMAKE_CXX_FLAGS "${TEST_CXX_FLAGS}")

message(STATUS "[LMAT] CMAKE_CXX_COMPILER_ID = ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "[LMAT] CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")


#==========================================================
#
#    third-party library dependencies
#
#==========================================================

set(CMAKE_MODULE_PATH "../cmake_modules" ${CMAKE_MODULE_PATH})

# Intel SVML 

set(SVML_PATH "$ENV{ICC_LIBPATH}")

if (SVML_PATH)

    message(STATUS "[LMAT] SVML_PATH = ${SVML_PATH}")

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

else (SVML_PATH)
    message(STATUS "[LMAT] SVML is not used")
endif (SVML_PATH)

# AMD LibM

set(LIBM_PATH "$ENV{AMDLIBM_PATH}")

if ( LIBM_PATH )

    message(STATUS "[LMAT] LIBM_PATH = ${LIBM_PATH}")

    set(AMDLIBM_PATH "$ENV{AMDLIBM_PATH}/lib/dynamic")

    find_library(AMDLIBM amdlibm
        HINTS ${AMDLIBM_PATH})
    
    if (AMDLIBM_NOT_FOUND)
        message(SEND_ERROR "AMD LibM is NOT found.")
    endif (AMDLIBM_NOT_FOUND)

else (LIBM_PATH)
    
    message(STATUS "[LMAT] AMD LibM is not used.") 

endif ( LIBM_PATH )


# Intel MKL

find_package(MKL)

message(STATUS "[LMAT] Intel MKL found: ${MKLROOT_PATH}")


#==========================================================
#
#    Header file groups (to be used as dependencies)
#
#==========================================================

set(INC ../light_mat)
set(EXECUTABLE_OUTPUT_PATH bin)

# config

set(CONFIG_HS
    ${INC}/config/platform_config.h
    ${INC}/config/user_config.h
    ${INC}/config/config.h)

# common module

set(BASIC_DEFS_HS_
    ${INC}/common/prim_types.h
    ${INC}/common/mask_type.h
    ${INC}/common/meta_base.h
    ${INC}/common/int_div.h
    ${INC}/common/range.h
    ${INC}/common/arg_check.h
    ${INC}/common/args_alg.h
    ${INC}/common/basic_defs.h)
    
set(BASIC_MEM_HS_
    ${INC}/common/internal/align_alloc.h
    ${INC}/common/memory.h
    ${INC}/common/memalloc.h
    ${INC}/common/block.h)
    
set(COMMON_HS 
    ${BASIC_DEFS_HS_}
    ${BASIC_MEM_HS_})
    
set(COMMON_HS_EX
    ${CONFIG_HS}
    ${COMMON_HS})
    
# math module
    
set(MATH_FUNBASE_HS_
    ${INC}/math/math_base.h
    ${INC}/math/math_constants.h
    ${INC}/math/fun_tags.h
    ${INC}/math/math_functors.h)
    
set(SIMD_BASE_HS_
    ${INC}/math/internal/numrepr_format.h
    ${INC}/math/simd_arch.h
    ${INC}/math/simd_base.h
    ${INC}/math/simd_debug.h)
    
set(SSE_HS_
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
    
set(AVX_HS_
    ${INC}/math/internal/avx_helpers.h
    ${INC}/math/internal/avx_math_emulate.h
    ${INC}/math/avx_packs.h
    ${INC}/math/avx_bpacks.h
    ${INC}/math/avx_ops.h
    ${INC}/math/avx_math.h
    ${INC}/math/avx_reduce.h
    ${INC}/math/avx.h)
    
set(MATH_HS
    ${MATH_FUNBASE_HS_}
    ${SIMD_BASE_HS_}
    ${SSE_HS_}
    ${AVX_HS_}
    ${INC}/math/simd.h)    
    
set(MATH_HS_EX
    ${CONFIG_HS}
    ${COMMON_HS}
    ${MATH_HS})
    
# matrix module    
        
set(MATRIX_BASE_HS_
    ${INC}/matrix/internal/matrix_layout_internal.h
    ${INC}/matrix/internal/matrix_iter_internal.h
    ${INC}/matrix/matrix_shape.h
    ${INC}/matrix/matrix_layout.h
    ${INC}/matrix/matrix_fwd.h 
    ${INC}/matrix/matrix_meta.h
    ${INC}/matrix/matrix_shape.h
    ${INC}/matrix/matrix_concepts.h
    ${INC}/matrix/matrix_properties.h
    ${INC}/matrix/matrix_iter.h
    ${INC}/matrix/matrix_base.h)
    
set(MATRIX_OPS_HS_
    ${INC}/matrix/internal/matrix_fill_internal.h
    ${INC}/matrix/internal/matrix_copy_internal.h
    ${INC}/matrix/matrix_fill.h
    ${INC}/matrix/matrix_copy.h
    ${INC}/matrix/matrix_print.h)
    
set(MATRIX_VIEWS_HS_
    ${INC}/matrix/internal/matrix_colviews_internal.h
    ${INC}/matrix/internal/matrix_rowviews_internal.h
    ${INC}/matrix/internal/matrix_matviews_internal.h
    ${INC}/matrix/internal/matrix_asvec_internal.h
    ${INC}/matrix/matrix_subviews.h
    ${INC}/matrix/matrix_asvec.h)
    
set(MATRIX_CLASS_HS_
    ${INC}/matrix/regular_mat_base.h
    ${INC}/matrix/dense_matrix.h
    ${INC}/matrix/ref_matrix.h
    ${INC}/matrix/ref_block.h
    ${INC}/matrix/ref_grid.h
    ${INC}/matrix/step_vecs.h
    ${INC}/matrix/dense_mutable_view.h
    ${INC}/matrix/matrix_classes.h)
    
set(MATRIX_MANIP_HS_
    ${INC}/matrix/internal/matrix_transpose_internal.h
    ${INC}/matrix/matrix_transpose.h
    ${INC}/matrix/matrix_select.h)
    
set(MATRIX_HS
    ${MATRIX_BASE_HS_}
    ${MATRIX_OPS_HS_}
    ${MATRIX_CLASS_HS_}
    ${MATRIX_VIEWS_HS_}
    ${MATRIX_MANIP_HS_})
    
set(MATRIX_HS_EX
    ${CONFIG_HS}
    ${COMMON_HS}
    ${MATRIX_HS})
    
# matrix evaluation module
    
set(MATRIX_EVAL_HS_
    ${INC}/mateval/internal/ewise_eval_internal.h
    ${INC}/mateval/mateval_fwd.h
    ${INC}/mateval/fun_maps.h
    ${INC}/mateval/common_kernels.h
    ${INC}/mateval/vec_accessors.h
    ${INC}/mateval/multicol_accessors.h
    ${INC}/mateval/macc_policy.h
    ${INC}/mateval/ewise_eval.h)
    
set(MAP_EXPR_HS_
    ${INC}/mateval/internal/map_expr_internal.h
    ${INC}/mateval/map_accessors.h
    ${INC}/mateval/map_expr.h
    ${INC}/mateval/mat_arith.h
    ${INC}/mateval/mat_emath.h
    ${INC}/mateval/mat_cast.h
    ${INC}/mateval/mat_pred.h)
    
set(OTHER_EXPR_HS_
    ${INC}/mateval/repvec_expr.h
    ${INC}/mateval/subs_expr.h
    ${INC}/mateval/mat_zip.h)  
    
set(MATRIX_REDUC_HS_
    ${INC}/mateval/internal/mat_fold_internal.h
    ${INC}/mateval/internal/mat_reduce_internal.h
    ${INC}/mateval/internal/mat_allany_internal.h
    ${INC}/mateval/mat_fold.h
    ${INC}/mateval/mat_reduce.h
    ${INC}/mateval/mat_enorms.h
    ${INC}/mateval/mat_minmax.h
    ${INC}/mateval/mat_allany.h
    ${INC}/mateval/mat_compare.h)   
    
set(MATRIX_ALG_HS_
    ${INC}/mateval/internal/matrix_find_internal.h
    ${INC}/mateval/matrix_find.h
    ${INC}/mateval/matrix_sort.h)  
    
set(MATEVAL_HS
    ${MATRIX_EVAL_HS_}
    ${MAP_EXPR_HS_}
    ${OTHER_EXPR_HS_}
    ${MATRIX_REDUC_HS_}
    ${MATRIX_ALG_HS_}
    ${INC}/mateval/mateval.h)
    
set(MATEVAL_HS_EX
    ${CONFIG_HS}
    ${COMMON_HS}
    ${MATH_HS}
    ${MATRIX_HS}
    ${MATEVAL_HS})
    
# linear algebra
    
set(LINALG_BASE_HS_
    ${INC}/linalg/linalg_fwd.h
    ${INC}/linalg/internal/linalg_aux.h)
    
set(BLAS_HS_
    ${INC}/linalg/blas_l1.h
    ${INC}/linalg/blas_l2.h
    ${INC}/linalg/blas_l3.h
    ${INC}/linalg/blas.h)    
    
set(LAPACK_HS_
    ${INC}/linalg/lapack_fwd.h
    ${INC}/linalg/lapack_lu.h
    ${INC}/linalg/lapack_chol.h
    ${INC}/linalg/lapack_qr.h
    ${INC}/linalg/lapack_syev.h
    ${INC}/linalg/lapack_svd.h)
    
set(LINALG_HS
    ${LINALG_BASE_HS_}
    ${BLAS_HS_}
    ${LAPACK_HS_})
    
set(LINALG_HS_EX
    ${CONFIG_HS}
    ${COMMON_HS}
    ${MATRIX_HS}
    ${LINALG_HS})
    
#==========================================================
#
#    Unit testing targets
#
#==========================================================
    
# main entry

add_library(test_main STATIC test_main.cpp)

# common module

set(COMMON_MEM_TEST_HS
    ${COMMON_HS_EX})

add_executable(test_memory ${COMMON_MEM_TEST_HS} common/test_memory.cpp)
add_executable(test_blocks ${COMMON_MEM_TEST_HS} common/test_blocks.cpp)

set(LMAT_COMMON_TESTS
    test_memory
    test_blocks)

# math module

set(SSE_TEST_HS
    ${COMMON_HS_EX}
    ${SIMD_BASE_HS_}
    ${SSE_HS_})

add_executable(test_sse_packs  ${SSE_TEST_HS} math/test_sse_packs.cpp)
add_executable(test_sse_bpacks ${SSE_TEST_HS} math/test_sse_bpacks.cpp)
add_executable(test_sse_ops    ${SSE_TEST_HS} math/test_sse_ops.cpp)
add_executable(test_sse_round  ${SSE_TEST_HS} math/test_sse_round.cpp)
add_executable(test_sse_reduce ${SSE_TEST_HS} math/test_sse_reduce.cpp)
if (SVML_PATH)
add_executable(test_sse_math_svml ${SSE_TEST_HS} math/test_sse_math.cpp)
endif (SVML_PATH)

set(AVX_TEST_HS
    ${COMMON_HS_EX}
    ${SIMD_BASE_HS_}
    ${AVX_HS_})

if (ALLOW_AVX)
add_executable(test_avx_packs  ${AVX_TEST_HS} math/test_avx_packs.cpp)
add_executable(test_avx_bpacks ${AVX_TEST_HS} math/test_avx_bpacks.cpp)
add_executable(test_avx_ops    ${AVX_TEST_HS} math/test_avx_ops.cpp)
add_executable(test_avx_round  ${AVX_TEST_HS} math/test_avx_round.cpp)
add_executable(test_avx_reduce ${AVX_TEST_HS} math/test_avx_reduce.cpp)
if (SVML_PATH)
add_executable(test_avx_math_svml ${AVX_TEST_HS} math/test_avx_math.cpp)
endif (SVML_PATH)
endif (ALLOW_AVX)

set(LMAT_SSE_TESTS
    test_sse_packs
    test_sse_bpacks
    test_sse_ops
    test_sse_round
    test_sse_reduce)

if (ALLOW_AVX)
set(LMAT_AVX_TESTS
    test_avx_packs
    test_avx_bpacks
    test_avx_ops
    test_avx_round
    test_avx_reduce)
endif (ALLOW_AVX)

if (SVML_PATH)
set(LMAT_SSE_TESTS
    ${LMAT_SSE_TESTS}
    test_sse_math_svml)
if (ALLOW_AVX)
set(LMAT_AVX_TESTS
    ${LMAT_AVX_TESTS}
    test_avx_math_svml)
endif (ALLOW_AVX)
endif (SVML_PATH)

if (ALLOW_AVX)
set(LMAT_MATH_TESTS
    ${LMAT_SSE_TESTS}
    ${LMAT_AVX_TESTS})
else (ALLOW_AVX)
set(LMAT_MATH_TESTS
    ${LMAT_SSE_TESTS})
endif (ALLOW_AVX) 

# matrix module

set(MATCLASS_TEST_HS
    ${COMMON_HS_EX}
    ${MATRIX_BASE_HS_}
    ${MATRIX_CLASS_HS_})

add_executable(test_dense_mat ${MATCLASS_TEST_HS} matrix/test_dense_mat.cpp)
add_executable(test_dense_vec ${MATCLASS_TEST_HS} matrix/test_dense_vec.cpp)     
add_executable(test_ref_mat   ${MATCLASS_TEST_HS} matrix/test_ref_mat.cpp)
add_executable(test_ref_vec   ${MATCLASS_TEST_HS} matrix/test_ref_vec.cpp)
add_executable(test_ref_block ${MATCLASS_TEST_HS} matrix/test_ref_block.cpp)
add_executable(test_ref_grid  ${MATCLASS_TEST_HS} matrix/test_ref_grid.cpp)
add_executable(test_step_vec  ${MATCLASS_TEST_HS} matrix/test_step_vec.cpp)

set(MATOPS_TEST_HS
    ${COMMON_HS_EX}
    ${MATRIX_BASE_HS_}
    ${MATRIX_CLASS_HS_}
    ${MATRIX_OPS_HS_})

add_executable(test_mat_props  ${MATOPS_TEST_HS} matrix/test_mat_props.cpp)
add_executable(test_mat_iter   ${MATOPS_TEST_HS} matrix/test_mat_iter.cpp)
add_executable(test_mat_fill   ${MATOPS_TEST_HS} matrix/test_mat_fill.cpp)
add_executable(test_mat_copy   ${MATOPS_TEST_HS} matrix/test_mat_copy.cpp)
add_executable(test_dense_eval ${MATOPS_TEST_HS}  matrix/test_dense_eval.cpp)

set(MATVIEWS_TEST_HS
    ${COMMON_HS_EX}
    ${MATRIX_BASE_HS_}
    ${MATRIX_CLASS_HS_}
    ${MATRIX_VIEWS_HS_})

add_executable(test_mat_vecviews ${MATVIEWS_TEST_HS} matrix/test_mat_vecviews.cpp)
add_executable(test_mat_matviews ${MATVIEWS_TEST_HS} matrix/test_mat_matviews.cpp)
add_executable(test_mat_asvec ${MATVIEWS_TEST_HS} matrix/test_mat_asvec.cpp) 

set(MATMANIP_TEST_HS
    ${COMMON_HS}
    ${MATRIX_BASE_HS_}
    ${MATRIX_CLASS_HS_}
    ${MATRIX_MANIP_HS_})
     
add_executable(test_mat_select ${MATMANIP_TEST_HS} matrix/test_mat_select.cpp)
add_executable(test_direct_trans ${MATMANIP_TEST_HS} matrix/test_direct_trans.cpp)    
add_executable(test_transpose_expr ${MATMANIP_TEST_HS} matrix/test_transpose_expr.cpp) 

set(LMAT_MATRIX_TESTS
    test_dense_mat
	test_dense_vec
	test_ref_mat
	test_ref_vec
	test_ref_block
	test_ref_grid
	test_step_vec
	test_mat_props
	test_mat_iter
	test_mat_fill
	test_mat_copy
	test_dense_eval
	test_mat_vecviews
	test_mat_matviews
	test_mat_asvec
	test_mat_select
	test_direct_trans
	test_transpose_expr)

# matrix evaluation module

set(MATEVAL_TEST_HS
    ${MATRIX_HS}
    ${MATRIX_EVAL_HS_})
    
add_executable(test_linear_ewise ${MATEVAL_TEST_HS} mateval/test_linear_ewise.cpp)    
add_executable(test_percol_ewise ${MATEVAL_TEST_HS} mateval/test_percol_ewise.cpp)
add_executable(test_ewise_accum ${MATEVAL_TEST_HS}  mateval/test_ewise_accum.cpp)

set(MAPEXPR_TEST_HS
    ${MATRIX_HS}
    ${MATRIX_EVAL_HS_}
    ${MAP_EXPR_HS_})
    
add_executable(test_map_expr ${MAPEXPR_TEST_HS} mateval/test_map_expr.cpp)
add_executable(test_mat_arith ${MAPEXPR_TEST_HS} mateval/test_mat_arith.cpp)
add_executable(test_mat_emath ${MAPEXPR_TEST_HS} mateval/test_mat_emath.cpp)
add_executable(test_mat_cast ${MAPEXPR_TEST_HS} mateval/test_mat_cast.cpp)
add_executable(test_mat_pred ${MAPEXPR_TEST_HS} mateval/test_mat_pred.cpp)

set(OTHEREXPR_TEST_HS
    ${MATRIX_HS}
    ${MATRIX_EVAL_HS_}
    ${OTHER_EXPR_HS_})

add_executable(test_repvecs ${OTHEREXPR_TEST_HS} mateval/test_repvecs.cpp)
add_executable(test_subs_expr ${OTHEREXPR_TEST_HS} mateval/test_subs_expr.cpp)
add_executable(test_mat_zip ${OTHEREXPR_TEST_HS} mateval/test_mat_zip.cpp)

set(MATREDUC_TEST_HS
    ${MATRIX_HS}
    ${MATRIX_EVAL_HS_}
    ${MAP_EXPR_HS_}
    ${MATRIX_REDUC_HS_})

add_executable(test_full_reduce    ${MATREDUC_TEST_HS} mateval/test_full_reduce.cpp)
add_executable(test_colwise_reduce ${MATREDUC_TEST_HS} mateval/test_colwise_reduce.cpp)
add_executable(test_rowwise_reduce ${MATREDUC_TEST_HS} mateval/test_rowwise_reduce.cpp)
add_executable(test_more_reduce ${MATREDUC_TEST_HS} mateval/test_more_reduce.cpp)
add_executable(test_mat_allany ${MATREDUC_TEST_HS} mateval/test_mat_allany.cpp)
add_executable(test_mat_compare ${MATREDUC_TEST_HS} mateval/test_mat_compare.cpp)

set(MATALG_TEST_HS
    ${MATRIX_HS}
    ${MATRIX_EVAL_HS_}
    ${MAP_EXPR_HS_}
    ${OTHER_EXPR_HS_}
    ${MATRIX_ALG_HS_})
    
add_executable(test_mat_find ${MATALG_TEST_HS} mateval/test_mat_find.cpp)
add_executable(test_mat_sort ${MATALG_TEST_HS} mateval/test_mat_sort.cpp)

set(MATEVAL_ALL_TEST_HS
    ${MATEVAL_HS_EX})

add_executable(test_cpd_ewise ${MATEVAL_ALL_TEST_HS} mateval/test_cpd_ewise.cpp)

set(LMAT_MATEVAL_TESTS
    test_linear_ewise
	test_percol_ewise
	test_ewise_accum
	test_map_expr
	test_mat_arith
	test_mat_emath
	test_mat_cast
	test_mat_pred
	test_repvecs
	test_subs_expr
	test_mat_zip
	test_full_reduce
	test_colwise_reduce
	test_rowwise_reduce
	test_more_reduce
	test_mat_allany
	test_mat_compare
	test_mat_find
	test_mat_sort
	test_cpd_ewise)


# linear algebra module

set(BLAS_TEST_HS
    ${MATRIX_HS}
    ${BLAS_HS_})

add_executable(test_blas_l1 ${BLAS_TEST_HS} linalg/test_blas_l1.cpp)
add_executable(test_blas_l2 ${BLAS_TEST_HS} linalg/test_blas_l2.cpp)
add_executable(test_blas_l3 ${BLAS_TEST_HS} linalg/test_blas_l3.cpp)

set(LAPACK_TEST_HS
    ${MATRIX_HS}
    ${BLAS_HS_}
    ${LAPACK_HS_})
    
add_executable(test_lapack_lu ${LAPACK_TEST_HS} linalg/test_lapack_lu.cpp)
add_executable(test_lapack_chol ${LAPACK_TEST_HS} linalg/test_lapack_chol.cpp)
add_executable(test_lapack_qr ${LAPACK_TEST_HS} linalg/test_lapack_qr.cpp)
add_executable(test_lapack_syev ${LAPACK_TEST_HS} linalg/test_lapack_syev.cpp)
add_executable(test_lapack_svd ${LAPACK_TEST_HS} linalg/test_lapack_svd.cpp)

set(LMAT_LINALG_TESTS
    test_blas_l1
    test_blas_l2
    test_blas_l3
    test_lapack_lu
    test_lapack_chol
    test_lapack_qr
    test_lapack_syev
    test_lapack_svd)
    
# all

set(LMAT_ALL_TESTS
    ${LMAT_COMMON_TESTS}
    ${LMAT_MATH_TESTS}
    ${LMAT_MATRIX_TESTS}
    ${LMAT_MATEVAL_TESTS}
    ${LMAT_LINALG_TESTS}
)


#==========================================================
#
#    Special treatment
#
#==========================================================

# suppress bogus free-nonheap-object complaint under GCC

if (${CMAKE_CXX_COMPILER_ID} MATCHES "GNU")
set_target_properties(test_dense_eval PROPERTIES COMPILE_FLAGS "-Wno-free-nonheap-object")
endif (${CMAKE_CXX_COMPILER_ID} MATCHES "GNU")
    
# Link to SVML

if (SVML_PATH)

set(TESTS_USING_SVML
    test_sse_math_svml
    test_mat_emath
    test_cpd_ewise)
    
if (ALLOW_AVX)
set(TESTS_USING_SVML
    ${TESTS_USING_SVML}
    test_avx_math_svml)
endif (ALLOW_AVX)

foreach (tname ${TESTS_USING_SVML}) 
set_target_properties(${tname}
	PROPERTIES
	COMPILE_FLAGS "-DLMAT_USE_INTEL_SVML"
)
if (MSVC)
target_link_libraries(${tname} ${SVML} ${LIBIRC})
else (MSVC)
target_link_libraries(${tname} ${SVML})
endif (MSVC)
endforeach (tname)
    
endif (SVML_PATH)
    
# Link to MKL

foreach(tname ${LMAT_LINALG_TESTS})
    target_link_libraries(${tname} ${MKL_LIBRARY})
endforeach(tname)


#==========================================================
#
#    treatments to all tests
#
#==========================================================

# Link to test_main
	
foreach(tname ${LMAT_ALL_TESTS})
	target_link_libraries(${tname} test_main)
endforeach(tname)	


# Add tests

foreach(tname ${LMAT_ALL_TESTS})
	add_test(NAME ${tname} COMMAND ${tname})
endforeach(tname)


