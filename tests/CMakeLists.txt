# The CMakeLists file for Light-Matrix (Testing)

cmake_minimum_required(VERSION 2.8)
project(LightMatrix)

enable_testing()

include_directories(..)
include_directories($ENV{LIGHT_TEST_HOME})

include("../cmake_modules/CompilerConfig.cmake")

#==========================================================
#
#    third-party library dependencies
#
#==========================================================

set(CMAKE_MODULE_PATH "../cmake_modules" ${CMAKE_MODULE_PATH})

# Intel library 

find_package(ICCLib)
if (ICCLIB_FOUND)
message(STATUS "[LMAT] ICC Library found: ${ICCLIB_PATH}")
else (ICCLIB_FOUND)
message(STATUS "[LMAT] ICC Library not found")
endif (ICCLIB_FOUND)

set(SVML_FOUND ICCLIB_FOUND)

# AMD LibM

find_package(LIBM)
if (LIBM_FOUND)
message(STATUS "[LMAT] AMD LibM found: ${LIBM_PATH}")
else (LIBM_FOUND)
message(STATUS "[LMAT] AMD LibM not found")
endif (LIBM_FOUND)

# Intel MKL

find_package(MKL)
if (MKL_FOUND)
message(STATUS "[LMAT] Intel MKL found: ${MKLROOT_PATH}")
else (MKL_FOUND)
message(STATUS "[LMAT] Intel MKL not found")
endif (MKL_FOUND)

set(BLAS_FOUND MKL_FOUND)
set(LAPACK_FOUND MKL_FOUND)


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
    ${INC}/common/preprocess_base.h
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
    
# simd module

set(SIMD_BASE_HS_
    ${INC}/simd/internal/numrepr_format.h
    ${INC}/simd/simd_arch.h
    ${INC}/simd/simd_base.h
    ${INC}/simd/simd_debug.h)
    
set(SSE_HS_
    ${INC}/simd/internal/sse_helpers.h
    ${INC}/simd/internal/sse_testz_impl.h
    ${INC}/simd/internal/sse2_round_impl.h
    ${INC}/simd/internal/sse_fpclass_impl.h
    ${INC}/simd/sse_packs.h
    ${INC}/simd/sse_bpacks.h
    ${INC}/simd/sse_arith.h
    ${INC}/simd/sse_pred.h
    ${INC}/simd/sse_reduce.h
    ${INC}/simd/sse.h)
    
set(AVX_HS_
    ${INC}/simd/internal/avx_helpers.h
    ${INC}/simd/avx_packs.h
    ${INC}/simd/avx_bpacks.h
    ${INC}/simd/avx_arith.h
    ${INC}/simd/avx_pred.h
    ${INC}/simd/avx_reduce.h
    ${INC}/simd/avx.h) 
    
set(SIMD_HS
    ${SIMD_BASE_HS_}
    ${SSE_HS_}
    ${AVX_HS_})    
    
set(SIMD_HS_EX
    ${CONFIG_HS}
    ${COMMON_HS}
    ${SIMD_HS})    
    
# math module
    
set(MATH_BASE_HS_
    ${INC}/math/internal/cmath_win32.h
    ${INC}/math/internal/svml_import.h
    ${INC}/math/internal/libm_simd_import.h
    ${INC}/math/math_base.h
    ${INC}/math/math_constants.h
    ${INC}/math/math.h
    ${INC}/math/simd_math.h)
    
set(MATH_SPECIAL_HS_
    ${INC}/math/internal/norminv_impl.h
    ${INC}/math/math_special.h)
    
set(MATH_FUNCTOR_HS_
    ${INC}/math/fun_tags.h
    ${INC}/math/functor_base.h
    ${INC}/math/basic_functors.h
    ${INC}/math/math_functors.h
    ${INC}/math/special_functors.h) 
    
set(MATH_HS
    ${MATH_BASE_HS_}
    ${MATH_SPECIAL_HS_}
    ${MATH_FUNCTOR_HS_})    
    
set(MATH_HS_EX
    ${CONFIG_HS}
    ${COMMON_HS}
    ${SIMD_HS}
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
    ${INC}/matrix/matrix_init.h
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
    ${INC}/matrix/matrix_classes.h
    ${INC}/matrix/matrix_inspect.h)
    
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
    ${INC}/mateval/common_kernels.h
    ${INC}/mateval/vec_accessors.h
    ${INC}/mateval/multicol_accessors.h
    ${INC}/mateval/macc_policy.h
    ${INC}/mateval/ewise_eval.h)
    
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
    ${INC}/mateval/matrix_sort.h
    ${INC}/mateval/matrix_ordstats.h)  
    
set(MATEVAL_HS
    ${MATRIX_EVAL_HS_}
    ${MATRIX_REDUC_HS_}
    ${MATRIX_ALG_HS_}
    ${INC}/mateval/mateval.h)
    
set(MATEVAL_HS_EX
    ${CONFIG_HS}
    ${COMMON_HS}
    ${SIMD_HS}
    ${MATH_HS}
    ${MATRIX_HS}
    ${MATEVAL_HS})
    
# matrix expression 

set(MAP_EXPR_HS_
    ${INC}/matexpr/internal/map_expr_internal.h
    ${INC}/matexpr/map_accessors.h
    ${INC}/matexpr/map_expr.h
    ${INC}/matexpr/map_expr_inspect.h
    ${INC}/matexpr/matfun_base.h
    ${INC}/matexpr/mat_arith.h
    ${INC}/matexpr/mat_emath.h
    ${INC}/matexpr/mat_special.h
    ${INC}/matexpr/mat_cast.h
    ${INC}/matexpr/mat_pred.h)
    
set(OTHER_EXPR_HS_
    ${INC}/matexpr/repvec_expr.h
    ${INC}/matexpr/subs_expr.h
    ${INC}/matexpr/mat_zip.h) 
    
set(MATEXPR_HS
    ${MAP_EXPR_HS_}
    ${OTHER_EXPR_HS_})

set(MATEXPR_HS_EX
    ${MATEVAL_HS_EX}
    ${MATEXPR_HS})
               
    
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
    
    
# random 

set(PRNG_HS_
    ${INC}/random/internal/rand_stream_internal.h
    ${INC}/random/internal/sfmt_params.h
    ${INC}/random/rand_stream.h
    ${INC}/random/stream_tracker.h
    ${INC}/random/sfmt.h)
    
set(DISTR_HS_
    ${INC}/random/distr_fwd.h
    ${INC}/random/internal/uniform_real_internal.h
    ${INC}/random/internal/normal_distr_internal.h
    ${INC}/random/internal/gamma_distr_internal.h
    ${INC}/random/uniform_int_distr.h
    ${INC}/random/sample_wor.h
    ${INC}/random/bernoulli_distr.h
    ${INC}/random/binomial_distr.h
    ${INC}/random/geometric_distr.h
    ${INC}/random/poisson_distr.h
    ${INC}/random/discrete_distr.h
    ${INC}/random/uniform_real_distr.h
    ${INC}/random/exponential_distr.h
    ${INC}/random/normal_distr.h
    ${INC}/random/gamma_distr.h
)
        
set(RANDOM_HS
    ${PRNG_HS_}
    ${DISTR_HS_}
    ${INC}/random/rand_accessors.h
    ${INC}/random/rand_expr.h)
    
set(RANDOM_HS_EX
    ${MATEXPR_HS_EX}
    ${RANDOM_HS})
        
    
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

# simd module

set(SSE_TEST_HS
    ${COMMON_HS_EX}
    ${SIMD_BASE_HS_}
    ${SSE_HS_})

add_executable(test_sse_packs  ${SSE_TEST_HS} simd/test_sse_packs.cpp)
add_executable(test_sse_bpacks ${SSE_TEST_HS} simd/test_sse_bpacks.cpp)
add_executable(test_sse_arith  ${SSE_TEST_HS} simd/test_sse_arith.cpp)
add_executable(test_sse_pred   ${SSE_TEST_HS} simd/test_sse_pred.cpp)
add_executable(test_sse_round  ${SSE_TEST_HS} simd/test_sse_round.cpp)
add_executable(test_sse_reduce ${SSE_TEST_HS} simd/test_sse_reduce.cpp)

set(AVX_TEST_HS
    ${COMMON_HS_EX}
    ${SIMD_BASE_HS_}
    ${AVX_HS_})

if (ALLOW_AVX)
add_executable(test_avx_packs  ${AVX_TEST_HS} simd/test_avx_packs.cpp)
add_executable(test_avx_bpacks ${AVX_TEST_HS} simd/test_avx_bpacks.cpp)
add_executable(test_avx_arith  ${SSE_TEST_HS} simd/test_avx_arith.cpp)
add_executable(test_avx_pred   ${SSE_TEST_HS} simd/test_avx_pred.cpp)
add_executable(test_avx_round  ${AVX_TEST_HS} simd/test_avx_round.cpp)
add_executable(test_avx_reduce ${AVX_TEST_HS} simd/test_avx_reduce.cpp)
endif (ALLOW_AVX)

set(LMAT_SSE_TESTS
    test_sse_packs
    test_sse_bpacks
    test_sse_arith
    test_sse_pred
    test_sse_round
    test_sse_reduce)

if (ALLOW_AVX)
set(LMAT_AVX_TESTS
    test_avx_packs
    test_avx_bpacks
    test_avx_arith
    test_avx_pred
    test_avx_round
    test_avx_reduce)
endif (ALLOW_AVX)

if (ALLOW_AVX)
set(LMAT_SIMD_TESTS
    ${LMAT_SSE_TESTS}
    ${LMAT_AVX_TESTS})
else (ALLOW_AVX)
set(LMAT_SIMD_TESTS
    ${LMAT_SSE_TESTS})
endif (ALLOW_AVX) 

# math module

if (SVML_FOUND)

add_executable(test_simd_math_svml ${MATH_HS_EX} math/test_simd_math.cpp)
add_executable(test_simd_special ${MATH_HS_EX} math/test_simd_special.cpp)

set(LMAT_MATH_TESTS
    test_simd_math_svml
    test_simd_special)

endif (SVML_FOUND)

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
	test_transpose_expr
	)

# matrix evaluation module

set(MATEVAL_TEST_HS
    ${MATRIX_HS}
    ${MATH_HS}
    ${MATRIX_EVAL_HS_})
    
add_executable(test_linear_ewise ${MATEVAL_TEST_HS} mateval/test_linear_ewise.cpp)    
add_executable(test_percol_ewise ${MATEVAL_TEST_HS} mateval/test_percol_ewise.cpp)
add_executable(test_ewise_accum ${MATEVAL_TEST_HS}  mateval/test_ewise_accum.cpp)

set(MATREDUC_TEST_HS
    ${MATRIX_HS}
    ${MATH_HS}
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
add_executable(test_mat_ordstat ${MATALG_TEST_HS} mateval/test_mat_ordstat.cpp)

set(LMAT_MATEVAL_TESTS
    test_linear_ewise
	test_percol_ewise
	test_ewise_accum
	test_full_reduce
	test_colwise_reduce
	test_rowwise_reduce
	test_more_reduce
	test_mat_allany
	test_mat_compare
	test_mat_find
	test_mat_sort
	test_mat_ordstat
	)


# matrix expressions

set(MAPEXPR_TEST_HS
    ${MATRIX_HS}
    ${MATH_HS}
    ${MATRIX_EVAL_HS_}
    ${MAP_EXPR_HS_})
    
add_executable(test_map_expr ${MAPEXPR_TEST_HS} matexpr/test_map_expr.cpp)
add_executable(test_mat_arith ${MAPEXPR_TEST_HS} matexpr/test_mat_arith.cpp)
add_executable(test_mat_emath ${MAPEXPR_TEST_HS} matexpr/test_mat_emath.cpp)
add_executable(test_mat_special ${MAPEXPR_TEST_HS} matexpr/test_mat_special.cpp)
add_executable(test_mat_cast ${MAPEXPR_TEST_HS} matexpr/test_mat_cast.cpp)
add_executable(test_mat_pred ${MAPEXPR_TEST_HS} matexpr/test_mat_pred.cpp)

set(OTHEREXPR_TEST_HS
    ${MATRIX_HS}
    ${MATH_HS}
    ${MATRIX_EVAL_HS_}
    ${OTHER_EXPR_HS_})

add_executable(test_repvecs ${OTHEREXPR_TEST_HS} matexpr/test_repvecs.cpp)
add_executable(test_subs_expr ${OTHEREXPR_TEST_HS} matexpr/test_subs_expr.cpp)
add_executable(test_mat_zip ${OTHEREXPR_TEST_HS} matexpr/test_mat_zip.cpp)

add_executable(test_cpd_ewise ${MATEXPR_HS_EX} matexpr/test_cpd_ewise.cpp)

set(LMAT_MATEXPR_TESTS
	test_map_expr
	test_mat_arith
	test_mat_emath
	test_mat_special
	test_mat_cast
	test_mat_pred
	test_repvecs
	test_subs_expr
	test_mat_zip
	test_cpd_ewise
	)


# linear algebra module

if (BLAS_FOUND)

set(BLAS_TEST_HS
    ${MATRIX_HS}
    ${BLAS_HS_})

add_executable(test_blas_l1 ${BLAS_TEST_HS} linalg/test_blas_l1.cpp)
add_executable(test_blas_l2 ${BLAS_TEST_HS} linalg/test_blas_l2.cpp)
add_executable(test_blas_l3 ${BLAS_TEST_HS} linalg/test_blas_l3.cpp)

set(LMAT_BLAS_TESTS
    test_blas_l1
    test_blas_l2
    test_blas_l3
)

else (BLAS_FOUND)
set(LMAT_BLAS_TESTS)
endif(BLAS_FOUND)

if (LAPACK_FOUND)

set(LAPACK_TEST_HS
    ${MATRIX_HS}
    ${BLAS_HS_}
    ${LAPACK_HS_})
    
add_executable(test_lapack_lu ${LAPACK_TEST_HS} linalg/test_lapack_lu.cpp)
add_executable(test_lapack_chol ${LAPACK_TEST_HS} linalg/test_lapack_chol.cpp)
add_executable(test_lapack_qr ${LAPACK_TEST_HS} linalg/test_lapack_qr.cpp)
add_executable(test_lapack_syev ${LAPACK_TEST_HS} linalg/test_lapack_syev.cpp)
add_executable(test_lapack_svd ${LAPACK_TEST_HS} linalg/test_lapack_svd.cpp)

set(LMAT_LAPACK_TESTS
    test_lapack_lu
    test_lapack_chol
    test_lapack_qr
    test_lapack_syev
    test_lapack_svd)

elseif(LAPACK_FOUND)
set(LMAT_LAPACK_TESTS)
endif (LAPACK_FOUND)

set(LMAT_LINALG_TESTS
    ${LMAT_BLAS_TESTS}
    ${LMAT_LAPACK_TESTS})
    
# random module

set(PRNG_TEST_HS
    ${COMMON_HS}
    ${MATH_HS_EX}
    ${PRNG_HS_})
    
add_executable(test_stracker ${PRNG_TEST_HS} random/test_stracker.cpp)
add_executable(test_sfmt ${PRNG_TEST_HS} random/test_sfmt.cpp)
  
# SFMT verification files
configure_file(random/data/sfmt/sfmt.000607.txt ${CMAKE_CURRENT_BINARY_DIR}/data/sfmt/sfmt.000607.txt COPYONLY)
configure_file(random/data/sfmt/sfmt.001279.txt ${CMAKE_CURRENT_BINARY_DIR}/data/sfmt/sfmt.001279.txt COPYONLY)
configure_file(random/data/sfmt/sfmt.002281.txt ${CMAKE_CURRENT_BINARY_DIR}/data/sfmt/sfmt.002281.txt COPYONLY)
configure_file(random/data/sfmt/sfmt.004253.txt ${CMAKE_CURRENT_BINARY_DIR}/data/sfmt/sfmt.004253.txt COPYONLY)
configure_file(random/data/sfmt/sfmt.011213.txt ${CMAKE_CURRENT_BINARY_DIR}/data/sfmt/sfmt.011213.txt COPYONLY)
configure_file(random/data/sfmt/sfmt.019937.txt ${CMAKE_CURRENT_BINARY_DIR}/data/sfmt/sfmt.019937.txt COPYONLY)
configure_file(random/data/sfmt/sfmt.044497.txt ${CMAKE_CURRENT_BINARY_DIR}/data/sfmt/sfmt.044497.txt COPYONLY)
configure_file(random/data/sfmt/sfmt.086243.txt ${CMAKE_CURRENT_BINARY_DIR}/data/sfmt/sfmt.086243.txt COPYONLY)
configure_file(random/data/sfmt/sfmt.132049.txt ${CMAKE_CURRENT_BINARY_DIR}/data/sfmt/sfmt.132049.txt COPYONLY)
configure_file(random/data/sfmt/sfmt.216091.txt ${CMAKE_CURRENT_BINARY_DIR}/data/sfmt/sfmt.216091.txt COPYONLY)

set(DISTR_TEST_HS
    random/distr_test_base.h
    ${COMMON_HS}
    ${MATH_HS_EX}
    ${PRNG_HS_}
    ${DISTR_HS_})

add_executable(test_uniform_int ${DISTR_TEST_HS} random/test_uniform_int.cpp)
add_executable(test_sample_wor ${DISTR_TEST_HS} random/test_sample_wor.cpp)
add_executable(test_bernoulli ${DISTR_TEST_HS} random/test_bernoulli.cpp)
add_executable(test_binomial ${DISTR_TEST_HS} random/test_binomial.cpp)
add_executable(test_geometric ${DISTR_TEST_HS} random/test_geometric.cpp)
add_executable(test_poissond ${DISTR_TEST_HS} random/test_poissond.cpp)
add_executable(test_discreted ${DISTR_TEST_HS} random/test_discreted.cpp)
  
add_executable(test_uniform_real ${DISTR_TEST_HS} random/test_uniform_real.cpp)
add_executable(test_exponential ${DISTR_TEST_HS} random/test_exponential.cpp)
add_executable(test_normal ${DISTR_TEST_HS} random/test_normal.cpp)
add_executable(test_gammad ${DISTR_TEST_HS} random/test_gammad.cpp)

add_executable(test_rand_expr ${RANDOM_HS_EX} random/test_rand_expr.cpp)
     
set(LMAT_RANDOM_TESTS
    test_stracker
    test_sfmt
    test_uniform_int
    test_sample_wor
    test_bernoulli
    test_binomial
    test_geometric
    test_poissond
    test_discreted
    test_uniform_real
    test_exponential
    test_normal
    test_gammad
    test_rand_expr)        

# all

set(LMAT_ALL_TESTS
    ${LMAT_COMMON_TESTS}
    ${LMAT_SIMD_TESTS}
    ${LMAT_MATH_TESTS}
    ${LMAT_MATRIX_TESTS}
    ${LMAT_MATEVAL_TESTS}
    ${LMAT_MATEXPR_TESTS}
    ${LMAT_LINALG_TESTS}
    ${LMAT_RANDOM_TESTS}
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

if (SVML_FOUND)

set(TESTS_USING_SVML
    test_simd_math_svml
    test_simd_special
    test_mat_emath
    test_mat_special
    test_cpd_ewise
    test_exponential
    test_normal
    test_rand_expr)

foreach (tname ${TESTS_USING_SVML}) 
    set_target_properties(${tname}
	    PROPERTIES
	    COMPILE_FLAGS "-DLMAT_USE_INTEL_SVML")
    target_link_libraries(${tname} ${SVML_LIBRARY})
endforeach (tname)
    
endif (SVML_FOUND)
    
# Link to MKL

if (MKL_FOUND)

foreach(tname ${LMAT_LINALG_TESTS})
    target_link_libraries(${tname} ${MKL_LIBRARY})
endforeach(tname)

endif (MKL_FOUND)

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


