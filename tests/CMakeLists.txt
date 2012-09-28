# The CMakeLists file for Light-Matrix (Testing)

cmake_minimum_required(VERSION 2.8)
project(LightMatrix)

enable_testing()

include_directories(..)
include_directories($ENV{LIGHT_TEST_HOME})

# Compiler configuration

set(CMAKE_BUILD_TYPE "Release")

if (MSVC)
	set(LANG_FLAGS "/arch:SSE2 /EHsc")
	set(WARNING_FLAGS "/W4")
else (MSVC)
	set(LANG_FLAGS "-std=c++0x -pedantic -march=native ")
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
    ${INC}/common/metaprog_base.h
    ${INC}/common/mask_type.h
    ${INC}/common/type_traits.h
    ${INC}/common/arg_check.h
    ${INC}/common/range.h
    ${INC}/common/basic_defs.h    
    ${INC}/common/memory.h
    ${INC}/common/block.h
    ${INC}/common/expr_base.h)
    
set(MATRIX_BASE_HS
    ${INC}/matrix/matrix_fwd.h 
    ${INC}/matrix/matrix_meta.h
    ${INC}/matrix/matrix_shape.h
    ${INC}/matrix/matrix_concepts.h
    ${INC}/matrix/matrix_properties.h
    ${INC}/matrix/bits/matrix_fill_internal.h
    ${INC}/matrix/bits/matrix_copy_internal.h
    ${INC}/matrix/bits/matrix_compare_internal.h
    ${INC}/matrix/bits/mat_transpose_impl.h
    ${INC}/matrix/matrix_fill.h
    ${INC}/matrix/matrix_copy.h
    ${INC}/matrix/matrix_compare.h
    ${INC}/matrix/matrix_print.h
    ${INC}/matrix/matrix_transpose.h
    ${INC}/matrix/matrix_subviews.h
    ${INC}/matrix/matrix_ewise_expr.h
    ${INC}/matrix/matrix_base.h)
    
set(MATRIX_CLASS_HS
    ${INC}/matrix/dense_matrix.h
    ${INC}/matrix/ref_matrix.h
    ${INC}/matrix/ref_matrix_ex.h
    ${INC}/matrix/const_matrix.h
    ${INC}/matrix/dense_mutable_view.h
    ${INC}/matrix/matrix_classes.h)
    
set(MATRIX_EVAL_BASE_HS
    ${INC}/math/math_base.h
    ${INC}/math/functor_base.h
    ${INC}/matrix/matrix_visitors.h
    ${INC}/matrix/matrix_visit_eval.h
    ${INC}/matrix/matrix_ewise_eval.h)   
    
set(MATRIX_EWISE_HS
    ${INC}/math/arith_functors.h
    ${INC}/math/emath_functors.h
    ${INC}/math/comp_functors.h
    ${INC}/math/mask_functors.h
    ${INC}/matrix/matrix_cast.h
    ${INC}/matrix/matrix_arith.h
    ${INC}/matrix/matrix_emath.h
    ${INC}/matrix/matrix_ecomp.h
    ${INC}/matrix/matrix_elogical.h)
    
set(MATRIX_REPEAT_HS
    ${INC}/matrix/bits/repeat_vecs_internal.h
    ${INC}/matrix/matrix_repeat_expr.h
    ${INC}/matrix/matrix_repeat_eval.h)   
    
set(MATRIX_REDUC_HS
    ${INC}/math/reduction_functors.h
    ${INC}/matrix/bits/matrix_reduce_internal.h
    ${INC}/matrix/bits/partial_reduce_internal.h
    ${INC}/matrix/matrix_reduce.h
    ${INC}/matrix/partial_reduce.h)    
    
set(BLAS_HS
    ${INC}/linalg/linalg_base.h
    ${INC}/linalg/blas_extern.h
    ${INC}/linalg/bits/matrix_blasl1_internal.h
    ${INC}/linalg/matrix_blasl1.h)    
              
set(MATRIX_HS
    ${COMMON_HS}
    ${MATRIX_BASE_HS}
    ${MATRIX_CLASS_HS})   
    
set(MATRIX_EVAL_HS
    ${COMMON_HS}
    ${MATRIX_BASE_HS}
    ${MATRIX_CLASS_HS}
    ${MATRIX_EVAL_BASE_HS}
    ${MATRIX_EWISE_HS}
    ${MATRIX_REPEAT_HS}
    ${MATRIX_REDUC_HS})    
    
# Executables

add_library(test_main STATIC test_main.cpp)

add_executable(test_blocks ${COMMON_HS} test_blocks.cpp)

set(DENSE_MAT_HS
    ${COMMON_HS}
    ${MATRIX_BASE_HS}
    ${INC}/matrix/dense_matrix.h)

add_executable(test_dense_mat ${DENSE_MAT_HS} test_dense_mat.cpp)
add_executable(test_dense_col ${DENSE_MAT_HS} test_dense_col.cpp)
add_executable(test_dense_row ${DENSE_MAT_HS} test_dense_row.cpp)  

set(CONST_MAT_HS
    ${COMMON_HS}
    ${MATRIX_BASE_HS}
    ${INC}/matrix/const_matrix.h
    ${INC}/matrix/dense_matrix.h)

add_executable(test_const_mat ${CONST_MAT_HS} test_const_mat.cpp)

set(REF_MAT_HS
    ${COMMON_HS}
    ${MATRIX_BASE_HS}
    ${INC}/matrix/ref_matrix.h
    ${INC}/matrix/ref_matrix_ex.h)
    
add_executable(test_ref_mat ${REF_MAT_HS} test_ref_mat.cpp)
add_executable(test_ref_col ${REF_MAT_HS} test_ref_col.cpp)
add_executable(test_ref_row ${REF_MAT_HS} test_ref_row.cpp)
add_executable(test_ref_mat_ex ${REF_MAT_HS} test_ref_mat_ex.cpp)

add_executable(test_mat_copy ${REF_MAT_HS} test_mat_copy.cpp)
add_executable(test_mat_compare ${REF_MAT_HS} test_mat_compare.cpp)

add_executable(test_dense_eval ${MATRIX_HS} test_dense_eval.cpp)
add_executable(test_mat_subviews ${MATRIX_HS} test_mat_subviews.cpp)
add_executable(test_mat_transpose ${MATRIX_HS} test_mat_transpose.cpp)

set(BASIC_MAT_EVAL_HS
    ${MATRIX_HS}
    ${MATRIX_EVAL_BASE_HS})
    
add_executable(test_dense_veval ${BASIC_MAT_EVAL_HS} test_dense_veval.cpp)

set(EWISE_MAT_EVAL_HS
    ${MATRIX_HS}
    ${MATRIX_EVAL_BASE_HS}
    ${MATRIX_EWISE_HS})

add_executable(test_dense_cast ${EWISE_MAT_EVAL_HS} test_dense_cast.cpp)   
add_executable(test_mat_arith ${EWISE_MAT_EVAL_HS} test_mat_arith.cpp)
add_executable(test_mat_emath ${EWISE_MAT_EVAL_HS} test_mat_emath.cpp)
add_executable(test_mat_ecomp ${EWISE_MAT_EVAL_HS} test_mat_ecomp.cpp)
add_executable(test_mat_elogical ${EWISE_MAT_EVAL_HS} test_mat_elogical.cpp)

set(MAT_REPEAT_HS
    ${MATRIX_HS}
    ${MATRIX_EVAL_BASE_HS}
    ${MATRIX_REPEAT_HS})
    
add_executable(test_vec_rep ${MAT_REPEAT_HS} test_vec_rep.cpp)

add_executable(test_cpd_ewise ${MATRIX_EVAL_HS} test_cpd_ewise.cpp)

set(MAT_REDUC_HS
    ${MATRIX_HS}
    ${MATRIX_EVAL_BASE_HS}
    ${MATRIX_REDUC_HS})
    
add_executable(test_mat_reduce ${MAT_REDUC_HS} test_mat_reduce.cpp) 
add_executable(test_partial_reduce ${MAT_REDUC_HS} test_partial_reduce.cpp) 

# Link to test_main

set(LMAT_TEST_NAMES
	test_blocks
	test_dense_mat
	test_dense_col
	test_dense_row
	test_const_mat
	test_ref_mat
	test_ref_col
	test_ref_row
	test_ref_mat_ex
	test_mat_copy
	test_mat_compare
	test_dense_eval
	test_mat_subviews
	test_mat_transpose
	test_dense_veval
	test_dense_cast
	test_mat_arith
	test_mat_emath
	test_mat_ecomp
	test_mat_elogical
	test_vec_rep
    test_cpd_ewise
	test_mat_reduce
	test_partial_reduce
)
	
foreach(tname ${LMAT_TEST_NAMES})
	target_link_libraries(${tname} test_main)
endforeach(tname)	


# Add tests

foreach(tname ${LMAT_TEST_NAMES})
	add_test(NAME ${tname} COMMAND ${tname})
endforeach(tname)





