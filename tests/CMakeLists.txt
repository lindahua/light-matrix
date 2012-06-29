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


# Header file groups (to be used as dependencies)

set(INC ../light_mat)

set(COMMON_BASIC_HS 
    ${INC}/config/user_config.h 
    ${INC}/config/platform_config.h
    ${INC}/config/config.h
    ${INC}/core/prim_types.h
    ${INC}/core/lang_base.h
    ${INC}/core/metaprog_base.h
    ${INC}/core/type_traits.h
    ${INC}/core/arg_check.h
    ${INC}/core/range.h
    ${INC}/core/basic_defs.h)
    
set(COMMON_ARRAY_HS
    ${INC}/core/mem_op.h
    ${INC}/core/mem_alloc.h
    ${INC}/core/array.h)
    
set(CORE_HS
    ${COMMON_BASIC_HS}
    ${COMMON_ARRAY_HS})
    
set(MATRIX_BASE_HS
    ${INC}/matrix/matrix_fwd.h 
    ${INC}/matrix/matrix_meta.h
    ${INC}/matrix/matrix_indexer.h
    ${INC}/matrix/matrix_concepts.h
    ${INC}/matrix/matrix_properties.h
    ${INC}/matrix/matrix_check.h
    ${INC}/matrix/matrix_generators.h
    ${INC}/matrix/matrix_copy.h
    ${INC}/matrix/bits/matrix_copy_internal.h
    ${INC}/matrix/matrix_print.h
    ${INC}/matrix/matrix_base.h)

# Executables

add_library(test_main STATIC test_main.cpp)

add_executable(test_arrays ${CORE_HS} test_arrays.cpp)

set(DENSE_MAT_HS
    ${MATRIX_BASE_HS}
    ${INC}/matrix/dense_matrix.h
    ${INC}/matrix/bits/dense_matrix_internal.h)

add_executable(test_dense_mat ${DENSE_MAT_HS} test_dense_mat.cpp)
add_executable(test_dense_col ${DENSE_MAT_HS} test_dense_col.cpp)
add_executable(test_dense_row ${DENSE_MAT_HS} test_dense_row.cpp)  

set(REF_MAT_HS
    ${MATRIX_BASE_HS}
    ${INC}/matrix/ref_matrix.h
    ${INC}/matrix/ref_matrix_ex.h
    ${INC}/matrix/bits/ref_matrix_internal.h 
    ${INC}/matrix/bits/ref_matrix_ex_internal.h)
    
add_executable(test_ref_mat ${REF_MAT_HS} test_ref_mat.cpp)
add_executable(test_ref_col ${REF_MAT_HS} test_ref_col.cpp)
add_executable(test_ref_row ${REF_MAT_HS} test_ref_row.cpp)
add_executable(test_ref_mat_ex ${REF_MAT_HS} test_ref_mat_ex.cpp)

add_executable(test_mat_copy ${REF_MAT_HS} test_mat_copy.cpp)

# Link to test_main

target_link_libraries(test_arrays test_main)

target_link_libraries(test_dense_mat test_main)
target_link_libraries(test_dense_col test_main)
target_link_libraries(test_dense_row test_main)

target_link_libraries(test_ref_mat test_main)
target_link_libraries(test_ref_col test_main)
target_link_libraries(test_ref_row test_main)
target_link_libraries(test_ref_mat_ex test_main)

target_link_libraries(test_mat_copy test_main)

# Set output directories

set(ALL_EXECUTABLES 
    test_arrays
    test_dense_mat
    test_dense_col
    test_dense_row
    test_ref_mat
    test_ref_col
    test_ref_row
    test_ref_mat_ex
    test_mat_copy)
    
set_target_properties(${ALL_EXECUTABLES}
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "bin")    


# Add Tests

add_test(NAME arrays COMMAND test_arrays)

add_test(NAME dense_mat COMMAND test_dense_mat)
add_test(NAME dense_col COMMAND test_dense_col)
add_test(NAME dense_row COMMAND test_dense_row)

add_test(NAME ref_mat COMMAND test_ref_mat)
add_test(NAME ref_col COMMAND test_ref_col)
add_test(NAME ref_row COMMAND test_ref_row)
add_test(NAME ref_mat_ex COMMAND test_ref_mat_ex)

add_test(NAME mat_copy COMMAND test_mat_copy)




