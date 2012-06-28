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
    ${INC}/core/basic_defs.h)
    
set(COMMON_ARRAY_HS
    ${INC}/core/mem_op.h
    ${INC}/core/mem_alloc.h
    ${INC}/core/array.h)
    
set(CORE_HS
    ${COMMON_BASIC_HS}
    ${COMMON_ARRAY_HS})

# Executables

add_library(test_main STATIC test_main.cpp)

add_executable(test_arrays ${CORE_HS} test_arrays.cpp)

target_link_libraries(test_arrays test_main)

set(ALL_EXECUTABLES 
    test_arrays)
    
set_target_properties(${ALL_EXECUTABLES}
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "bin")    


# Add Tests

add_test(NAME arrays COMMAND test_arrays)




