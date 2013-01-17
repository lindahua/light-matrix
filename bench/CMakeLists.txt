# The CMakeLists file for Light-Matrix (benchmark)

cmake_minimum_required(VERSION 2.8)
project(LightMatrix)

enable_testing()

include_directories(..)
include_directories($ENV{LIGHT_TEST_HOME})

set(EXECUTABLE_OUTPUT_PATH bin)

# Compiler configuration

include("../cmake_modules/CompilerConfig.cmake")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lrt")
endif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")

message(STATUS "[LMAT] Updated CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")


# third-party libraries

set(CMAKE_MODULE_PATH "../cmake_modules" ${CMAKE_MODULE_PATH})

# Intel library 

find_package(ICCLib)
if (ICCLIB_FOUND)
message(STATUS "[LMAT] ICC Library found: ${ICCLIB_PATH}")
else (ICCLIB_FOUND)
message(STATUS "[LMAT] ICC Library not found")
endif (ICCLIB_FOUND)

set(SVML_FOUND ICCLIB_FOUND)

# Add executables

add_executable(bench_copy ${COMMON_HS} bench_copy.cpp)
add_executable(bench_arith ${COMMON_HS} bench_arith.cpp)

if (SVML_FOUND)
add_executable(bench_math_svml ${COMMON_HS} bench_math.cpp)
endif (SVML_FOUND)

add_executable(bench_reduction ${COMMON_HS} bench_reduction.cpp)
add_executable(bench_prng ${COMMON_HS} bench_prng.cpp)

# Special Linking

set(BENCH_ON_SVML
    bench_math_svml
    bench_prng)

if (SVML_FOUND)
foreach (tname ${BENCH_ON_SVML}) 
    set_target_properties(${tname}
	    PROPERTIES
	    COMPILE_FLAGS "-DLMAT_USE_INTEL_SVML")
    target_link_libraries(${tname} ${SVML_LIBRARY})
endforeach (tname)
endif (SVML_FOUND)
