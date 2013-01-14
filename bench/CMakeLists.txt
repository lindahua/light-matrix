# The CMakeLists file for Light-Matrix (benchmark)

cmake_minimum_required(VERSION 2.8)
project(LightMatrix)

enable_testing()

include_directories(..)
include_directories($ENV{LIGHT_TEST_HOME})

set(EXECUTABLE_OUTPUT_PATH bin)

# Compiler configuration

include("../cmake_modules/CompilerConfig.cmake")

if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lrt")
endif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")


# Add executables

add_executable(bench_copy ${COMMON_HS} bench_copy.cpp)






