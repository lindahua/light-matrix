# The CMakeLists file for Light-Matrix (benchmark)

cmake_minimum_required(VERSION 2.8)
project(LightMatrix)

enable_testing()

include_directories(..)
include_directories($ENV{LIGHT_TEST_HOME})


# Compiler configuration

set(CMAKE_BUILD_TYPE "Release")

if (MSVC)
	set(LANG_FLAGS "/arch:AVX /EHsc /DLMAT_DIAGNOSIS_LEVEL=0")
	set(WARNING_FLAGS "/W4")
else (MSVC)
	set(LANG_FLAGS "-std=c++0x -pedantic -mavx -DLMAT_DIAGNOSIS_LEVEL=0")
	set(WARNING_FLAGS "-Wall -Wextra -Wconversion -Wformat -Wno-unused-parameter ")
endif (MSVC)

if (${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
	set(LANG_FLAGS "${LANG_FLAGS} -stdlib=libc++ -Qunused-arguments")
	set(CMAKE_CXX_COMPILER "clang++")
endif (${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")

set(CMAKE_CXX_FLAGS "${LANG_FLAGS} ${WARNING_FLAGS}")

set(EXECUTABLE_OUTPUT_PATH bin)

# Third-party libraries

# real-time library

if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lrt")
endif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")

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

# Add executables

add_executable(bench_copy ${COMMON_HS} bench_copy.cpp)






