# Configuration of compiler options

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
