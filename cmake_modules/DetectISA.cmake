# Detection of Instruction Set Architectures

set(TARGET_ISA "$ENV{LMAT_TARGET_ISA}")

if (TARGET_ISA)
    message(STATUS "[LMAT] TARGET_ISA = ${TARGET_ISA}")
else (TARGET_ISA)
    set(TARGET_ISA "avx")
    message(STATUS "[LMAT] Environment variable TARGET_ISA is not set, using default:")
    message(STATUS "[LMAT] TARGET_ISA = ${TARGET_ISA}")
endif (TARGET_ISA)


# set ALLOW_* variables

if (${TARGET_ISA} STREQUAL "sse2")
    set(ALLOW_SSE2   "yes")
    set(ALLOW_SSE3   "no")
    set(ALLOW_SSSE3  "no")
    set(ALLOW_SSE4_1 "no")
    set(ALLOW_SSE4_2 "no")
    set(ALLOW_AVX    "no")
endif (${TARGET_ISA} STREQUAL "sse2")

if (${TARGET_ISA} STREQUAL "sse3")
    set(ALLOW_SSE2   "yes")
    set(ALLOW_SSE3   "yes")
    set(ALLOW_SSSE3  "no")
    set(ALLOW_SSE4_1 "no")
    set(ALLOW_SSE4_2 "no")
    set(ALLOW_AVX    "no")
endif (${TARGET_ISA} STREQUAL "sse3")

if (${TARGET_ISA} STREQUAL "ssse3")
    set(ALLOW_SSE2   "yes")
    set(ALLOW_SSE3   "yes")
    set(ALLOW_SSSE3  "yes")
    set(ALLOW_SSE4_1 "no")
    set(ALLOW_SSE4_2 "no")
    set(ALLOW_AVX    "no")
endif (${TARGET_ISA} STREQUAL "ssse3")

if (${TARGET_ISA} STREQUAL "sse4.1")
    set(ALLOW_SSE2   "yes")
    set(ALLOW_SSE3   "yes")
    set(ALLOW_SSSE3  "yes")
    set(ALLOW_SSE4_1 "yes")
    set(ALLOW_SSE4_2 "no")
    set(ALLOW_AVX    "no")
endif (${TARGET_ISA} STREQUAL "sse4.1")

if (${TARGET_ISA} STREQUAL "sse4.2")
    set(ALLOW_SSE2   "yes")
    set(ALLOW_SSE3   "yes")
    set(ALLOW_SSSE3  "yes")
    set(ALLOW_SSE4_1 "yes")
    set(ALLOW_SSE4_2 "yes")
    set(ALLOW_AVX    "no")
endif (${TARGET_ISA} STREQUAL "sse4.2")

if (${TARGET_ISA} STREQUAL "avx")
    set(ALLOW_SSE2   "yes")
    set(ALLOW_SSE3   "yes")
    set(ALLOW_SSSE3  "yes")
    set(ALLOW_SSE4_1 "yes")
    set(ALLOW_SSE4_2 "yes")
    set(ALLOW_AVX    "yes")
endif (${TARGET_ISA} STREQUAL "avx")


# set compiler arch flags

if (MSVC)
    if (ALLOW_AVX)
        set(ARCH_FLAG "/arch:AVX")
    else(ALLOW_AVX)
        set(ARCH_FLAG "/arch:SSE2")
    endif (ALLOW_AVX)
else (MSVC)
    set(ARCH_FLAG "-m${TARGET_ISA}")
endif (MSVC)

message(STATUS "[LMAT] ARCH_FLAG = ${ARCH_FLAG}")


