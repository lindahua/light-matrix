/**
 * @file simd_arch.h
 *
 * @brief SIMD instruction architecture detection
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_SIMD_ARCH_H_
#define LIGHTMAT_SIMD_ARCH_H_

#include <light_mat/config/config.h>

#ifndef LMAT_SIMD_LEVEL
#if defined ( __AVX2__ )
#define LMAT_SIMD_LEVEL 8
#elif defined ( __AVX__ )
#define LMAT_SIMD_LEVEL 7
#elif defined ( __SSE4_2__ )
#define LMAT_SIMD_LEVEL 6
#elif defined ( __SSE4_1__ )
#define LMAT_SIMD_LEVEL 5
#elif defined ( __SSSE3__ )
#define LMAT_SIMD_LEVEL 4
#elif defined ( __SSE3__ )
#define LMAT_SIMD_LEVEL 3
#elif defined ( __SSE2__ ) || defined ( __x86_64__ )
#define LMAT_SIMD_LEVEL 2
#elif defined ( __SSE__ )
#define LMAT_SIMD_LEVEL 1
#elif defined ( _M_IX86_FP )
#define LMAT_SIMD_LEVEL _M_IX86_FP
#else
#define LMAT_SIMD_LEVEL 0
#endif // instruction set defines
#endif

#if LMAT_SIMD_LEVEL >= 1
#define LMAT_HAS_SSE
#endif

#if LMAT_SIMD_LEVEL >= 2
#define LMAT_HAS_SSE2
#endif

#if LMAT_SIMD_LEVEL >= 3
#define LMAT_HAS_SSE3
#endif

#if LMAT_SIMD_LEVEL >= 4
#define LMAT_HAS_SSSE3
#endif

#if LMAT_SIMD_LEVEL >= 5
#define LMAT_HAS_SSE4_1
#endif

#if LMAT_SIMD_LEVEL >= 6
#define LMAT_HAS_SSE4_2
#endif

#if LMAT_SIMD_LEVEL >= 7
#define LMAT_HAS_AVX
#endif

#if LMAT_SIMD_LEVEL >= 8
#define LMAT_HAS_AVX2
#endif


#if (!defined(LMAT_HAS_SSE2))
#error LightMatrix requires at least SSE2 support.
#endif


#endif /* SIMD_ARCH_H_ */
