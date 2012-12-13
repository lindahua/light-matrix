/**
 * @file simd_arch.h
 *
 * @brief SIMD instruction architecture detection
 *
 * @author Dahua Lin
 */

#ifndef LIGHTMAT_SIMD_ARCH_H_
#define LIGHTMAT_SIMD_ARCH_H_

#include <light_mat/config/config.h>

#ifndef LMAT_SIMD
#if defined ( __AVX2__ )
#define LMAT_SIMD 8
#elif defined ( __AVX__ )
#define LMAT_SIMD 7
#elif defined ( __SSE4_2__ )
#define LMAT_SIMD 6
#elif defined ( __SSE4_1__ )
#define LMAT_SIMD 5
#elif defined ( __SSSE3__ )
#define LMAT_SIMD 4
#elif defined ( __SSE3__ )
#define LMAT_SIMD 3
#elif defined ( __SSE2__ ) || defined ( __x86_64__ )
#define LMAT_SIMD 2
#elif defined ( __SSE__ )
#define LMAT_SIMD 1
#elif defined ( _M_IX86_FP )
#define LMAT_SIMD _M_IX86_FP
#else
#define LMAT_SIMD 0
#endif // instruction set defines
#endif


#endif /* SIMD_ARCH_H_ */
