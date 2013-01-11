/**
 * @file simd_special.h
 *
 * @brief SIMD special functions
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_SIMD_SPECIAL_H_
#define LIGHTMAT_SIMD_SPECIAL_H_

#include <light_mat/math/sse_math.h>
#ifdef LMAT_HAS_AVX
#include <light_mat/math/avx_math.h>
#endif
#include <light_mat/math/special.h>


/************************************************
 *
 *  SIMD support declaration
 *
 ***********************************************/

namespace lmat { namespace meta {

#ifdef LMAT_HAS_EXTERN_SSE_ERF
	LMAT_DEFINE_HAS_SSE_SUPPORT( erf_ )
	LMAT_DEFINE_HAS_SSE_SUPPORT( erfc_ )
#endif

#ifdef LMAT_HAS_EXTERN_AVX_ERF
	LMAT_DEFINE_HAS_AVX_SUPPORT( erf_ )
	LMAT_DEFINE_HAS_AVX_SUPPORT( erfc_ )
#endif

#ifdef LMAT_HAS_EXTERN_SSE_GAMMA
	LMAT_DEFINE_HAS_SSE_SUPPORT( lgamma_ )
	LMAT_DEFINE_HAS_SSE_SUPPORT( tgamma_ )
#endif

#ifdef LMAT_HAS_EXTERN_AVX_GAMMA
	LMAT_DEFINE_HAS_AVX_SUPPORT( lgamma_ )
	LMAT_DEFINE_HAS_AVX_SUPPORT( tgamma_ )
#endif

} }

namespace lmat { namespace math {

	/************************************************
	 *
	 *  emulation by scalar-function
	 *
	 ***********************************************/

	namespace internal
	{
		LMAT_DEFINE_SSE_MATH_EMULATE_1( erf )
		LMAT_DEFINE_SSE_MATH_EMULATE_1( erfc )

		LMAT_DEFINE_SSE_MATH_EMULATE_1( lgamma )
		LMAT_DEFINE_SSE_MATH_EMULATE_1( tgamma )

		LMAT_DEFINE_AVX_MATH_EMULATE_1( erf )
		LMAT_DEFINE_AVX_MATH_EMULATE_1( erfc )

		LMAT_DEFINE_AVX_MATH_EMULATE_1( lgamma )
		LMAT_DEFINE_AVX_MATH_EMULATE_1( tgamma )
	}


	/************************************************
	 *
	 *  functions for end-users
	 *
	 ***********************************************/

#ifdef LMAT_HAS_EXTERN_SSE_ERF
	LMAT_ACTIVATE_SSE_EXTERN_1( erf )
	LMAT_ACTIVATE_SSE_EXTERN_1( erfc )
#elif (defined(LMAT_ENABLE_SIMD_EMULATE))
	LMAT_ACTIVATE_SSE_MATH_EMULATE_1( erf )
	LMAT_ACTIVATE_SSE_MATH_EMULATE_1( erfc )
#endif

#ifdef LMAT_HAS_EXTERN_SSE_GAMMA
	LMAT_ACTIVATE_SSE_EXTERN_1( lgamma )
	LMAT_ACTIVATE_SSE_EXTERN_1( tgamma )
#elif (defined(LMAT_ENABLE_SIMD_EMULATE))
	LMAT_ACTIVATE_SSE_MATH_EMULATE_1( lgamma )
	LMAT_ACTIVATE_SSE_MATH_EMULATE_1( tgamma )
#endif

#ifdef LMAT_HAS_EXTERN_AVX_ERF
	LMAT_ACTIVATE_AVX_EXTERN_1( erf )
	LMAT_ACTIVATE_AVX_EXTERN_1( erfc )
#elif (defined(LMAT_ENABLE_SIMD_EMULATE))
	LMAT_ACTIVATE_AVX_MATH_EMULATE_1( erf )
	LMAT_ACTIVATE_AVX_MATH_EMULATE_1( erfc )
#endif

#ifdef LMAT_HAS_EXTERN_AVX_GAMMA
	LMAT_ACTIVATE_AVX_EXTERN_1( lgamma )
	LMAT_ACTIVATE_AVX_EXTERN_1( tgamma )
#elif (defined(LMAT_ENABLE_SIMD_EMULATE))
	LMAT_ACTIVATE_AVX_MATH_EMULATE_1( lgamma )
	LMAT_ACTIVATE_AVX_MATH_EMULATE_1( tgamma )
#endif


} }


#endif
