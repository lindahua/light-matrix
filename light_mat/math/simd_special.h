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

#ifdef LMAT_HAS_EXTERN_SSE_ERFINV
	LMAT_DEFINE_HAS_SSE_SUPPORT( erfinv_ )
	LMAT_DEFINE_HAS_SSE_SUPPORT( erfcinv_ )
#endif

#ifdef LMAT_HAS_EXTERN_AVX_ERFINV
	LMAT_DEFINE_HAS_AVX_SUPPORT( erfinv_ )
	LMAT_DEFINE_HAS_AVX_SUPPORT( erfcinv_ )
#endif

#ifdef LMAT_HAS_EXTERN_SSE_NORMINV
	LMAT_DEFINE_HAS_SSE_SUPPORT( norminv_ )
#endif

#ifdef LMAT_HAS_EXTERN_AVX_NORMINV
	LMAT_DEFINE_HAS_AVX_SUPPORT( norminv_ )
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
		LMAT_DEFINE_AVX_MATH_EMULATE_1( erf )
		LMAT_DEFINE_SSE_MATH_EMULATE_1( erfc )
		LMAT_DEFINE_AVX_MATH_EMULATE_1( erfc )

		LMAT_DEFINE_SSE_MATH_EMULATE_1( lgamma )
		LMAT_DEFINE_AVX_MATH_EMULATE_1( lgamma )
		LMAT_DEFINE_SSE_MATH_EMULATE_1( tgamma )
		LMAT_DEFINE_AVX_MATH_EMULATE_1( tgamma )

		LMAT_DEFINE_SSE_MATH_EMULATE_1( norminv )
		LMAT_DEFINE_AVX_MATH_EMULATE_1( norminv )
	}


	/************************************************
	 *
	 *  functions for end-users
	 *
	 ***********************************************/

	// erf & erfc

#ifdef LMAT_HAS_EXTERN_SSE_ERF
	LMAT_ACTIVATE_SSE_EXTERN_1( erf )
	LMAT_ACTIVATE_SSE_EXTERN_1( erfc )
#elif (defined(LMAT_ENABLE_SIMD_EMULATE))
	LMAT_ACTIVATE_SSE_MATH_EMULATE_1( erf )
	LMAT_ACTIVATE_SSE_MATH_EMULATE_1( erfc )
#endif

#ifdef LMAT_HAS_EXTERN_AVX_ERF
	LMAT_ACTIVATE_AVX_EXTERN_1( erf )
	LMAT_ACTIVATE_AVX_EXTERN_1( erfc )
#elif (defined(LMAT_ENABLE_SIMD_EMULATE))
	LMAT_ACTIVATE_AVX_MATH_EMULATE_1( erf )
	LMAT_ACTIVATE_AVX_MATH_EMULATE_1( erfc )
#endif

	// lgamma & tgamma

#ifdef LMAT_HAS_EXTERN_SSE_GAMMA
	LMAT_ACTIVATE_SSE_EXTERN_1( lgamma )
	LMAT_ACTIVATE_SSE_EXTERN_1( tgamma )
#elif (defined(LMAT_ENABLE_SIMD_EMULATE))
	LMAT_ACTIVATE_SSE_MATH_EMULATE_1( lgamma )
	LMAT_ACTIVATE_SSE_MATH_EMULATE_1( tgamma )
#endif

#ifdef LMAT_HAS_EXTERN_AVX_GAMMA
	LMAT_ACTIVATE_AVX_EXTERN_1( lgamma )
	LMAT_ACTIVATE_AVX_EXTERN_1( tgamma )
#elif (defined(LMAT_ENABLE_SIMD_EMULATE))
	LMAT_ACTIVATE_AVX_MATH_EMULATE_1( lgamma )
	LMAT_ACTIVATE_AVX_MATH_EMULATE_1( tgamma )
#endif

	// norminv

#ifdef LMAT_HAS_EXTERN_SSE_NORMINV
	LMAT_ENSURE_INLINE
	inline sse_f32pk norminv(const sse_f32pk& a)
	{
		return LMAT_SSE_F(cdfnorminv)(a);
	}

	LMAT_ENSURE_INLINE
	inline sse_f64pk norminv(const sse_f64pk& a)
	{
		return LMAT_SSE_D(cdfnorminv)(a);
	}
#endif

#ifdef LMAT_HAS_EXTERN_AVX_NORMINV
	LMAT_ENSURE_INLINE
	inline avx_f32pk norminv(const avx_f32pk& a)
	{
		return LMAT_AVX_F(cdfnorminv)(a);
	}

	LMAT_ENSURE_INLINE
	inline avx_f64pk norminv(const avx_f64pk& a)
	{
		return LMAT_AVX_D(cdfnorminv)(a);
	}
#endif

} }


#endif
