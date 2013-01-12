/**
 * @file sse_pred.h
 *
 * @brief SSE predicate functions
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_SSE_PRED_H_
#define LIGHTMAT_SSE_PRED_H_

#include <light_mat/simd/sse_packs.h>
#include <light_mat/simd/sse_bpacks.h>
#include "internal/sse_fpclass_impl.h"

namespace lmat { namespace meta {

	LMAT_DEFINE_HAS_SSE_SUPPORT( eq_ )
	LMAT_DEFINE_HAS_SSE_SUPPORT( ne_ )
	LMAT_DEFINE_HAS_SSE_SUPPORT( gt_ )
	LMAT_DEFINE_HAS_SSE_SUPPORT( ge_ )
	LMAT_DEFINE_HAS_SSE_SUPPORT( lt_ )
	LMAT_DEFINE_HAS_SSE_SUPPORT( le_ )

	LMAT_DEFINE_HAS_SSE_SUPPORT( logical_not_ )
	LMAT_DEFINE_HAS_SSE_SUPPORT( logical_and_ )
	LMAT_DEFINE_HAS_SSE_SUPPORT( logical_or_ )
	LMAT_DEFINE_HAS_SSE_SUPPORT( logical_eq_ )
	LMAT_DEFINE_HAS_SSE_SUPPORT( logical_ne_ )

	LMAT_DEFINE_HAS_SSE_SUPPORT( signbit_ )
	LMAT_DEFINE_HAS_SSE_SUPPORT( isfinite_ )
	LMAT_DEFINE_HAS_SSE_SUPPORT( isinf_ )
	LMAT_DEFINE_HAS_SSE_SUPPORT( isnan_ )

} }


namespace lmat
{

	/********************************************
	 *
	 *  comparison operator
	 *
	 ********************************************/

	LMAT_ENSURE_INLINE
	inline sse_f32bpk operator == (const sse_f32pk& a, const sse_f32pk& b)
	{
		return _mm_cmpeq_ps(a, b);
	}

	LMAT_ENSURE_INLINE
	inline sse_f64bpk operator == (const sse_f64pk& a, const sse_f64pk& b)
	{
		return _mm_cmpeq_pd(a, b);
	}

	LMAT_ENSURE_INLINE
	inline sse_f32bpk operator != (const sse_f32pk& a, const sse_f32pk& b)
	{
		return _mm_cmpneq_ps(a, b);
	}

	LMAT_ENSURE_INLINE
	inline sse_f64bpk operator != (const sse_f64pk& a, const sse_f64pk& b)
	{
		return _mm_cmpneq_pd(a, b);
	}

	LMAT_ENSURE_INLINE
	inline sse_f32bpk operator > (const sse_f32pk& a, const sse_f32pk& b)
	{
		return _mm_cmpgt_ps(a, b);
	}

	LMAT_ENSURE_INLINE
	inline sse_f64bpk operator > (const sse_f64pk& a, const sse_f64pk& b)
	{
		return _mm_cmpgt_pd(a, b);
	}

	LMAT_ENSURE_INLINE
	inline sse_f32bpk operator >= (const sse_f32pk& a, const sse_f32pk& b)
	{
		return _mm_cmpge_ps(a, b);
	}

	LMAT_ENSURE_INLINE
	inline sse_f64bpk operator >= (const sse_f64pk& a, const sse_f64pk& b)
	{
		return _mm_cmpge_pd(a, b);
	}

	LMAT_ENSURE_INLINE
	inline sse_f32bpk operator < (const sse_f32pk& a, const sse_f32pk& b)
	{
		return _mm_cmplt_ps(a, b);
	}

	LMAT_ENSURE_INLINE
	inline sse_f64bpk operator < (const sse_f64pk& a, const sse_f64pk& b)
	{
		return _mm_cmplt_pd(a, b);
	}


	LMAT_ENSURE_INLINE
	inline sse_f32bpk operator <= (const sse_f32pk& a, const sse_f32pk& b)
	{
		return _mm_cmple_ps(a, b);
	}

	LMAT_ENSURE_INLINE
	inline sse_f64bpk operator <= (const sse_f64pk& a, const sse_f64pk& b)
	{
		return _mm_cmple_pd(a, b);
	}


	/********************************************
	 *
	 *  logical operations
	 *
	 ********************************************/

	LMAT_ENSURE_INLINE
	inline sse_f32bpk operator ~ (const sse_f32bpk& a)
	{
		return _mm_castsi128_ps(
				_mm_cmpeq_epi32(_mm_castps_si128(a), _mm_setzero_si128()));
	}

	LMAT_ENSURE_INLINE
	inline sse_f64bpk operator ~ (const sse_f64bpk& a)
	{
		return _mm_castsi128_pd(
				_mm_cmpeq_epi64(_mm_castpd_si128(a), _mm_setzero_si128()));
	}

	LMAT_ENSURE_INLINE
	inline sse_f32bpk operator & (const sse_f32bpk& a, const sse_f32bpk& b)
	{
		return _mm_and_ps(a, b);
	}

	LMAT_ENSURE_INLINE
	inline sse_f64bpk operator & (const sse_f64bpk& a, const sse_f64bpk& b)
	{
		return _mm_and_pd(a, b);
	}

	LMAT_ENSURE_INLINE
	inline sse_f32bpk operator | (const sse_f32bpk& a, const sse_f32bpk& b)
	{
		return _mm_or_ps(a, b);
	}

	LMAT_ENSURE_INLINE
	inline sse_f64bpk operator | (const sse_f64bpk& a, const sse_f64bpk& b)
	{
		return _mm_or_pd(a, b);
	}

	LMAT_ENSURE_INLINE
	inline sse_f32bpk operator == (const sse_f32bpk& a, const sse_f32bpk& b)
	{
		return _mm_castsi128_ps(
				_mm_cmpeq_epi32(_mm_castps_si128(a), _mm_castps_si128(b)));
	}

	LMAT_ENSURE_INLINE
	inline sse_f64bpk operator == (const sse_f64bpk& a, const sse_f64bpk& b)
	{
		return _mm_castsi128_pd(
				_mm_cmpeq_epi64(_mm_castpd_si128(a), _mm_castpd_si128(b)));
	}

	LMAT_ENSURE_INLINE
	inline sse_f32bpk operator != (const sse_f32bpk& a, const sse_f32bpk& b)
	{
		return _mm_xor_ps(a, b);
	}

	LMAT_ENSURE_INLINE
	inline sse_f64bpk operator != (const sse_f64bpk& a, const sse_f64bpk& b)
	{
		return _mm_xor_pd(a, b);
	}


	LMAT_ENSURE_INLINE
	inline sse_f32bpk& operator &= (sse_f32bpk& a, const sse_f32bpk& b)
	{
		a = a & b;
		return a;
	}

	LMAT_ENSURE_INLINE
	inline sse_f64bpk& operator &= (sse_f64bpk& a, const sse_f64bpk& b)
	{
		a = a & b;
		return a;
	}

	LMAT_ENSURE_INLINE
	inline sse_f32bpk& operator |= (sse_f32bpk& a, const sse_f32bpk& b)
	{
		a = a | b;
		return a;
	}

	LMAT_ENSURE_INLINE
	inline sse_f64bpk& operator |= (sse_f64bpk& a, const sse_f64bpk& b)
	{
		a = a | b;
		return a;
	}
}


namespace lmat { namespace math {

	/********************************************
	 *
	 *  FP classification
	 *
	 ********************************************/

	LMAT_ENSURE_INLINE
	inline sse_f32bpk signbit(const sse_f32pk& a)
	{
		return internal::sse_is_neg_ps(a);
	}

	LMAT_ENSURE_INLINE
	inline sse_f64bpk signbit(const sse_f64pk& a)
	{
		return internal::sse_is_neg_pd(a);
	}


	LMAT_ENSURE_INLINE
	inline sse_f32bpk isfinite(const sse_f32pk& a)
	{
		return internal::sse_is_finite_ps(a);
	}

	LMAT_ENSURE_INLINE
	inline sse_f64bpk isfinite(const sse_f64pk& a)
	{
		return internal::sse_is_finite_pd(a);
	}


	LMAT_ENSURE_INLINE
	inline sse_f32bpk isinf(const sse_f32pk& a)
	{
		return internal::sse_is_inf_ps(a);
	}

	LMAT_ENSURE_INLINE
	inline sse_f64bpk isinf(const sse_f64pk& a)
	{
		return internal::sse_is_inf_pd(a);
	}


	LMAT_ENSURE_INLINE
	inline sse_f32bpk isnan(const sse_f32pk& a)
	{
		return internal::sse_is_nan_ps(a);
	}

	LMAT_ENSURE_INLINE
	inline sse_f64bpk isnan(const sse_f64pk& a)
	{
		return internal::sse_is_nan_pd(a);
	}

} }

#endif
