/**
 * @file avx_pred.h
 *
 * @brief AVX predicates
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_AVX_PRED_H_
#define LIGHTMAT_AVX_PRED_H_

#include <light_mat/simd/avx_packs.h>
#include <light_mat/simd/avx_bpacks.h>
#include "internal/sse_fpclass_impl.h"

namespace lmat { namespace meta {

	// comparison

	LMAT_DEFINE_HAS_AVX_SUPPORT( eq_ )
	LMAT_DEFINE_HAS_AVX_SUPPORT( ne_ )
	LMAT_DEFINE_HAS_AVX_SUPPORT( gt_ )
	LMAT_DEFINE_HAS_AVX_SUPPORT( ge_ )
	LMAT_DEFINE_HAS_AVX_SUPPORT( lt_ )
	LMAT_DEFINE_HAS_AVX_SUPPORT( le_ )

	LMAT_DEFINE_HAS_AVX_SUPPORT( logical_not_ )
	LMAT_DEFINE_HAS_AVX_SUPPORT( logical_and_ )
	LMAT_DEFINE_HAS_AVX_SUPPORT( logical_or_ )
	LMAT_DEFINE_HAS_AVX_SUPPORT( logical_eq_ )
	LMAT_DEFINE_HAS_AVX_SUPPORT( logical_ne_ )

	// numeric predicates

	LMAT_DEFINE_HAS_AVX_SUPPORT( signbit_ )
	LMAT_DEFINE_HAS_AVX_SUPPORT( isfinite_ )
	LMAT_DEFINE_HAS_AVX_SUPPORT( isinf_ )
	LMAT_DEFINE_HAS_AVX_SUPPORT( isnan_ )

} }


namespace lmat
{

	/********************************************
	 *
	 *  comparison operator
	 *
	 ********************************************/

	LMAT_ENSURE_INLINE
	inline avx_f32bpk operator == (const avx_f32pk& a, const avx_f32pk& b)
	{
		return _mm256_cmp_ps(a, b, _CMP_EQ_OQ);
	}

	LMAT_ENSURE_INLINE
	inline avx_f64bpk operator == (const avx_f64pk& a, const avx_f64pk& b)
	{
		return _mm256_cmp_pd(a, b, _CMP_EQ_OQ);
	}

	LMAT_ENSURE_INLINE
	inline avx_f32bpk operator != (const avx_f32pk& a, const avx_f32pk& b)
	{
		return _mm256_cmp_ps(a, b, _CMP_NEQ_OQ);
	}

	LMAT_ENSURE_INLINE
	inline avx_f64bpk operator != (const avx_f64pk& a, const avx_f64pk& b)
	{
		return _mm256_cmp_pd(a, b, _CMP_NEQ_OQ);
	}

	LMAT_ENSURE_INLINE
	inline avx_f32bpk operator > (const avx_f32pk& a, const avx_f32pk& b)
	{
		return _mm256_cmp_ps(a, b, _CMP_GT_OQ);
	}

	LMAT_ENSURE_INLINE
	inline avx_f64bpk operator > (const avx_f64pk& a, const avx_f64pk& b)
	{
		return _mm256_cmp_pd(a, b, _CMP_GT_OQ);
	}

	LMAT_ENSURE_INLINE
	inline avx_f32bpk operator >= (const avx_f32pk& a, const avx_f32pk& b)
	{
		return _mm256_cmp_ps(a, b, _CMP_GE_OQ);
	}

	LMAT_ENSURE_INLINE
	inline avx_f64bpk operator >= (const avx_f64pk& a, const avx_f64pk& b)
	{
		return _mm256_cmp_pd(a, b, _CMP_GE_OQ);
	}

	LMAT_ENSURE_INLINE
	inline avx_f32bpk operator < (const avx_f32pk& a, const avx_f32pk& b)
	{
		return _mm256_cmp_ps(a, b, _CMP_LT_OQ);
	}

	LMAT_ENSURE_INLINE
	inline avx_f64bpk operator < (const avx_f64pk& a, const avx_f64pk& b)
	{
		return _mm256_cmp_pd(a, b, _CMP_LT_OQ);
	}


	LMAT_ENSURE_INLINE
	inline avx_f32bpk operator <= (const avx_f32pk& a, const avx_f32pk& b)
	{
		return _mm256_cmp_ps(a, b, _CMP_LE_OQ);
	}

	LMAT_ENSURE_INLINE
	inline avx_f64bpk operator <= (const avx_f64pk& a, const avx_f64pk& b)
	{
		return _mm256_cmp_pd(a, b, _CMP_LE_OQ);
	}


	/********************************************
	 *
	 *  logical operations
	 *
	 ********************************************/

	LMAT_ENSURE_INLINE
	inline avx_f32bpk operator ~ (const avx_f32bpk& a)
	{
		return _mm256_xor_ps(a, _mm256_castsi256_ps(_mm256_set1_epi32(-1)));
	}

	LMAT_ENSURE_INLINE
	inline avx_f64bpk operator ~ (const avx_f64bpk& a)
	{
		return _mm256_xor_pd(a, _mm256_castsi256_pd(_mm256_set1_epi32(-1)));
	}

	LMAT_ENSURE_INLINE
	inline avx_f32bpk operator & (const avx_f32bpk& a, const avx_f32bpk& b)
	{
		return _mm256_and_ps(a, b);
	}

	LMAT_ENSURE_INLINE
	inline avx_f64bpk operator & (const avx_f64bpk& a, const avx_f64bpk& b)
	{
		return _mm256_and_pd(a, b);
	}

	LMAT_ENSURE_INLINE
	inline avx_f32bpk operator | (const avx_f32bpk& a, const avx_f32bpk& b)
	{
		return _mm256_or_ps(a, b);
	}

	LMAT_ENSURE_INLINE
	inline avx_f64bpk operator | (const avx_f64bpk& a, const avx_f64bpk& b)
	{
		return _mm256_or_pd(a, b);
	}

	LMAT_ENSURE_INLINE
	inline avx_f32bpk operator != (const avx_f32bpk& a, const avx_f32bpk& b)
	{
		return _mm256_xor_ps(a, b);
	}

	LMAT_ENSURE_INLINE
	inline avx_f64bpk operator != (const avx_f64bpk& a, const avx_f64bpk& b)
	{
		return _mm256_xor_pd(a, b);
	}

	LMAT_ENSURE_INLINE
	inline avx_f32bpk operator == (const avx_f32bpk& a, const avx_f32bpk& b)
	{
		return ~(a != b);
	}

	LMAT_ENSURE_INLINE
	inline avx_f64bpk operator == (const avx_f64bpk& a, const avx_f64bpk& b)
	{
		return ~(a != b);
	}

	LMAT_ENSURE_INLINE
	inline avx_f32bpk& operator &= (avx_f32bpk& a, const avx_f32bpk& b)
	{
		a = a & b;
		return a;
	}

	LMAT_ENSURE_INLINE
	inline avx_f64bpk& operator &= (avx_f64bpk& a, const avx_f64bpk& b)
	{
		a = a & b;
		return a;
	}

	LMAT_ENSURE_INLINE
	inline avx_f32bpk& operator |= (avx_f32bpk& a, const avx_f32bpk& b)
	{
		a = a | b;
		return a;
	}

	LMAT_ENSURE_INLINE
	inline avx_f64bpk& operator |= (avx_f64bpk& a, const avx_f64bpk& b)
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
	inline avx_f32bpk signbit(const avx_f32pk& a)
	{
		return internal::combine_m128(
				internal::sse_is_neg_ps(a.get_low()),
				internal::sse_is_neg_ps(a.get_high()));
	}

	LMAT_ENSURE_INLINE
	inline avx_f64bpk signbit(const avx_f64pk& a)
	{
		return internal::combine_m128d(
				internal::sse_is_neg_pd(a.get_low()),
				internal::sse_is_neg_pd(a.get_high()));
	}


	LMAT_ENSURE_INLINE
	inline avx_f32bpk isfinite(const avx_f32pk& a)
	{
		return internal::combine_m128(
				internal::sse_is_finite_ps(a.get_low()),
				internal::sse_is_finite_ps(a.get_high()));
	}

	LMAT_ENSURE_INLINE
	inline avx_f64bpk isfinite(const avx_f64pk& a)
	{
		return internal::combine_m128d(
				internal::sse_is_finite_pd(a.get_low()),
				internal::sse_is_finite_pd(a.get_high()));
	}


	LMAT_ENSURE_INLINE
	inline avx_f32bpk isinf(const avx_f32pk& a)
	{
		return internal::combine_m128(
				internal::sse_is_inf_ps(a.get_low()),
				internal::sse_is_inf_ps(a.get_high()));
	}

	LMAT_ENSURE_INLINE
	inline avx_f64bpk isinf(const avx_f64pk& a)
	{
		return internal::combine_m128d(
				internal::sse_is_inf_pd(a.get_low()),
				internal::sse_is_inf_pd(a.get_high()));
	}


	LMAT_ENSURE_INLINE
	inline avx_f32bpk isnan(const avx_f32pk& a)
	{
		return internal::combine_m128(
				internal::sse_is_nan_ps(a.get_low()),
				internal::sse_is_nan_ps(a.get_high()));
	}

	LMAT_ENSURE_INLINE
	inline avx_f64bpk isnan(const avx_f64pk& a)
	{
		return internal::combine_m128d(
				internal::sse_is_nan_pd(a.get_low()),
				internal::sse_is_nan_pd(a.get_high()));
	}

} }

#endif
