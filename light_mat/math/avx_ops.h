/**
 * @file avx_ops.h
 *
 * Basic operations on AVX packs
 * 
 * @author Dahua Lin 
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_AVX_OPS_H_
#define LIGHTMAT_AVX_OPS_H_

#include <light_mat/math/avx_packs.h>
#include <light_mat/math/avx_bpacks.h>
#include "internal/numrepr_format.h"


namespace lmat { namespace math {

	/********************************************
	 *
	 *  Floating-point arithmetics
	 *
	 ********************************************/

	LMAT_ENSURE_INLINE
	inline avx_f32pk operator + (const avx_f32pk& a, const avx_f32pk& b)
	{
		return _mm256_add_ps(a.v, b.v);
	}

	LMAT_ENSURE_INLINE
	inline avx_f64pk operator + (const avx_f64pk& a, const avx_f64pk& b)
	{
		return _mm256_add_pd(a.v, b.v);
	}

	LMAT_ENSURE_INLINE
	inline avx_f32pk operator - (const avx_f32pk& a, const avx_f32pk& b)
	{
		return _mm256_sub_ps(a.v, b.v);
	}

	LMAT_ENSURE_INLINE
	inline avx_f64pk operator - (const avx_f64pk& a, const avx_f64pk& b)
	{
		return _mm256_sub_pd(a.v, b.v);
	}

	LMAT_ENSURE_INLINE
	inline avx_f32pk operator * (const avx_f32pk& a, const avx_f32pk& b)
	{
		return _mm256_mul_ps(a.v, b.v);
	}

	LMAT_ENSURE_INLINE
	inline avx_f64pk operator * (const avx_f64pk& a, const avx_f64pk& b)
	{
		return _mm256_mul_pd(a.v, b.v);
	}

	LMAT_ENSURE_INLINE
	inline avx_f32pk operator / (const avx_f32pk& a, const avx_f32pk& b)
	{
		return _mm256_div_ps(a.v, b.v);
	}

	LMAT_ENSURE_INLINE
	inline avx_f64pk operator / (const avx_f64pk& a, const avx_f64pk& b)
	{
		return _mm256_div_pd(a.v, b.v);
	}

	LMAT_ENSURE_INLINE
	inline avx_f32pk operator - (const avx_f32pk& a)
	{
		typedef internal::num_fmt<float> fmt;
		return _mm256_xor_ps(_mm256_castsi256_ps(_mm256_set1_epi32(fmt::sign_bit)), a.v);
	}

	LMAT_ENSURE_INLINE
	inline avx_f64pk operator - (const avx_f64pk& a)
	{
		typedef internal::num_fmt<double> fmt;
		return _mm256_xor_pd(_mm256_castsi256_pd(_mm256_set1_epi64x(fmt::sign_bit)), a.v);
	}

	LMAT_ENSURE_INLINE
	inline avx_f32pk abs(const avx_f32pk& a)
	{
		typedef internal::num_fmt<float> fmt;
		return _mm256_andnot_ps(_mm256_castsi256_ps(_mm256_set1_epi32(fmt::sign_bit)), a.v);
	}

	LMAT_ENSURE_INLINE
	inline avx_f64pk abs(const avx_f64pk& a)
	{
		typedef internal::num_fmt<double> fmt;
		return _mm256_andnot_pd(_mm256_castsi256_pd(_mm256_set1_epi64x(fmt::sign_bit)), a.v);
	}

	LMAT_ENSURE_INLINE
	inline avx_f32pk& operator += (avx_f32pk& a, const avx_f32pk& b)
	{
		a.v = _mm256_add_ps(a.v, b.v);
		return a;
	}

	LMAT_ENSURE_INLINE
	inline avx_f64pk& operator += (avx_f64pk& a, const avx_f64pk& b)
	{
		a.v = _mm256_add_pd(a.v, b.v);
		return a;
	}

	LMAT_ENSURE_INLINE
	inline avx_f32pk& operator -= (avx_f32pk& a, const avx_f32pk& b)
	{
		a.v = _mm256_sub_ps(a.v, b.v);
		return a;
	}

	LMAT_ENSURE_INLINE
	inline avx_f64pk& operator -= (avx_f64pk& a, const avx_f64pk& b)
	{
		a.v = _mm256_sub_pd(a.v, b.v);
		return a;
	}

	LMAT_ENSURE_INLINE
	inline avx_f32pk& operator *= (avx_f32pk& a, const avx_f32pk& b)
	{
		a.v = _mm256_mul_ps(a.v, b.v);
		return a;
	}

	LMAT_ENSURE_INLINE
	inline avx_f64pk& operator *= (avx_f64pk& a, const avx_f64pk& b)
	{
		a.v = _mm256_mul_pd(a.v, b.v);
		return a;
	}

	LMAT_ENSURE_INLINE
	inline avx_f32pk& operator /= (avx_f32pk& a, const avx_f32pk& b)
	{
		a.v = _mm256_div_ps(a.v, b.v);
		return a;
	}

	LMAT_ENSURE_INLINE
	inline avx_f64pk& operator /= (avx_f64pk& a, const avx_f64pk& b)
	{
		a.v = _mm256_div_pd(a.v, b.v);
		return a;
	}


	/********************************************
	 *
	 *  Floating-point min and max
	 *
	 ********************************************/

	LMAT_ENSURE_INLINE
	inline avx_f32pk (min)(const avx_f32pk& a, const avx_f32pk& b)
	{
		return _mm256_min_ps(a.v, b.v);
	}

	LMAT_ENSURE_INLINE
	inline avx_f64pk (min)(const avx_f64pk& a, const avx_f64pk& b)
	{
		return _mm256_min_pd(a.v, b.v);
	}

	LMAT_ENSURE_INLINE
	inline avx_f32pk (max)(const avx_f32pk& a, const avx_f32pk& b)
	{
		return _mm256_max_ps(a.v, b.v);
	}
	LMAT_ENSURE_INLINE
	inline avx_f64pk (max)(const avx_f64pk& a, const avx_f64pk& b)
	{
		return _mm256_max_pd(a.v, b.v);
	}

	/********************************************
	 *
	 *  Simple power functions
	 *
	 ********************************************/

	LMAT_ENSURE_INLINE
	inline avx_f32pk sqr(const avx_f32pk& a)
	{
		return _mm256_mul_ps(a.v, a.v);
	}

	LMAT_ENSURE_INLINE
	inline avx_f64pk sqr(const avx_f64pk& a)
	{
		return _mm256_mul_pd(a.v, a.v);
	}

	LMAT_ENSURE_INLINE
	inline avx_f32pk cube(const avx_f32pk& a)
	{
		return _mm256_mul_ps(_mm256_mul_ps(a.v, a.v), a.v);
	}

	LMAT_ENSURE_INLINE
	inline avx_f64pk cube(const avx_f64pk& a)
	{
		return _mm256_mul_pd(_mm256_mul_pd(a.v, a.v), a.v);
	}

	LMAT_ENSURE_INLINE
	inline avx_f32pk sqrt(const avx_f32pk& a)
	{
		return _mm256_sqrt_ps(a.v);
	}

	LMAT_ENSURE_INLINE
	inline avx_f64pk sqrt(const avx_f64pk& a)
	{
		return _mm256_sqrt_pd(a.v);
	}

	LMAT_ENSURE_INLINE
	inline avx_f32pk rcp(const avx_f32pk& a)
	{
		return _mm256_div_ps(_mm256_set1_ps(1.0f), a.v);
	}

	LMAT_ENSURE_INLINE
	inline avx_f32pk approx_rcp(const avx_f32pk& a)
	{
		return _mm256_rcp_ps(a.v);
	}

	LMAT_ENSURE_INLINE
	inline avx_f64pk rcp(const avx_f64pk& a)
	{
		return _mm256_div_pd(_mm256_set1_pd(1.0), a.v);
	}

	LMAT_ENSURE_INLINE
	inline avx_f32pk rsqrt(const avx_f32pk& a)
	{
		return _mm256_div_ps(_mm256_set1_ps(1.0f), _mm256_sqrt_ps(a.v));
	}

	LMAT_ENSURE_INLINE
	inline avx_f32pk approx_rsqrt(const avx_f32pk& a)
	{
		return _mm256_rsqrt_ps(a.v);
	}

	LMAT_ENSURE_INLINE
	inline avx_f64pk rsqrt(const avx_f64pk& a)
	{
		return _mm256_div_pd(_mm256_set1_pd(1.0), _mm256_sqrt_pd(a.v));
	}


	/********************************************
	 *
	 *  comparison operator
	 *
	 ********************************************/

	/********************************************
	 *
	 *  comparison operator
	 *
	 ********************************************/

	LMAT_ENSURE_INLINE
	inline avx_f32bpk operator == (const avx_f32pk& a, const avx_f32pk& b)
	{
		return _mm256_cmp_ps(a.v, b.v, _CMP_EQ_OQ);
	}

	LMAT_ENSURE_INLINE
	inline avx_f64bpk operator == (const avx_f64pk& a, const avx_f64pk& b)
	{
		return _mm256_cmp_pd(a.v, b.v, _CMP_EQ_OQ);
	}

	LMAT_ENSURE_INLINE
	inline avx_f32bpk operator != (const avx_f32pk& a, const avx_f32pk& b)
	{
		return _mm256_cmp_ps(a.v, b.v, _CMP_NEQ_OQ);
	}

	LMAT_ENSURE_INLINE
	inline avx_f64bpk operator != (const avx_f64pk& a, const avx_f64pk& b)
	{
		return _mm256_cmp_pd(a.v, b.v, _CMP_NEQ_OQ);
	}

	LMAT_ENSURE_INLINE
	inline avx_f32bpk operator > (const avx_f32pk& a, const avx_f32pk& b)
	{
		return _mm256_cmp_ps(a.v, b.v, _CMP_GT_OQ);
	}

	LMAT_ENSURE_INLINE
	inline avx_f64bpk operator > (const avx_f64pk& a, const avx_f64pk& b)
	{
		return _mm256_cmp_pd(a.v, b.v, _CMP_GT_OQ);
	}

	LMAT_ENSURE_INLINE
	inline avx_f32bpk operator >= (const avx_f32pk& a, const avx_f32pk& b)
	{
		return _mm256_cmp_ps(a.v, b.v, _CMP_GE_OQ);
	}

	LMAT_ENSURE_INLINE
	inline avx_f64bpk operator >= (const avx_f64pk& a, const avx_f64pk& b)
	{
		return _mm256_cmp_pd(a.v, b.v, _CMP_GE_OQ);
	}

	LMAT_ENSURE_INLINE
	inline avx_f32bpk operator < (const avx_f32pk& a, const avx_f32pk& b)
	{
		return _mm256_cmp_ps(a.v, b.v, _CMP_LT_OQ);
	}

	LMAT_ENSURE_INLINE
	inline avx_f64bpk operator < (const avx_f64pk& a, const avx_f64pk& b)
	{
		return _mm256_cmp_pd(a.v, b.v, _CMP_LT_OQ);
	}


	LMAT_ENSURE_INLINE
	inline avx_f32bpk operator <= (const avx_f32pk& a, const avx_f32pk& b)
	{
		return _mm256_cmp_ps(a.v, b.v, _CMP_LE_OQ);
	}

	LMAT_ENSURE_INLINE
	inline avx_f64bpk operator <= (const avx_f64pk& a, const avx_f64pk& b)
	{
		return _mm256_cmp_pd(a.v, b.v, _CMP_LE_OQ);
	}


	/********************************************
	 *
	 *  logical operations
	 *
	 ********************************************/

	LMAT_ENSURE_INLINE
	inline avx_f32bpk operator ~ (const avx_f32bpk& a)
	{
		return _mm256_xor_ps(a.v, _mm256_castsi256_ps(_mm256_set1_epi32(-1)));
	}

	LMAT_ENSURE_INLINE
	inline avx_f64bpk operator ~ (const avx_f64bpk& a)
	{
		return _mm256_xor_pd(a.v, _mm256_castsi256_pd(_mm256_set1_epi32(-1)));
	}

	LMAT_ENSURE_INLINE
	inline avx_f32bpk operator & (const avx_f32bpk& a, const avx_f32bpk& b)
	{
		return _mm256_and_ps(a.v, b.v);
	}

	LMAT_ENSURE_INLINE
	inline avx_f64bpk operator & (const avx_f64bpk& a, const avx_f64bpk& b)
	{
		return _mm256_and_pd(a.v, b.v);
	}

	LMAT_ENSURE_INLINE
	inline avx_f32bpk operator | (const avx_f32bpk& a, const avx_f32bpk& b)
	{
		return _mm256_or_ps(a.v, b.v);
	}

	LMAT_ENSURE_INLINE
	inline avx_f64bpk operator | (const avx_f64bpk& a, const avx_f64bpk& b)
	{
		return _mm256_or_pd(a.v, b.v);
	}

	LMAT_ENSURE_INLINE
	inline avx_f32bpk operator != (const avx_f32bpk& a, const avx_f32bpk& b)
	{
		return _mm256_xor_ps(a.v, b.v);
	}

	LMAT_ENSURE_INLINE
	inline avx_f64bpk operator != (const avx_f64bpk& a, const avx_f64bpk& b)
	{
		return _mm256_xor_pd(a.v, b.v);
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


	/********************************************
	 *
	 *  FP classification
	 *
	 ********************************************/

	LMAT_ENSURE_INLINE
	inline avx_f32bpk is_neg(const avx_f32pk& a)
	{
		return internal::combine_m128(
				internal::sse_is_neg_ps(a.get_low()),
				internal::sse_is_neg_ps(a.get_high()));
	}

	LMAT_ENSURE_INLINE
	inline avx_f64bpk is_neg(const avx_f64pk& a)
	{
		return internal::combine_m128d(
				internal::sse_is_neg_pd(a.get_low()),
				internal::sse_is_neg_pd(a.get_high()));
	}


	LMAT_ENSURE_INLINE
	inline avx_f32bpk is_finite(const avx_f32pk& a)
	{
		return internal::combine_m128(
				internal::sse_is_finite_ps(a.get_low()),
				internal::sse_is_finite_ps(a.get_high()));
	}

	LMAT_ENSURE_INLINE
	inline avx_f64bpk is_finite(const avx_f64pk& a)
	{
		return internal::combine_m128d(
				internal::sse_is_finite_pd(a.get_low()),
				internal::sse_is_finite_pd(a.get_high()));
	}


	LMAT_ENSURE_INLINE
	inline avx_f32bpk is_inf(const avx_f32pk& a)
	{
		return internal::combine_m128(
				internal::sse_is_inf_ps(a.get_low()),
				internal::sse_is_inf_ps(a.get_high()));
	}

	LMAT_ENSURE_INLINE
	inline avx_f64bpk is_inf(const avx_f64pk& a)
	{
		return internal::combine_m128d(
				internal::sse_is_inf_pd(a.get_low()),
				internal::sse_is_inf_pd(a.get_high()));
	}


	LMAT_ENSURE_INLINE
	inline avx_f32bpk is_nan(const avx_f32pk& a)
	{
		return internal::combine_m128(
				internal::sse_is_nan_ps(a.get_low()),
				internal::sse_is_nan_ps(a.get_high()));
	}

	LMAT_ENSURE_INLINE
	inline avx_f64bpk is_nan(const avx_f64pk& a)
	{
		return internal::combine_m128d(
				internal::sse_is_nan_pd(a.get_low()),
				internal::sse_is_nan_pd(a.get_high()));
	}


} }

#endif 
