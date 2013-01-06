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


#define LMAT_DEFINE_HAS_AVX_SUPPORT( FTag ) \
	template<> struct meta::has_simd_support<FTag, float, avx_t> { static const bool value = true; }; \
	template<> struct meta::has_simd_support<FTag, double, avx_t> { static const bool value = true; };

namespace lmat
{
	LMAT_DEFINE_HAS_AVX_SUPPORT( add_ )
	LMAT_DEFINE_HAS_AVX_SUPPORT( sub_ )
	LMAT_DEFINE_HAS_AVX_SUPPORT( mul_ )
	LMAT_DEFINE_HAS_AVX_SUPPORT( div_ )
	LMAT_DEFINE_HAS_AVX_SUPPORT( neg_ )

	LMAT_DEFINE_HAS_AVX_SUPPORT( eq_ )
	LMAT_DEFINE_HAS_AVX_SUPPORT( ne_ )
	LMAT_DEFINE_HAS_AVX_SUPPORT( gt_ )
	LMAT_DEFINE_HAS_AVX_SUPPORT( ge_ )
	LMAT_DEFINE_HAS_AVX_SUPPORT( lt_ )
	LMAT_DEFINE_HAS_AVX_SUPPORT( le_ )

	LMAT_DEFINE_HAS_AVX_SUPPORT( logical_and_ )
	LMAT_DEFINE_HAS_AVX_SUPPORT( logical_or_ )
	LMAT_DEFINE_HAS_AVX_SUPPORT( logical_eq_ )
	LMAT_DEFINE_HAS_AVX_SUPPORT( logical_ne_ )

	LMAT_DEFINE_HAS_AVX_SUPPORT( abs_ )
	LMAT_DEFINE_HAS_AVX_SUPPORT( fma_ )
	LMAT_DEFINE_HAS_AVX_SUPPORT( sqr_ )
	LMAT_DEFINE_HAS_AVX_SUPPORT( cube_ )

	LMAT_DEFINE_HAS_AVX_SUPPORT( rcp_ )
	LMAT_DEFINE_HAS_AVX_SUPPORT( sqrt_ )
	LMAT_DEFINE_HAS_AVX_SUPPORT( rsqrt_ )

	LMAT_DEFINE_HAS_AVX_SUPPORT( min_ )
	LMAT_DEFINE_HAS_AVX_SUPPORT( max_ )
	LMAT_DEFINE_HAS_AVX_SUPPORT( clamp_ )
	LMAT_DEFINE_HAS_AVX_SUPPORT( cond_ )

	LMAT_DEFINE_HAS_AVX_SUPPORT( floor_ )
	LMAT_DEFINE_HAS_AVX_SUPPORT( ceil_ )
	LMAT_DEFINE_HAS_AVX_SUPPORT( round_ )
	LMAT_DEFINE_HAS_AVX_SUPPORT( trunc_ )

	LMAT_DEFINE_HAS_AVX_SUPPORT( signbit_ )
	LMAT_DEFINE_HAS_AVX_SUPPORT( isfinite_ )
	LMAT_DEFINE_HAS_AVX_SUPPORT( isinf_ )
	LMAT_DEFINE_HAS_AVX_SUPPORT( isnan_ )
}

namespace lmat { namespace math {

	/********************************************
	 *
	 *  Floating-point arithmetics
	 *
	 ********************************************/

	LMAT_ENSURE_INLINE
	inline avx_f32pk operator + (const avx_f32pk& a, const avx_f32pk& b)
	{
		return _mm256_add_ps(a, b);
	}

	LMAT_ENSURE_INLINE
	inline avx_f64pk operator + (const avx_f64pk& a, const avx_f64pk& b)
	{
		return _mm256_add_pd(a, b);
	}

	LMAT_ENSURE_INLINE
	inline avx_f32pk operator - (const avx_f32pk& a, const avx_f32pk& b)
	{
		return _mm256_sub_ps(a, b);
	}

	LMAT_ENSURE_INLINE
	inline avx_f64pk operator - (const avx_f64pk& a, const avx_f64pk& b)
	{
		return _mm256_sub_pd(a, b);
	}

	LMAT_ENSURE_INLINE
	inline avx_f32pk operator * (const avx_f32pk& a, const avx_f32pk& b)
	{
		return _mm256_mul_ps(a, b);
	}

	LMAT_ENSURE_INLINE
	inline avx_f64pk operator * (const avx_f64pk& a, const avx_f64pk& b)
	{
		return _mm256_mul_pd(a, b);
	}

	LMAT_ENSURE_INLINE
	inline avx_f32pk operator / (const avx_f32pk& a, const avx_f32pk& b)
	{
		return _mm256_div_ps(a, b);
	}

	LMAT_ENSURE_INLINE
	inline avx_f64pk operator / (const avx_f64pk& a, const avx_f64pk& b)
	{
		return _mm256_div_pd(a, b);
	}

	LMAT_ENSURE_INLINE
	inline avx_f32pk operator - (const avx_f32pk& a)
	{
		typedef internal::num_fmt<float> fmt;
		return _mm256_xor_ps(_mm256_castsi256_ps(_mm256_set1_epi32(fmt::sign_bit)), a);
	}

	LMAT_ENSURE_INLINE
	inline avx_f64pk operator - (const avx_f64pk& a)
	{
		typedef internal::num_fmt<double> fmt;
		return _mm256_xor_pd(_mm256_castsi256_pd(_mm256_set1_epi64x(fmt::sign_bit)), a);
	}

	LMAT_ENSURE_INLINE
	inline avx_f32pk abs(const avx_f32pk& a)
	{
		typedef internal::num_fmt<float> fmt;
		return _mm256_andnot_ps(_mm256_castsi256_ps(_mm256_set1_epi32(fmt::sign_bit)), a);
	}

	LMAT_ENSURE_INLINE
	inline avx_f64pk abs(const avx_f64pk& a)
	{
		typedef internal::num_fmt<double> fmt;
		return _mm256_andnot_pd(_mm256_castsi256_pd(_mm256_set1_epi64x(fmt::sign_bit)), a);
	}

	LMAT_ENSURE_INLINE
	inline avx_f32pk fma(const avx_f32pk& x, const avx_f32pk& y, const avx_f32pk& z)
	{
		return _mm256_add_ps(_mm256_mul_ps(x, y), z);
	}

	LMAT_ENSURE_INLINE
	inline avx_f64pk fma(const avx_f64pk& x, const avx_f64pk& y, const avx_f64pk& z)
	{
		return _mm256_add_pd(_mm256_mul_pd(x, y), z);
	}

	LMAT_ENSURE_INLINE
	inline avx_f32pk& operator += (avx_f32pk& a, const avx_f32pk& b)
	{
		a = _mm256_add_ps(a, b);
		return a;
	}

	LMAT_ENSURE_INLINE
	inline avx_f64pk& operator += (avx_f64pk& a, const avx_f64pk& b)
	{
		a = _mm256_add_pd(a, b);
		return a;
	}

	LMAT_ENSURE_INLINE
	inline avx_f32pk& operator -= (avx_f32pk& a, const avx_f32pk& b)
	{
		a = _mm256_sub_ps(a, b);
		return a;
	}

	LMAT_ENSURE_INLINE
	inline avx_f64pk& operator -= (avx_f64pk& a, const avx_f64pk& b)
	{
		a = _mm256_sub_pd(a, b);
		return a;
	}

	LMAT_ENSURE_INLINE
	inline avx_f32pk& operator *= (avx_f32pk& a, const avx_f32pk& b)
	{
		a = _mm256_mul_ps(a, b);
		return a;
	}

	LMAT_ENSURE_INLINE
	inline avx_f64pk& operator *= (avx_f64pk& a, const avx_f64pk& b)
	{
		a = _mm256_mul_pd(a, b);
		return a;
	}

	LMAT_ENSURE_INLINE
	inline avx_f32pk& operator /= (avx_f32pk& a, const avx_f32pk& b)
	{
		a = _mm256_div_ps(a, b);
		return a;
	}

	LMAT_ENSURE_INLINE
	inline avx_f64pk& operator /= (avx_f64pk& a, const avx_f64pk& b)
	{
		a = _mm256_div_pd(a, b);
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
		return _mm256_min_ps(a, b);
	}

	LMAT_ENSURE_INLINE
	inline avx_f64pk (min)(const avx_f64pk& a, const avx_f64pk& b)
	{
		return _mm256_min_pd(a, b);
	}

	LMAT_ENSURE_INLINE
	inline avx_f32pk (max)(const avx_f32pk& a, const avx_f32pk& b)
	{
		return _mm256_max_ps(a, b);
	}

	LMAT_ENSURE_INLINE
	inline avx_f64pk (max)(const avx_f64pk& a, const avx_f64pk& b)
	{
		return _mm256_max_pd(a, b);
	}

	LMAT_ENSURE_INLINE
	inline avx_f32pk clamp(const avx_f32pk& x, const avx_f32pk& lb, const avx_f32pk& ub)
	{
		return (min)((max)(x, lb), ub);
	}

	LMAT_ENSURE_INLINE
	inline avx_f64pk clamp(const avx_f64pk& x, const avx_f64pk& lb, const avx_f64pk& ub)
	{
		return (min)((max)(x, lb), ub);
	}

	/********************************************
	 *
	 *  Simple power functions
	 *
	 ********************************************/

	LMAT_ENSURE_INLINE
	inline avx_f32pk sqr(const avx_f32pk& a)
	{
		return _mm256_mul_ps(a, a);
	}

	LMAT_ENSURE_INLINE
	inline avx_f64pk sqr(const avx_f64pk& a)
	{
		return _mm256_mul_pd(a, a);
	}

	LMAT_ENSURE_INLINE
	inline avx_f32pk diff_sqr(const avx_f32pk& a, const avx_f32pk& b)
	{
		return sqr(a - b);
	}

	LMAT_ENSURE_INLINE
	inline avx_f64pk diff_sqr(const avx_f64pk& a, const avx_f64pk& b)
	{
		return sqr(a - b);
	}

	LMAT_ENSURE_INLINE
	inline avx_f32pk diff_abs(const avx_f32pk& a, const avx_f32pk& b)
	{
		return abs(a - b);
	}

	LMAT_ENSURE_INLINE
	inline avx_f64pk diff_abs(const avx_f64pk& a, const avx_f64pk& b)
	{
		return abs(a - b);
	}

	LMAT_ENSURE_INLINE
	inline avx_f32pk cube(const avx_f32pk& a)
	{
		return _mm256_mul_ps(_mm256_mul_ps(a, a), a);
	}

	LMAT_ENSURE_INLINE
	inline avx_f64pk cube(const avx_f64pk& a)
	{
		return _mm256_mul_pd(_mm256_mul_pd(a, a), a);
	}

	LMAT_ENSURE_INLINE
	inline avx_f32pk sqrt(const avx_f32pk& a)
	{
		return _mm256_sqrt_ps(a);
	}

	LMAT_ENSURE_INLINE
	inline avx_f64pk sqrt(const avx_f64pk& a)
	{
		return _mm256_sqrt_pd(a);
	}

	LMAT_ENSURE_INLINE
	inline avx_f32pk rcp(const avx_f32pk& a)
	{
		return _mm256_div_ps(_mm256_set1_ps(1.0f), a);
	}

	LMAT_ENSURE_INLINE
	inline avx_f32pk approx_rcp(const avx_f32pk& a)
	{
		return _mm256_rcp_ps(a);
	}

	LMAT_ENSURE_INLINE
	inline avx_f64pk rcp(const avx_f64pk& a)
	{
		return _mm256_div_pd(_mm256_set1_pd(1.0), a);
	}

	LMAT_ENSURE_INLINE
	inline avx_f32pk rsqrt(const avx_f32pk& a)
	{
		return _mm256_div_ps(_mm256_set1_ps(1.0f), _mm256_sqrt_ps(a));
	}

	LMAT_ENSURE_INLINE
	inline avx_f32pk approx_rsqrt(const avx_f32pk& a)
	{
		return _mm256_rsqrt_ps(a);
	}

	LMAT_ENSURE_INLINE
	inline avx_f64pk rsqrt(const avx_f64pk& a)
	{
		return _mm256_div_pd(_mm256_set1_pd(1.0), _mm256_sqrt_pd(a));
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


	/********************************************
	 *
	 *  conditional
	 *
	 ********************************************/

	LMAT_ENSURE_INLINE
	inline avx_f32pk cond(const avx_f32bpk& b, const avx_f32pk& x, const avx_f32pk& y)
	{
		return _mm256_blendv_ps(y, x, b);
	}

	LMAT_ENSURE_INLINE
	inline avx_f64pk cond(const avx_f64bpk& b, const avx_f64pk& x, const avx_f64pk& y)
	{
		return _mm256_blendv_pd(y, x, b);
	}


	/********************************************
	 *
	 *  rounding
	 *
	 ********************************************/

	LMAT_ENSURE_INLINE
	inline avx_f32pk round(const avx_f32pk& a)
	{
		return _mm256_round_ps(a, 0);
	}

	LMAT_ENSURE_INLINE
	inline avx_f64pk round(const avx_f64pk& a)
	{
		return _mm256_round_pd(a, 0);
	}

	LMAT_ENSURE_INLINE
	inline avx_f32pk floor(const avx_f32pk& a)
	{
		return _mm256_round_ps(a, 1);
	}

	LMAT_ENSURE_INLINE
	inline avx_f64pk floor(const avx_f64pk& a)
	{
		return _mm256_round_pd(a, 1);
	}

	LMAT_ENSURE_INLINE
	inline avx_f32pk ceil(const avx_f32pk& a)
	{
		return _mm256_round_ps(a, 2);
	}

	LMAT_ENSURE_INLINE
	inline avx_f64pk ceil(const avx_f64pk& a)
	{
		return _mm256_round_pd(a, 2);
	}

	LMAT_ENSURE_INLINE
	inline avx_f32pk trunc(const avx_f32pk& a)
	{
		return _mm256_round_ps(a, 3);
	}

	LMAT_ENSURE_INLINE
	inline avx_f64pk trunc(const avx_f64pk& a)
	{
		return _mm256_round_pd(a, 3);
	}


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
