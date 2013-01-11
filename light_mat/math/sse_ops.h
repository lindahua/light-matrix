/**
 * @file sse_ops.h
 *
 * @brief Basic operations on SSE packs
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_SSE_OPS_H_
#define LIGHTMAT_SSE_OPS_H_

#include <light_mat/math/sse_packs.h>
#include <light_mat/math/sse_bpacks.h>
#include "internal/sse2_round_impl.h"


#define LMAT_DEFINE_HAS_SSE_SUPPORT( FTag ) \
	template<> struct has_simd_support<ftags::FTag, float, sse_t> : public true_ { }; \
	template<> struct has_simd_support<ftags::FTag, double, sse_t> : public true_ { };

namespace lmat { namespace meta {

	// arithmetics

	LMAT_DEFINE_HAS_SSE_SUPPORT( add_ )
	LMAT_DEFINE_HAS_SSE_SUPPORT( sub_ )
	LMAT_DEFINE_HAS_SSE_SUPPORT( mul_ )
	LMAT_DEFINE_HAS_SSE_SUPPORT( div_ )
	LMAT_DEFINE_HAS_SSE_SUPPORT( neg_ )
	LMAT_DEFINE_HAS_SSE_SUPPORT( fma_ )

	LMAT_DEFINE_HAS_SSE_SUPPORT( min_ )
	LMAT_DEFINE_HAS_SSE_SUPPORT( max_ )
	LMAT_DEFINE_HAS_SSE_SUPPORT( clamp_ )
	LMAT_DEFINE_HAS_SSE_SUPPORT( cond_ )

	// comparison

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

	// simple power functions

	LMAT_DEFINE_HAS_SSE_SUPPORT( abs_ )
	LMAT_DEFINE_HAS_SSE_SUPPORT( sqr_ )
	LMAT_DEFINE_HAS_SSE_SUPPORT( cube_ )

	LMAT_DEFINE_HAS_SSE_SUPPORT( rcp_ )
	LMAT_DEFINE_HAS_SSE_SUPPORT( sqrt_ )
	LMAT_DEFINE_HAS_SSE_SUPPORT( rsqrt_ )

	// rounding

	LMAT_DEFINE_HAS_SSE_SUPPORT( floor_ )
	LMAT_DEFINE_HAS_SSE_SUPPORT( ceil_ )
	LMAT_DEFINE_HAS_SSE_SUPPORT( round_ )
	LMAT_DEFINE_HAS_SSE_SUPPORT( trunc_ )

	// numeric predicates

	LMAT_DEFINE_HAS_SSE_SUPPORT( signbit_ )
	LMAT_DEFINE_HAS_SSE_SUPPORT( isfinite_ )
	LMAT_DEFINE_HAS_SSE_SUPPORT( isinf_ )
	LMAT_DEFINE_HAS_SSE_SUPPORT( isnan_ )

} }

namespace lmat { namespace math {


	/********************************************
	 *
	 *  Floating-point arithmetics
	 *
	 ********************************************/

	LMAT_ENSURE_INLINE
	inline sse_f32pk operator + (const sse_f32pk& a, const sse_f32pk& b)
	{
		return _mm_add_ps(a, b);
	}

	LMAT_ENSURE_INLINE
	inline sse_f64pk operator + (const sse_f64pk& a, const sse_f64pk& b)
	{
		return _mm_add_pd(a, b);
	}

	LMAT_ENSURE_INLINE
	inline sse_f32pk operator - (const sse_f32pk& a, const sse_f32pk& b)
	{
		return _mm_sub_ps(a, b);
	}

	LMAT_ENSURE_INLINE
	inline sse_f64pk operator - (const sse_f64pk& a, const sse_f64pk& b)
	{
		return _mm_sub_pd(a, b);
	}

	LMAT_ENSURE_INLINE
	inline sse_f32pk operator * (const sse_f32pk& a, const sse_f32pk& b)
	{
		return _mm_mul_ps(a, b);
	}

	LMAT_ENSURE_INLINE
	inline sse_f64pk operator * (const sse_f64pk& a, const sse_f64pk& b)
	{
		return _mm_mul_pd(a, b);
	}

	LMAT_ENSURE_INLINE
	inline sse_f32pk operator / (const sse_f32pk& a, const sse_f32pk& b)
	{
		return _mm_div_ps(a, b);
	}

	LMAT_ENSURE_INLINE
	inline sse_f64pk operator / (const sse_f64pk& a, const sse_f64pk& b)
	{
		return _mm_div_pd(a, b);
	}

	LMAT_ENSURE_INLINE
	inline sse_f32pk operator - (const sse_f32pk& a)
	{
		return _mm_xor_ps(
				_mm_castsi128_ps(internal::sse_signmask_ps()), a);
	}

	LMAT_ENSURE_INLINE
	inline sse_f64pk operator - (const sse_f64pk& a)
	{
		return _mm_xor_pd(
				_mm_castsi128_pd(internal::sse_signmask_pd()), a);
	}

	LMAT_ENSURE_INLINE
	inline sse_f32pk abs(const sse_f32pk& a)
	{
		return _mm_andnot_ps(
				_mm_castsi128_ps(internal::sse_signmask_ps()), a);
	}

	LMAT_ENSURE_INLINE
	inline sse_f64pk abs(const sse_f64pk& a)
	{
		return _mm_andnot_pd(
				_mm_castsi128_pd(internal::sse_signmask_pd()), a);
	}

	LMAT_ENSURE_INLINE
	inline sse_f32pk fma(const sse_f32pk& x, const sse_f32pk& y, const sse_f32pk& z)
	{
		return _mm_add_ps(_mm_mul_ps(x, y), z);
	}

	LMAT_ENSURE_INLINE
	inline sse_f64pk fma(const sse_f64pk& x, const sse_f64pk& y, const sse_f64pk& z)
	{
		return _mm_add_pd(_mm_mul_pd(x, y), z);
	}


	LMAT_ENSURE_INLINE
	inline sse_f32pk& operator += (sse_f32pk& a, const sse_f32pk& b)
	{
		a = _mm_add_ps(a, b);
		return a;
	}

	LMAT_ENSURE_INLINE
	inline sse_f64pk& operator += (sse_f64pk& a, const sse_f64pk& b)
	{
		a = _mm_add_pd(a, b);
		return a;
	}

	LMAT_ENSURE_INLINE
	inline sse_f32pk& operator -= (sse_f32pk& a, const sse_f32pk& b)
	{
		a = _mm_sub_ps(a, b);
		return a;
	}

	LMAT_ENSURE_INLINE
	inline sse_f64pk& operator -= (sse_f64pk& a, const sse_f64pk& b)
	{
		a = _mm_sub_pd(a, b);
		return a;
	}

	LMAT_ENSURE_INLINE
	inline sse_f32pk& operator *= (sse_f32pk& a, const sse_f32pk& b)
	{
		a = _mm_mul_ps(a, b);
		return a;
	}

	LMAT_ENSURE_INLINE
	inline sse_f64pk& operator *= (sse_f64pk& a, const sse_f64pk& b)
	{
		a = _mm_mul_pd(a, b);
		return a;
	}

	LMAT_ENSURE_INLINE
	inline sse_f32pk& operator /= (sse_f32pk& a, const sse_f32pk& b)
	{
		a = _mm_div_ps(a, b);
		return a;
	}

	LMAT_ENSURE_INLINE
	inline sse_f64pk& operator /= (sse_f64pk& a, const sse_f64pk& b)
	{
		a = _mm_div_pd(a, b);
		return a;
	}


	/********************************************
	 *
	 *  Floating-point min and max
	 *
	 ********************************************/

	LMAT_ENSURE_INLINE
	inline sse_f32pk (min)(const sse_f32pk& a, const sse_f32pk& b)
	{
		return _mm_min_ps(a, b);
	}

	LMAT_ENSURE_INLINE
	inline sse_f64pk (min)(const sse_f64pk& a, const sse_f64pk& b)
	{
		return _mm_min_pd(a, b);
	}

	LMAT_ENSURE_INLINE
	inline sse_f32pk (max)(const sse_f32pk& a, const sse_f32pk& b)
	{
		return _mm_max_ps(a, b);
	}
	LMAT_ENSURE_INLINE
	inline sse_f64pk (max)(const sse_f64pk& a, const sse_f64pk& b)
	{
		return _mm_max_pd(a, b);
	}

	LMAT_ENSURE_INLINE
	inline sse_f32pk clamp(const sse_f32pk& x, const sse_f32pk& lb, const sse_f32pk& ub)
	{
		return (min)((max)(x, lb), ub);
	}

	LMAT_ENSURE_INLINE
	inline sse_f64pk clamp(const sse_f64pk& x, const sse_f64pk& lb, const sse_f64pk& ub)
	{
		return (min)((max)(x, lb), ub);
	}

	/********************************************
	 *
	 *  Simple power functions
	 *
	 ********************************************/

	LMAT_ENSURE_INLINE
	inline sse_f32pk sqr(const sse_f32pk& a)
	{
		return _mm_mul_ps(a, a);
	}

	LMAT_ENSURE_INLINE
	inline sse_f64pk sqr(const sse_f64pk& a)
	{
		return _mm_mul_pd(a, a);
	}

	LMAT_ENSURE_INLINE
	inline sse_f32pk diff_sqr(const sse_f32pk& a, const sse_f32pk& b)
	{
		return sqr(a - b);
	}

	LMAT_ENSURE_INLINE
	inline sse_f64pk diff_sqr(const sse_f64pk& a, const sse_f64pk& b)
	{
		return sqr(a - b);
	}

	LMAT_ENSURE_INLINE
	inline sse_f32pk diff_abs(const sse_f32pk& a, const sse_f32pk& b)
	{
		return abs(a - b);
	}

	LMAT_ENSURE_INLINE
	inline sse_f64pk diff_abs(const sse_f64pk& a, const sse_f64pk& b)
	{
		return abs(a - b);
	}

	LMAT_ENSURE_INLINE
	inline sse_f32pk cube(const sse_f32pk& a)
	{
		return _mm_mul_ps(_mm_mul_ps(a, a), a);
	}

	LMAT_ENSURE_INLINE
	inline sse_f64pk cube(const sse_f64pk& a)
	{
		return _mm_mul_pd(_mm_mul_pd(a, a), a);
	}

	LMAT_ENSURE_INLINE
	inline sse_f32pk sqrt(const sse_f32pk& a)
	{
		return _mm_sqrt_ps(a);
	}

	LMAT_ENSURE_INLINE
	inline sse_f64pk sqrt(const sse_f64pk& a)
	{
		return _mm_sqrt_pd(a);
	}

	LMAT_ENSURE_INLINE
	inline sse_f32pk rcp(const sse_f32pk& a)
	{
		return _mm_div_ps(_mm_set1_ps(1.0f), a);
	}

	LMAT_ENSURE_INLINE
	inline sse_f32pk approx_rcp(const sse_f32pk& a)
	{
		return _mm_rcp_ps(a);
	}

	LMAT_ENSURE_INLINE
	inline sse_f64pk rcp(const sse_f64pk& a)
	{
		return _mm_div_pd(_mm_set1_pd(1.0), a);
	}

	LMAT_ENSURE_INLINE
	inline sse_f32pk rsqrt(const sse_f32pk& a)
	{
		return _mm_div_ps(_mm_set1_ps(1.0f), _mm_sqrt_ps(a));
	}

	LMAT_ENSURE_INLINE
	inline sse_f32pk approx_rsqrt(const sse_f32pk& a)
	{
		return _mm_rsqrt_ps(a);
	}

	LMAT_ENSURE_INLINE
	inline sse_f64pk rsqrt(const sse_f64pk& a)
	{
		return _mm_div_pd(_mm_set1_pd(1.0), _mm_sqrt_pd(a));
	}


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


	/********************************************
	 *
	 *  blending
	 *
	 ********************************************/

	namespace internal
	{
		LMAT_ENSURE_INLINE
		inline sse_f32pk cond_sse2(const sse_f32bpk& b, const sse_f32pk& x, const sse_f32pk& y)
		{
			return _mm_or_ps(_mm_and_ps(b, x),  _mm_and_ps(~b, y));
		}

		LMAT_ENSURE_INLINE
		inline sse_f64pk cond_sse2(const sse_f64bpk& b, const sse_f64pk& x, const sse_f64pk& y)
		{
			return _mm_or_pd(_mm_and_pd(b, x),  _mm_and_pd(~b, y));
		}
	}

	LMAT_ENSURE_INLINE
	inline sse_f32pk cond(const sse_f32bpk& b, const sse_f32pk& x, const sse_f32pk& y)
	{
#ifdef LMAT_HAS_SSE4_1
		return _mm_blendv_ps(y, x, b);
#else
		return cond_sse2(b, x, y);
#endif
	}

	LMAT_ENSURE_INLINE
	inline sse_f64pk cond(const sse_f64bpk& b, const sse_f64pk& x, const sse_f64pk& y)
	{
#ifdef LMAT_HAS_SSE4_1
		return _mm_blendv_pd(y, x, b);
#else
		return cond_sse2(b, x, y);
#endif
	}



	/********************************************
	 *
	 *  rounding
	 *
	 ********************************************/

	LMAT_ENSURE_INLINE
	inline sse_f32pk round(const sse_f32pk& a)
	{
#ifdef LMAT_HAS_SSE4_1
		return _mm_round_ps(a, 0);
#else
		return internal::round_sse2(a);
#endif
	}

	LMAT_ENSURE_INLINE
	inline sse_f64pk round(const sse_f64pk& a)
	{
#ifdef LMAT_HAS_SSE4_1
		return _mm_round_pd(a, 0);
#else
		return internal::round_sse2(a);
#endif
	}

	LMAT_ENSURE_INLINE
	inline sse_f32pk floor(const sse_f32pk& a)
	{
#ifdef LMAT_HAS_SSE4_1
		return _mm_round_ps(a, 1);
#else
		return internal::floor_sse2(a);
#endif
	}

	LMAT_ENSURE_INLINE
	inline sse_f64pk floor(const sse_f64pk& a)
	{
#ifdef LMAT_HAS_SSE4_1
		return _mm_round_pd(a, 1);
#else
		return internal::floor_sse2(a);
#endif
	}

	LMAT_ENSURE_INLINE
	inline sse_f32pk ceil(const sse_f32pk& a)
	{
#ifdef LMAT_HAS_SSE4_1
		return _mm_round_ps(a, 2);
#else
		return internal::ceil_sse2(a);
#endif
	}

	LMAT_ENSURE_INLINE
	inline sse_f64pk ceil(const sse_f64pk& a)
	{
#ifdef LMAT_HAS_SSE4_1
		return _mm_round_pd(a, 2);
#else
		return internal::ceil_sse2(a);
#endif
	}

	LMAT_ENSURE_INLINE
	inline sse_f32pk trunc(const sse_f32pk& a)
	{
#ifdef LMAT_HAS_SSE4_1
		return _mm_round_ps(a, 3);
#else
		return internal::trunc_sse2(a);
#endif
	}

	LMAT_ENSURE_INLINE
	inline sse_f64pk trunc(const sse_f64pk& a)
	{
#ifdef LMAT_HAS_SSE4_1
		return _mm_round_pd(a, 3);
#else
		return internal::trunc_sse2(a);
#endif
	}



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

#endif /* SSE_ARITH_H_ */
