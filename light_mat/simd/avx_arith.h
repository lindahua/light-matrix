/**
 * @file avx_arith.h
 *
 * @brief Arithmetics on AVX packs
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_AVX_ARITH_H_
#define LIGHTMAT_AVX_ARITH_H_

#include <light_mat/simd/avx_packs.h>
#include <light_mat/simd/avx_bpacks.h>

namespace lmat { namespace meta {

	// arithmetics

	LMAT_DEFINE_HAS_AVX_SUPPORT( add_ )
	LMAT_DEFINE_HAS_AVX_SUPPORT( sub_ )
	LMAT_DEFINE_HAS_AVX_SUPPORT( mul_ )
	LMAT_DEFINE_HAS_AVX_SUPPORT( div_ )
	LMAT_DEFINE_HAS_AVX_SUPPORT( neg_ )
	LMAT_DEFINE_HAS_AVX_SUPPORT( fma_ )

	LMAT_DEFINE_HAS_AVX_SUPPORT( min_ )
	LMAT_DEFINE_HAS_AVX_SUPPORT( max_ )
	LMAT_DEFINE_HAS_AVX_SUPPORT( clamp_ )
	LMAT_DEFINE_HAS_AVX_SUPPORT( cond_ )

	// simple power functions

	LMAT_DEFINE_HAS_AVX_SUPPORT( abs_ )
	LMAT_DEFINE_HAS_AVX_SUPPORT( sqr_ )
	LMAT_DEFINE_HAS_AVX_SUPPORT( cube_ )

	LMAT_DEFINE_HAS_AVX_SUPPORT( rcp_ )
	LMAT_DEFINE_HAS_AVX_SUPPORT( sqrt_ )
	LMAT_DEFINE_HAS_AVX_SUPPORT( rsqrt_ )

	// rounding

	LMAT_DEFINE_HAS_AVX_SUPPORT( floor_ )
	LMAT_DEFINE_HAS_AVX_SUPPORT( ceil_ )
	LMAT_DEFINE_HAS_AVX_SUPPORT( round_ )
	LMAT_DEFINE_HAS_AVX_SUPPORT( trunc_ )

} }


namespace lmat
{

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

}


namespace lmat {  namespace math {

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

} }

#endif
