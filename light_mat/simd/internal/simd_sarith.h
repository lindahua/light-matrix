/**
 * @file simd_sarith.h
 *
 * @brief Single-scalar arithmetics on SIMD packs
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_SIMD_SARITH_H_
#define LIGHTMAT_SIMD_SARITH_H_

#include <light_mat/simd/simd_packs.h>
#include <light_mat/simd/sse_reduce.h>

namespace lmat { namespace internal {

	// SSE

	LMAT_ENSURE_INLINE
	inline sse_f32pk sadd(const sse_f32pk& a, const sse_f32pk& b)
	{
		return _mm_add_ps(a, b);
	}

	LMAT_ENSURE_INLINE
	inline sse_f32pk ssub(const sse_f32pk& a, const sse_f32pk& b)
	{
		return _mm_sub_ps(a, b);
	}

	LMAT_ENSURE_INLINE
	inline sse_f32pk smul(const sse_f32pk& a, const sse_f32pk& b)
	{
		return _mm_mul_ps(a, b);
	}

	LMAT_ENSURE_INLINE
	inline sse_f64pk sadd(const sse_f64pk& a, const sse_f64pk& b)
	{
		return _mm_add_pd(a, b);
	}

	LMAT_ENSURE_INLINE
	inline sse_f64pk ssub(const sse_f64pk& a, const sse_f64pk& b)
	{
		return _mm_sub_pd(a, b);
	}

	LMAT_ENSURE_INLINE
	inline sse_f64pk smul(const sse_f64pk& a, const sse_f64pk& b)
	{
		return _mm_mul_pd(a, b);
	}


	// AVX

	LMAT_ENSURE_INLINE
	inline avx_f32pk sadd(const avx_f32pk& a, const avx_f32pk& b)
	{
		return _mm256_add_ps(a, b);
	}

	LMAT_ENSURE_INLINE
	inline avx_f32pk ssub(const avx_f32pk& a, const avx_f32pk& b)
	{
		return _mm256_sub_ps(a, b);
	}

	LMAT_ENSURE_INLINE
	inline avx_f32pk smul(const avx_f32pk& a, const avx_f32pk& b)
	{
		return _mm256_mul_ps(a, b);
	}

	LMAT_ENSURE_INLINE
	inline avx_f64pk sadd(const avx_f64pk& a, const avx_f64pk& b)
	{
		return _mm256_add_pd(a, b);
	}

	LMAT_ENSURE_INLINE
	inline avx_f64pk ssub(const avx_f64pk& a, const avx_f64pk& b)
	{
		return _mm256_sub_pd(a, b);
	}

	LMAT_ENSURE_INLINE
	inline avx_f64pk smul(const avx_f64pk& a, const avx_f64pk& b)
	{
		return _mm256_mul_pd(a, b);
	}


	// Part sum

	LMAT_ENSURE_INLINE
	inline float sum_2(const sse_f32pk& a)
	{
		__m128 t = _mm_add_ss(a, _mm_shuffle_ps(a, a, 1));
		return _mm_cvtss_f32(t);
	}

	LMAT_ENSURE_INLINE
	inline float sum_2(const avx_f32pk& a)
	{
		return sum_2(sse_f32pk(a.get_low()));
	}

	LMAT_ENSURE_INLINE
	inline double sum_2(const avx_f64pk& a)
	{
		return sum(sse_f64pk(a.get_low()));
	}

	LMAT_ENSURE_INLINE
	inline float sum_3(const sse_f32pk& a)
	{
		__m128 t1 = _mm_add_ss(a, _mm_movehl_ps(a, a));
		__m128 t2 = _mm_add_ss(t1, _mm_shuffle_ps(t1, t1, 1));
		return _mm_cvtss_f32(t2);
	}

	LMAT_ENSURE_INLINE
	inline float sum_3(const avx_f32pk& a)
	{
		return sum_3(sse_f32pk(a.get_low()));
	}

	LMAT_ENSURE_INLINE
	inline double sum_3(const avx_f64pk& a)
	{
		__m128d t = _mm_add_sd(a.get_low(), a.get_high());
		return sum(sse_f64pk(t));
	}

	LMAT_ENSURE_INLINE
	inline float sum_4(const avx_f32pk& a)
	{
		return sum(sse_f32pk(a.get_low()));
	}


} }

#endif
