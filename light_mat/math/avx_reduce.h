/**
 * @file avx_reduce.h
 *
 * Reduction on AVX packs
 * 
 * @author Dahua Lin 
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_AVX_REDUCE_H_
#define LIGHTMAT_AVX_REDUCE_H_

#include <light_mat/math/avx_packs.h>
#include <light_mat/math/avx_bpacks.h>
#include <light_mat/math/sse_reduce.h>

namespace lmat { namespace math {

	// numeric reduction

	LMAT_ENSURE_INLINE
	inline float sum(const avx_f32pk& a)
	{
		sse_f32pk t = _mm_add_ps(a.get_low(), a.get_high());
		return sum(t);
	}

	LMAT_ENSURE_INLINE
	inline double sum(const avx_f64pk& a)
	{
		sse_f64pk t = _mm_add_pd(a.get_low(), a.get_high());
		return sum(t);
	}

	LMAT_ENSURE_INLINE
	inline float maximum(const avx_f32pk& a)
	{
		sse_f32pk t = _mm_max_ps(a.get_low(), a.get_high());
		return maximum(t);
	}

	LMAT_ENSURE_INLINE
	inline double maximum(const avx_f64pk& a)
	{
		sse_f64pk t = _mm_max_pd(a.get_low(), a.get_high());
		return maximum(t);
	}

	LMAT_ENSURE_INLINE
	inline float minimum(const avx_f32pk& a)
	{
		sse_f32pk t = _mm_min_ps(a.get_low(), a.get_high());
		return minimum(t);
	}

	LMAT_ENSURE_INLINE
	inline double minimum(const avx_f64pk& a)
	{
		sse_f64pk t = _mm_min_pd(a.get_low(), a.get_high());
		return minimum(t);
	}


	// all & any

	LMAT_ENSURE_INLINE
	inline bool all_true(const avx_f32bpk& a)
	{
		return _mm256_testc_ps(a, avx_f32bpk::all_true());
	}

	LMAT_ENSURE_INLINE
	inline bool all_true(const avx_f64bpk& a)
	{
		return _mm256_testc_pd(a, avx_f64bpk::all_true());
	}

	LMAT_ENSURE_INLINE
	inline bool all_false(const avx_f32bpk& a)
	{
		return _mm256_testz_ps(a, avx_f32bpk::all_true());
	}

	LMAT_ENSURE_INLINE
	inline bool all_false(const avx_f64bpk& a)
	{
		return _mm256_testz_pd(a, avx_f64bpk::all_true());
	}


	LMAT_ENSURE_INLINE
	inline bool any_true(const avx_f32bpk& a)
	{
		return !all_false(a);
	}

	LMAT_ENSURE_INLINE
	inline bool any_true(const avx_f64bpk& a)
	{
		return !all_false(a);
	}

	LMAT_ENSURE_INLINE
	inline bool any_false(const avx_f32bpk& a)
	{
		return !all_true(a);
	}

	LMAT_ENSURE_INLINE
	inline bool any_false(const avx_f64bpk& a)
	{
		return !all_true(a);
	}


} }

#endif 
