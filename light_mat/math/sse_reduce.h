/**
 * @file sse_reduce.h
 *
 * @brief SSE reduction
 *
 * @author Dahua Lin
 */

#ifndef LIGHTMAT_SSE_REDUCE_H_
#define LIGHTMAT_SSE_REDUCE_H_

#include <light_mat/math/sse_packs.h>
#include <light_mat/math/sse_bpacks.h>
#include "internal/sse_testz_impl.h"

namespace lmat { namespace math {

	// sum

	LMAT_ENSURE_INLINE
	inline float sum(const sse_f32pk& a)
	{
		__m128 t1 = _mm_add_ps(a, _mm_movehl_ps(a, a));
		__m128 t2 = _mm_add_ss(t1, _mm_shuffle_ps(t1, t1, 1));
		return _mm_cvtss_f32(t2);
	}

	LMAT_ENSURE_INLINE
	inline double sum(const sse_f64pk& a)
	{
		__m128d t1 = _mm_add_pd(a, _mm_unpackhi_pd(a, a));
		return _mm_cvtsd_f64(t1);
	}


	// max

	LMAT_ENSURE_INLINE
	inline float maximum(const sse_f32pk& a)
	{
		__m128 t1 = _mm_max_ps(a, _mm_movehl_ps(a, a));
		__m128 t2 = _mm_max_ss(t1, _mm_shuffle_ps(t1, t1, 1));
		return _mm_cvtss_f32(t2);
	}

	LMAT_ENSURE_INLINE
	inline double maximum(const sse_f64pk& a)
	{
		__m128d t1 = _mm_max_pd(a, _mm_unpackhi_pd(a, a));
		return _mm_cvtsd_f64(t1);
	}


	// min

	LMAT_ENSURE_INLINE
	inline float minimum(const sse_f32pk& a)
	{
		__m128 t1 = _mm_min_ps(a, _mm_movehl_ps(a, a));
		__m128 t2 = _mm_min_ss(t1, _mm_shuffle_ps(t1, t1, 1));
		return _mm_cvtss_f32(t2);
	}

	LMAT_ENSURE_INLINE
	inline double minimum(const sse_f64pk& a)
	{
		__m128d t1 = _mm_min_pd(a, _mm_unpackhi_pd(a, a));
		return _mm_cvtsd_f64(t1);
	}


	// all & any

	LMAT_ENSURE_INLINE
	inline bool all_true(const sse_f32bpk& a)
	{
		return internal::testc(_mm_castps_si128(a));
	}

	LMAT_ENSURE_INLINE
	inline bool all_true(const sse_f64bpk& a)
	{
		return internal::testc(_mm_castpd_si128(a));
	}

	LMAT_ENSURE_INLINE
	inline bool all_false(const sse_f32bpk& a)
	{
		return internal::testz(_mm_castps_si128(a));
	}

	LMAT_ENSURE_INLINE
	inline bool all_false(const sse_f64bpk& a)
	{
		return internal::testz(_mm_castpd_si128(a));
	}


	LMAT_ENSURE_INLINE
	inline bool any_true(const sse_f32bpk& a)
	{
		return !all_false(a);
	}

	LMAT_ENSURE_INLINE
	inline bool any_true(const sse_f64bpk& a)
	{
		return !all_false(a);
	}

	LMAT_ENSURE_INLINE
	inline bool any_false(const sse_f32bpk& a)
	{
		return !all_true(a);
	}

	LMAT_ENSURE_INLINE
	inline bool any_false(const sse_f64bpk& a)
	{
		return !all_true(a);
	}

} }

#endif /* SSE_REDUCE_H_ */
