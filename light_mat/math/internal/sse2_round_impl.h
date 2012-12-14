/**
 * @file sse_round_impl.h
 *
 * @brief Internal implementation of SSE rounding
 *
 * @author Dahua Lin
 */

#ifndef LIGHTMAT_SSE_ROUND_IMPL_H_
#define LIGHTMAT_SSE_ROUND_IMPL_H_

#include <light_mat/math/simd_base.h>

namespace lmat {  namespace math { namespace internal {

	LMAT_ENSURE_INLINE
	inline sse_f32pk floor_sse2(const sse_f32pk& a)
	{
		__m128 t = _mm_cvtepi32_ps(_mm_cvttps_epi32(a.v));
		__m128 b = _mm_and_ps(_mm_cmpgt_ps(t, a.v), _mm_set1_ps(1.0f));

		return _mm_sub_ps(t, b);
	}

	LMAT_ENSURE_INLINE
	inline sse_f64pk floor_sse2(const sse_f64pk& a)
	{
		__m128d t = _mm_cvtepi32_pd(_mm_cvttpd_epi32(a.v));
		__m128d b = _mm_and_pd(_mm_cmpgt_pd(t, a.v), _mm_set1_pd(1.0));

		return _mm_sub_pd(t, b);
	}

	LMAT_ENSURE_INLINE
	inline sse_f32pk ceil_sse2(const sse_f32pk& a)
	{
		__m128 t = _mm_cvtepi32_ps(_mm_cvttps_epi32(a.v));
		__m128 b = _mm_and_ps(_mm_cmplt_ps(t, a.v), _mm_set1_ps(1.0f));

		return _mm_add_ps(t, b);
	}

	LMAT_ENSURE_INLINE
	inline sse_f64pk ceil_sse2(const sse_f64pk& a)
	{
		__m128d t = _mm_cvtepi32_pd(_mm_cvttpd_epi32(a.v));
		__m128d b = _mm_and_pd(_mm_cmplt_pd(t, a.v), _mm_set1_pd(1.0));

		return _mm_add_pd(t, b);
	}

	LMAT_ENSURE_INLINE
	inline sse_f32pk trunc_sse2(const sse_f32pk& a)
	{
		return _mm_cvtepi32_ps(_mm_cvttps_epi32(a.v));
	}

	LMAT_ENSURE_INLINE
	inline sse_f64pk trunc_sse2(const sse_f64pk& a)
	{
		return _mm_cvtepi32_pd(_mm_cvttpd_epi32(a.v));
	}


	// round: using magic-number method from Agner Fog.

	LMAT_ENSURE_INLINE
	inline sse_f32pk round_sse2(const sse_f32pk& a)
	{
		__m128 sign  = _mm_and_ps(a.v, _mm_castsi128_ps(_mm_set1_epi32((int)0x80000000)));
		__m128 magic = _mm_castsi128_ps(_mm_set1_epi32((int)0x4b000000));
		__m128 smagic = _mm_or_ps(magic, sign);

		return _mm_sub_ps(_mm_add_ps(a.v, smagic), smagic);
	}

	LMAT_ENSURE_INLINE
	inline sse_f64pk round_sse2(const sse_f64pk& a)
	{
		__m128d signmsk = _mm_castsi128_pd(_mm_setr_epi32(0, (int)0x80000000, 0, (int)0x80000000));
		__m128d sign    = _mm_and_pd(a.v, signmsk);
		__m128d magic   = _mm_castsi128_pd(_mm_setr_epi32(0, (int)0x43300000, 0, (int)0x43300000));
		__m128d smagic  = _mm_or_pd(magic, sign);

		return _mm_sub_pd(_mm_add_pd(a.v, smagic), smagic);
	}

} } }

#endif /* SSE_ROUND_IMPL_H_ */
