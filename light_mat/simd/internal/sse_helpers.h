/**
 * @file sse_helpers.h
 *
 * @brief Auxiliary devices for SSE operations
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_SSE_HELPERS_H_
#define LIGHTMAT_SSE_HELPERS_H_

#include <light_mat/simd/simd_base.h>
#include "numrepr_format.h"

namespace lmat { namespace internal {

	// bitwise not

	LMAT_ENSURE_INLINE
	inline __m128i sse_bitwise_not(const __m128i& a)
	{
		return _mm_xor_si128(a, _mm_set1_epi32(-1));
	}


	// partial load

	LMAT_ENSURE_INLINE
	inline __m128 sse_loadpart_f32(siz_<1>, const float *p)
	{
		return _mm_load_ss(p);
	}

	LMAT_ENSURE_INLINE
	inline __m128 sse_loadpart_f32(siz_<2>, const float *p)
	{
		return _mm_castpd_ps(_mm_load_sd((double*)p));
	}

	LMAT_ENSURE_INLINE
	inline __m128 sse_loadpart_f32(siz_<3>, const float *p)
	{
		return _mm_movelh_ps(
				_mm_castpd_ps(_mm_load_sd((double*)p)), _mm_load_ss(p + 2));
	}

	LMAT_ENSURE_INLINE
	inline __m128 sse_loadpart_f32(siz_<4>, const float *p)
	{
		return _mm_loadu_ps(p);
	}


	LMAT_ENSURE_INLINE
	inline __m128d sse_loadpart_f64(siz_<1>, const double *p)
	{
		return _mm_load_sd(p);
	}


	LMAT_ENSURE_INLINE
	inline __m128d sse_loadpart_f64(siz_<2>, const double *p)
	{
		return _mm_loadu_pd(p);
	}


	// partial store

	LMAT_ENSURE_INLINE
	inline void sse_storepart_f32(siz_<1>, float *p, const __m128& v)
	{
		_mm_store_ss(p, v);
	}

	LMAT_ENSURE_INLINE
	inline void sse_storepart_f32(siz_<2>, float *p, const __m128& v)
	{
		_mm_store_sd((double*)p, _mm_castps_pd(v));
	}

	LMAT_ENSURE_INLINE
	inline void sse_storepart_f32(siz_<3>, float *p, const __m128& v)
	{
		_mm_store_sd((double*)p, _mm_castps_pd(v));
		_mm_store_ss(p + 2, _mm_movehl_ps(v, v));
	}

	LMAT_ENSURE_INLINE
	inline void sse_storepart_f32(siz_<4>, float *p, const __m128& v)
	{
		_mm_storeu_ps(p, v);
	}

	LMAT_ENSURE_INLINE
	inline void sse_storepart_f64(siz_<1>, double *p, const __m128d& v)
	{
		_mm_store_sd(p, v);
	}

	LMAT_ENSURE_INLINE
	inline void sse_storepart_f64(siz_<2>, double *p, const __m128d& v)
	{
		_mm_storeu_pd(p, v);
	}


	// extract scalar

	LMAT_ENSURE_INLINE
	inline float sse_extract_f32(const __m128& v, pos_<0> )
	{
		return _mm_cvtss_f32(v);
	}

	LMAT_ENSURE_INLINE
	inline float sse_extract_f32(const __m128& v, pos_<1> )
	{
		return _mm_cvtss_f32(
				_mm_castsi128_ps(_mm_srli_si128(_mm_castps_si128(v), 4)));
	}

	LMAT_ENSURE_INLINE
	inline float sse_extract_f32(const __m128& v, pos_<2> )
	{
		return _mm_cvtss_f32(
				_mm_castsi128_ps(_mm_srli_si128(_mm_castps_si128(v), 8)));
	}

	LMAT_ENSURE_INLINE
	inline float sse_extract_f32(const __m128& v, pos_<3> )
	{
		return _mm_cvtss_f32(
				_mm_castsi128_ps(_mm_srli_si128(_mm_castps_si128(v), 12)));
	}

	LMAT_ENSURE_INLINE
	inline double sse_extract_f64(const __m128d& v, pos_<0> )
	{
		return _mm_cvtsd_f64(v);
	}

	LMAT_ENSURE_INLINE
	inline double sse_extract_f64(const __m128d& v, pos_<1> )
	{
		return _mm_cvtsd_f64(
				_mm_castsi128_pd(_mm_srli_si128(_mm_castpd_si128(v), 8)));
	}

	// broadcasting

    LMAT_ENSURE_INLINE __m128 sse_broadcast_f32(__m128 v, pos_<0>)
    {
    	return _mm_shuffle_ps(v, v, 0);
    }

    LMAT_ENSURE_INLINE __m128 sse_broadcast_f32(__m128 v, pos_<1>)
    {
    	return _mm_shuffle_ps(v, v, 0x55);
    }

    LMAT_ENSURE_INLINE __m128 sse_broadcast_f32(__m128 v, pos_<2>)
    {
    	return _mm_shuffle_ps(v, v, 0xaa);
    }

    LMAT_ENSURE_INLINE __m128 sse_broadcast_f32(__m128 v, pos_<3>)
    {
    	return _mm_shuffle_ps(v, v, 0xff);
    }

    LMAT_ENSURE_INLINE __m128d sse_broadcast_f64(__m128d v, pos_<0>)
    {
    	return _mm_unpacklo_pd(v, v);
    }

    LMAT_ENSURE_INLINE __m128d sse_broadcast_f64(__m128d v, pos_<1>)
    {
    	return _mm_unpackhi_pd(v, v);
    }



	// conditional operator (SSE2)

	LMAT_ENSURE_INLINE
	inline __m128 cond_sse2(const __m128& b, const __m128& x, const __m128& y)
	{
		return _mm_or_ps(_mm_and_ps(b, x),  _mm_andnot_ps(b, y));
	}

	LMAT_ENSURE_INLINE
	inline __m128d cond_sse2(const __m128d& b, const __m128d& x, const __m128d& y)
	{
		return _mm_or_pd(_mm_and_pd(b, x),  _mm_andnot_pd(b, y));
	}

	// sign-mask

	LMAT_ENSURE_INLINE
	inline __m128i sse_signmask_ps()
	{
		typedef num_fmt<float> fmt;
		return _mm_set1_epi32(fmt::sign_bit);
	}

	LMAT_ENSURE_INLINE
	inline __m128i sse_signmask_pd()
	{
		typedef num_fmt<double> fmt;
		return _mm_set1_epi64x(fmt::sign_bit);
	}
} }


#endif /* SSE_HELPERS_H_ */
