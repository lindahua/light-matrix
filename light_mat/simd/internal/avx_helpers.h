/**
 * @file avx_helpers.h
 *
 * @brief
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_AVX_HELPERS_H_
#define LIGHTMAT_AVX_HELPERS_H_

#include "sse_helpers.h"

namespace lmat { namespace internal {

	LMAT_ENSURE_INLINE
	inline __m256 combine_m128(const __m128& lo, const __m128& hi)
	{
		return _mm256_insertf128_ps(_mm256_castps128_ps256(lo), hi, 1);
	}

	LMAT_ENSURE_INLINE
	inline __m256d combine_m128d(const __m128d& lo, const __m128d& hi)
	{
		return _mm256_insertf128_pd(_mm256_castpd128_pd256(lo), hi, 1);
	}


	// part mask

	LMAT_ENSURE_INLINE
	inline __m256i avx_part_mask_32(siz_<0> )
	{
		return _mm256_setzero_si256();
	}

	LMAT_ENSURE_INLINE
	inline __m256i avx_part_mask_32(siz_<1> )
	{
		return _mm256_setr_epi32(-1, 0, 0, 0, 0, 0, 0, 0);
	}

	LMAT_ENSURE_INLINE
	inline __m256i avx_part_mask_32(siz_<2> )
	{
		return _mm256_setr_epi32(-1, -1, 0, 0, 0, 0, 0, 0);
	}

	LMAT_ENSURE_INLINE
	inline __m256i avx_part_mask_32(siz_<3> )
	{
		return _mm256_setr_epi32(-1, -1, -1, 0, 0, 0, 0, 0);
	}

	LMAT_ENSURE_INLINE
	inline __m256i avx_part_mask_32(siz_<4> )
	{
		return _mm256_setr_epi32(-1, -1, -1, -1, 0, 0, 0, 0);
	}

	LMAT_ENSURE_INLINE
	inline __m256i avx_part_mask_32(siz_<5> )
	{
		return _mm256_setr_epi32(-1, -1, -1, -1, -1, 0, 0, 0);
	}

	LMAT_ENSURE_INLINE
	inline __m256i avx_part_mask_32(siz_<6> )
	{
		return _mm256_setr_epi32(-1, -1, -1, -1, -1, -1, 0, 0);
	}

	LMAT_ENSURE_INLINE
	inline __m256i avx_part_mask_32(siz_<7> )
	{
		return _mm256_setr_epi32(-1, -1, -1, -1, -1, -1, -1, 0);
	}

	LMAT_ENSURE_INLINE
	inline __m256i avx_part_mask_32(siz_<8> )
	{
		return _mm256_set1_epi32(-1);
	}


	LMAT_ENSURE_INLINE
	inline __m256i avx_part_mask_64(siz_<0> )
	{
		return _mm256_setzero_si256();
	}

	LMAT_ENSURE_INLINE
	inline __m256i avx_part_mask_64(siz_<1> )
	{
		return _mm256_setr_epi32(-1, -1, 0, 0, 0, 0, 0, 0);
	}

	LMAT_ENSURE_INLINE
	inline __m256i avx_part_mask_64(siz_<2> )
	{
		return _mm256_setr_epi32(-1, -1, -1, -1, 0, 0, 0, 0);
	}

	LMAT_ENSURE_INLINE
	inline __m256i avx_part_mask_64(siz_<3> )
	{
		return _mm256_setr_epi32(-1, -1, -1, -1, -1, -1, 0, 0);
	}

	LMAT_ENSURE_INLINE
	inline __m256i avx_part_mask_64(siz_<4> )
	{
		return _mm256_set1_epi32(-1);
	}


	// AVX extraction

	LMAT_ENSURE_INLINE
	inline float avx_extract_f32(__m256 v, pos_<0> p)
	{
		return sse_extract_f32(_mm256_castps256_ps128(v), pos_<0>());
	}

	LMAT_ENSURE_INLINE
	inline float avx_extract_f32(__m256 v, pos_<1> p)
	{
		return sse_extract_f32(_mm256_castps256_ps128(v), pos_<1>());
	}

	LMAT_ENSURE_INLINE
	inline float avx_extract_f32(__m256 v, pos_<2> p)
	{
		return sse_extract_f32(_mm256_castps256_ps128(v), pos_<2>());
	}

	LMAT_ENSURE_INLINE
	inline float avx_extract_f32(__m256 v, pos_<3> p)
	{
		return sse_extract_f32(_mm256_castps256_ps128(v), pos_<3>());
	}

	LMAT_ENSURE_INLINE
	inline float avx_extract_f32(__m256 v, pos_<4> p)
	{
		return sse_extract_f32(_mm256_extractf128_ps(v, 1), pos_<0>());
	}

	LMAT_ENSURE_INLINE
	inline float avx_extract_f32(__m256 v, pos_<5> p)
	{
		return sse_extract_f32(_mm256_extractf128_ps(v, 1), pos_<1>());
	}

	LMAT_ENSURE_INLINE
	inline float avx_extract_f32(__m256 v, pos_<6> p)
	{
		return sse_extract_f32(_mm256_extractf128_ps(v, 1), pos_<2>());
	}

	LMAT_ENSURE_INLINE
	inline float avx_extract_f32(__m256 v, pos_<7> p)
	{
		return sse_extract_f32(_mm256_extractf128_ps(v, 1), pos_<3>());
	}


	LMAT_ENSURE_INLINE
	inline double avx_extract_f64(__m256d v, pos_<0> p)
	{
		return sse_extract_f64(_mm256_castpd256_pd128(v), pos_<0>());
	}

	LMAT_ENSURE_INLINE
	inline double avx_extract_f64(__m256d v, pos_<1> p)
	{
		return sse_extract_f64(_mm256_castpd256_pd128(v), pos_<1>());
	}

	LMAT_ENSURE_INLINE
	inline double avx_extract_f64(__m256d v, pos_<2> p)
	{
		return sse_extract_f64(_mm256_extractf128_ps(v, 1), pos_<0>());
	}

	LMAT_ENSURE_INLINE
	inline double avx_extract_f64(__m256d v, pos_<3> p)
	{
		return sse_extract_f64(_mm256_extractf128_ps(v, 1), pos_<1>());
	}


	// AVX-broadcasting

	LMAT_ENSURE_INLINE
	inline __m256 avx_broadcast_f32(__m256 v, pos_<0> p)
	{
		__m256 t = _mm256_shuffle_ps(v, v, 0x00);
		return _mm256_permute2f128_ps(t, t, 0x00);
	}

	LMAT_ENSURE_INLINE
	inline __m256 avx_broadcast_f32(__m256 v, pos_<1> p)
	{
		__m256 t = _mm256_shuffle_ps(v, v, 0x55);
		return _mm256_permute2f128_ps(t, t, 0x00);
	}

	LMAT_ENSURE_INLINE
	inline __m256 avx_broadcast_f32(__m256 v, pos_<2> p)
	{
		__m256 t = _mm256_shuffle_ps(v, v, 0xaa);
		return _mm256_permute2f128_ps(t, t, 0x00);
	}

	LMAT_ENSURE_INLINE
	inline __m256 avx_broadcast_f32(__m256 v, pos_<3> p)
	{
		__m256 t = _mm256_shuffle_ps(v, v, 0xff);
		return _mm256_permute2f128_ps(t, t, 0x00);
	}

	LMAT_ENSURE_INLINE
	inline __m256 avx_broadcast_f32(__m256 v, pos_<4> p)
	{
		__m256 t = _mm256_shuffle_ps(v, v, 0x00);
		return _mm256_permute2f128_ps(t, t, 0x11);
	}

	LMAT_ENSURE_INLINE
	inline __m256 avx_broadcast_f32(__m256 v, pos_<5> p)
	{
		__m256 t = _mm256_shuffle_ps(v, v, 0x55);
		return _mm256_permute2f128_ps(t, t, 0x11);
	}

	LMAT_ENSURE_INLINE
	inline __m256 avx_broadcast_f32(__m256 v, pos_<6> p)
	{
		__m256 t = _mm256_shuffle_ps(v, v, 0xaa);
		return _mm256_permute2f128_ps(t, t, 0x11);
	}

	LMAT_ENSURE_INLINE
	inline __m256 avx_broadcast_f32(__m256 v, pos_<7> p)
	{
		__m256 t = _mm256_shuffle_ps(v, v, 0xff);
		return _mm256_permute2f128_ps(t, t, 0x11);
	}


	LMAT_ENSURE_INLINE
	inline __m256d avx_broadcast_f64(__m256d v, pos_<0> p)
	{
		__m256d t = _mm256_unpacklo_pd(v, v);
		return _mm256_permute2f128_ps(t, t, 0x00);
	}

	LMAT_ENSURE_INLINE
	inline __m256d avx_broadcast_f64(__m256d v, pos_<1> p)
	{
		__m256d t = _mm256_unpackhi_pd(v, v);
		return _mm256_permute2f128_ps(t, t, 0x00);
	}

	LMAT_ENSURE_INLINE
	inline __m256d avx_broadcast_f64(__m256d v, pos_<2> p)
	{
		__m256d t = _mm256_unpacklo_pd(v, v);
		return _mm256_permute2f128_ps(t, t, 0x11);
	}

	LMAT_ENSURE_INLINE
	inline __m256d avx_broadcast_f64(__m256d v, pos_<3> p)
	{
		__m256d t = _mm256_unpackhi_pd(v, v);
		return _mm256_permute2f128_ps(t, t, 0x11);
	}



} }

#endif /* AVX_HELPERS_H_ */
