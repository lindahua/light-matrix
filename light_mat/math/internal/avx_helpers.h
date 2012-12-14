/**
 * @file avx_helpers.h
 *
 * @brief
 *
 * @author Dahua Lin
 */

#ifndef LIGHTMAT_AVX_HELPERS_H_
#define LIGHTMAT_AVX_HELPERS_H_

#include "sse_helpers.h"

namespace lmat { namespace math { namespace internal {


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


	LMAT_ENSURE_INLINE
	inline __m256 mm256_bitwise_eq_ps(const __m256& a, const __m256& b)
	{
		__m128i ai_lo = _mm_castps_si128(_mm256_castps256_ps128(a));
		__m128i ai_hi = _mm_castps_si128(_mm256_extractf128_ps(a, 1));

		__m128i bi_lo = _mm_castps_si128(_mm256_castps256_ps128(b));
		__m128i bi_hi = _mm_castps_si128(_mm256_extractf128_ps(b, 1));

		return combine_m128(
				_mm_castsi128_ps(_mm_cmpeq_epi32(ai_lo, bi_lo)),
				_mm_castsi128_ps(_mm_cmpeq_epi32(ai_hi, bi_hi)));
	}






	LMAT_ENSURE_INLINE
	inline __m256i avx_part_mask_32(unsigned int n)
	{
		__m256i m;

		switch (n)
		{
		case 1:
			m = _mm256_setr_epi32(-1, 0, 0, 0, 0, 0, 0, 0);
			break;
		case 2:
			m = _mm256_setr_epi32(-1, -1, 0, 0, 0, 0, 0, 0);
			break;
		case 3:
			m = _mm256_setr_epi32(-1, -1, -1, 0, 0, 0, 0, 0);
			break;
		case 4:
			m = _mm256_setr_epi32(-1, -1, -1, -1, 0, 0, 0, 0);
			break;
		case 5:
			m = _mm256_setr_epi32(-1, -1, -1, -1, -1, 0, 0, 0);
			break;
		case 6:
			m = _mm256_setr_epi32(-1, -1, -1, -1, -1, -1, 0, 0);
			break;
		case 7:
			m = _mm256_setr_epi32(-1, -1, -1, -1, -1, -1, -1, 0);
			break;
		case 8:
			m = _mm256_set1_epi32(-1);
			break;
		default:
			m = _mm256_setzero_si256();
			break;
		}

		return m;
	}


	LMAT_ENSURE_INLINE
	inline __m256i avx_part_mask_64(unsigned int n)
	{
		__m256i m;

		switch (n)
		{
		case 1:
			m = _mm256_setr_epi32(-1, -1, 0, 0, 0, 0, 0, 0);
			break;
		case 2:
			m = _mm256_setr_epi32(-1, -1, -1, -1, 0, 0, 0, 0);
			break;
		case 3:
			m = _mm256_setr_epi32(-1, -1, -1, -1, -1, -1, 0, 0);
			break;
		case 4:
			m = _mm256_set1_epi32(-1);
			break;
		default:
			m = _mm256_setzero_si256();
			break;
		}

		return m;
	}


} } }

#endif /* AVX_HELPERS_H_ */
