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
