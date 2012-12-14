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
