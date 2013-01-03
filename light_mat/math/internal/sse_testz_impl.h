/**
 * @file sse_testz_impl.h
 *
 * @brief SSE implementation of testing whether all bits are zeros/ones
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_SSE_TESTZ_IMPL_H_
#define LIGHTMAT_SSE_TESTZ_IMPL_H_

#include <light_mat/math/simd_base.h>

namespace lmat { namespace math { namespace internal {

	LMAT_ENSURE_INLINE
	inline bool testz_sse2(const __m128i& a) // returns true when all bits in a are 0s
	{
	    __m128i t = _mm_or_si128(a, _mm_unpackhi_epi64(a, a));
	    return _mm_cvtsi128_si64(t) == 0;
	}

	LMAT_ENSURE_INLINE
	inline bool testc_sse2(const __m128i& a) // returns true when all bits in a set 1s
	{
	    __m128i t = _mm_and_si128(a, _mm_unpackhi_epi64(a, a));
	    return _mm_cvtsi128_si64(t) == int64_t(-1);
	}

#if (defined(LMAT_HAS_SSE4_1))

	LMAT_ENSURE_INLINE
	inline bool testz_sse4(const __m128i& a)
	{
		return (bool)_mm_testz_si128(a, _mm_setr_epi32(-1, -1, -1, -1));
	}

	LMAT_ENSURE_INLINE
	inline bool testc_sse4(const __m128i& a)
	{
		return (bool)_mm_testc_si128(a, _mm_setr_epi32(-1, -1, -1, -1));
	}


#endif

	LMAT_ENSURE_INLINE
	inline bool testz(const __m128i& a)
	{
#if (defined(LMAT_HAS_SSE4_1))
		return testz_sse4(a);
#else
		return testz_sse2(a);
#endif
	}

	LMAT_ENSURE_INLINE
	inline bool testc(const __m128i& a)
	{
#if (defined(LMAT_HAS_SSE4_1))
		return testc_sse4(a);
#else
		return testc_sse2(a);
#endif
	}



} } }

#endif /* SSE2_ALLANY_IMPL_H_ */
