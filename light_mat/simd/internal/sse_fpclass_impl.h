/**
 * @file sse_fpclass_impl.h
 *
 * @brief Implementation of FP classification on SSE
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_SSE_FPCLASS_IMPL_H_
#define LIGHTMAT_SSE_FPCLASS_IMPL_H_

#include <light_mat/simd/simd_base.h>
#include "numrepr_format.h"

namespace lmat { namespace internal {

	LMAT_ENSURE_INLINE
	inline __m128i sse_expmask_ps()
	{
		typedef num_fmt<float> fmt;
		return _mm_set1_epi32(fmt::exponent_bits);
	}

	LMAT_ENSURE_INLINE
	inline __m128i sse_expmask_pd()
	{
		typedef num_fmt<double> fmt;
		return _mm_set1_epi64x(fmt::exponent_bits);
	}


	LMAT_ENSURE_INLINE
	inline __m128 sse_is_neg_ps(const __m128& a)
	{
		__m128i ai = _mm_castps_si128(a);

		__m128i sgn_m = sse_signmask_ps();
		__m128i sgn   = _mm_and_si128(sgn_m, ai);
		return _mm_castsi128_ps(_mm_cmpeq_epi32(sgn, sgn_m));
	}

	LMAT_ENSURE_INLINE
	inline __m128d sse_is_neg_pd(const __m128d& a)
	{
		__m128i ai = _mm_castpd_si128(a);

		__m128i sgn_m = sse_signmask_pd();
		__m128i sgn   = _mm_and_si128(sgn_m, ai);
		return _mm_castsi128_pd(_mm_cmpeq_epi64(sgn, sgn_m));
	}


	LMAT_ENSURE_INLINE
	inline __m128 sse_is_finite_ps(const __m128& a)
	{
		__m128i ai = _mm_castps_si128(a);

		__m128i exp_m = sse_expmask_ps();
		__m128i exp = _mm_and_si128(exp_m, ai);
		__m128i not_finite = _mm_cmpeq_epi32(exp, exp_m);
		return _mm_castsi128_ps(sse_bitwise_not(not_finite));
	}

	LMAT_ENSURE_INLINE
	inline __m128d sse_is_finite_pd(const __m128d& a)
	{
		__m128i ai = _mm_castpd_si128(a);

		__m128i exp_m = sse_expmask_pd();
		__m128i exp = _mm_and_si128(exp_m, ai);
		__m128i not_finite = _mm_cmpeq_epi64(exp, exp_m);
		return _mm_castsi128_pd(sse_bitwise_not(not_finite));
	}


	LMAT_ENSURE_INLINE
	inline __m128 sse_is_inf_ps(const __m128& a)
	{
		__m128i ai = _mm_castps_si128(a);
		ai = _mm_andnot_si128(sse_signmask_ps(), ai);

		return _mm_castsi128_ps(_mm_cmpeq_epi32(ai, sse_expmask_ps()));
	}

	LMAT_ENSURE_INLINE
	inline __m128d sse_is_inf_pd(const __m128d& a)
	{
		__m128i ai = _mm_castpd_si128(a);
		ai = _mm_andnot_si128(sse_signmask_pd(), ai);

		return _mm_castsi128_pd(_mm_cmpeq_epi64(ai, sse_expmask_pd()));
	}


	LMAT_ENSURE_INLINE
	inline __m128 sse_is_nan_ps(const __m128& a)
	{
		typedef num_fmt<float> fmt;

		__m128i ai = _mm_castps_si128(a);
		__m128i e_msk = _mm_set1_epi32(fmt::exponent_bits);
		__m128i s_msk = _mm_set1_epi32(fmt::mantissa_bits);

		__m128i e = _mm_and_si128(ai, e_msk);
		__m128i s = _mm_and_si128(ai, s_msk);

		return _mm_castsi128_ps(_mm_andnot_si128(
				_mm_cmpeq_epi32(s, _mm_setzero_si128()),
				_mm_cmpeq_epi32(e, e_msk)));
	}

	LMAT_ENSURE_INLINE
	inline __m128d sse_is_nan_pd(const __m128d& a)
	{
		typedef num_fmt<double> fmt;

		__m128i ai = _mm_castpd_si128(a);
		__m128i e_msk = _mm_set1_epi64x(fmt::exponent_bits);
		__m128i s_msk = _mm_set1_epi64x(fmt::mantissa_bits);

		__m128i e = _mm_and_si128(ai, e_msk);
		__m128i s = _mm_and_si128(ai, s_msk);

		return _mm_castsi128_pd(_mm_andnot_si128(
				_mm_cmpeq_epi64(s, _mm_setzero_si128()),
				_mm_cmpeq_epi64(e, e_msk)));
	}

} }

#endif
