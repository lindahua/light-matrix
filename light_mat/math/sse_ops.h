/**
 * @file sse_ops.h
 *
 * @brief Basic operations on SSE packs
 *
 * @author Dahua Lin
 */

#ifndef LIGHTMAT_SSE_OPS_H_
#define LIGHTMAT_SSE_OPS_H_

#include <light_mat/math/sse_packs.h>
#include <light_mat/math/sse_bpacks.h>
#include "internal/numrepr_format.h"
#include "internal/sse2_round_impl.h"

namespace lmat { namespace math {

	/********************************************
	 *
	 *  Floating-point arithmetics
	 *
	 ********************************************/

	LMAT_ENSURE_INLINE
	inline sse_f32pk operator + (const sse_f32pk& a, const sse_f32pk& b)
	{
		return _mm_add_ps(a.v, b.v);
	}

	LMAT_ENSURE_INLINE
	inline sse_f64pk operator + (const sse_f64pk& a, const sse_f64pk& b)
	{
		return _mm_add_pd(a.v, b.v);
	}

	LMAT_ENSURE_INLINE
	inline sse_f32pk operator - (const sse_f32pk& a, const sse_f32pk& b)
	{
		return _mm_sub_ps(a.v, b.v);
	}

	LMAT_ENSURE_INLINE
	inline sse_f64pk operator - (const sse_f64pk& a, const sse_f64pk& b)
	{
		return _mm_sub_pd(a.v, b.v);
	}

	LMAT_ENSURE_INLINE
	inline sse_f32pk operator * (const sse_f32pk& a, const sse_f32pk& b)
	{
		return _mm_mul_ps(a.v, b.v);
	}

	LMAT_ENSURE_INLINE
	inline sse_f64pk operator * (const sse_f64pk& a, const sse_f64pk& b)
	{
		return _mm_mul_pd(a.v, b.v);
	}

	LMAT_ENSURE_INLINE
	inline sse_f32pk operator / (const sse_f32pk& a, const sse_f32pk& b)
	{
		return _mm_div_ps(a.v, b.v);
	}

	LMAT_ENSURE_INLINE
	inline sse_f64pk operator / (const sse_f64pk& a, const sse_f64pk& b)
	{
		return _mm_div_pd(a.v, b.v);
	}

	LMAT_ENSURE_INLINE
	inline sse_f32pk operator - (const sse_f32pk& a)
	{
		typedef internal::num_fmt<float> fmt;
		return _mm_xor_ps(_mm_castsi128_ps(_mm_set1_epi32(fmt::sign_bit)), a.v);
	}

	LMAT_ENSURE_INLINE
	inline sse_f64pk operator - (const sse_f64pk& a)
	{
		typedef internal::num_fmt<double> fmt;
		return _mm_xor_pd(_mm_castsi128_pd(_mm_set1_epi64x(fmt::sign_bit)), a.v);
	}

	LMAT_ENSURE_INLINE
	inline sse_f32pk abs(const sse_f32pk& a)
	{
		typedef internal::num_fmt<float> fmt;
		return _mm_andnot_ps(_mm_castsi128_ps(_mm_set1_epi32(fmt::sign_bit)), a.v);
	}

	LMAT_ENSURE_INLINE
	inline sse_f64pk abs(const sse_f64pk& a)
	{
		typedef internal::num_fmt<double> fmt;
		return _mm_andnot_pd(_mm_castsi128_pd(_mm_set1_epi64x(fmt::sign_bit)), a.v);
	}

	LMAT_ENSURE_INLINE
	inline sse_f32pk& operator += (sse_f32pk& a, const sse_f32pk& b)
	{
		a.v = _mm_add_ps(a.v, b.v);
		return a;
	}

	LMAT_ENSURE_INLINE
	inline sse_f64pk& operator += (sse_f64pk& a, const sse_f64pk& b)
	{
		a.v = _mm_add_pd(a.v, b.v);
		return a;
	}

	LMAT_ENSURE_INLINE
	inline sse_f32pk& operator -= (sse_f32pk& a, const sse_f32pk& b)
	{
		a.v = _mm_sub_ps(a.v, b.v);
		return a;
	}

	LMAT_ENSURE_INLINE
	inline sse_f64pk& operator -= (sse_f64pk& a, const sse_f64pk& b)
	{
		a.v = _mm_sub_pd(a.v, b.v);
		return a;
	}

	LMAT_ENSURE_INLINE
	inline sse_f32pk& operator *= (sse_f32pk& a, const sse_f32pk& b)
	{
		a.v = _mm_mul_ps(a.v, b.v);
		return a;
	}

	LMAT_ENSURE_INLINE
	inline sse_f64pk& operator *= (sse_f64pk& a, const sse_f64pk& b)
	{
		a.v = _mm_mul_pd(a.v, b.v);
		return a;
	}

	LMAT_ENSURE_INLINE
	inline sse_f32pk& operator /= (sse_f32pk& a, const sse_f32pk& b)
	{
		a.v = _mm_div_ps(a.v, b.v);
		return a;
	}

	LMAT_ENSURE_INLINE
	inline sse_f64pk& operator /= (sse_f64pk& a, const sse_f64pk& b)
	{
		a.v = _mm_div_pd(a.v, b.v);
		return a;
	}


	/********************************************
	 *
	 *  Floating-point min and max
	 *
	 ********************************************/

	LMAT_ENSURE_INLINE
	inline sse_f32pk (min)(const sse_f32pk& a, const sse_f32pk& b)
	{
		return _mm_min_ps(a.v, b.v);
	}

	LMAT_ENSURE_INLINE
	inline sse_f64pk (min)(const sse_f64pk& a, const sse_f64pk& b)
	{
		return _mm_min_pd(a.v, b.v);
	}

	LMAT_ENSURE_INLINE
	inline sse_f32pk (max)(const sse_f32pk& a, const sse_f32pk& b)
	{
		return _mm_max_ps(a.v, b.v);
	}
	LMAT_ENSURE_INLINE
	inline sse_f64pk (max)(const sse_f64pk& a, const sse_f64pk& b)
	{
		return _mm_max_pd(a.v, b.v);
	}

	/********************************************
	 *
	 *  Simple power functions
	 *
	 ********************************************/

	LMAT_ENSURE_INLINE
	inline sse_f32pk sqr(const sse_f32pk& a)
	{
		return _mm_mul_ps(a.v, a.v);
	}

	LMAT_ENSURE_INLINE
	inline sse_f64pk sqr(const sse_f64pk& a)
	{
		return _mm_mul_pd(a.v, a.v);
	}

	LMAT_ENSURE_INLINE
	inline sse_f32pk cube(const sse_f32pk& a)
	{
		return _mm_mul_ps(_mm_mul_ps(a.v, a.v), a.v);
	}

	LMAT_ENSURE_INLINE
	inline sse_f64pk cube(const sse_f64pk& a)
	{
		return _mm_mul_pd(_mm_mul_pd(a.v, a.v), a.v);
	}

	LMAT_ENSURE_INLINE
	inline sse_f32pk sqrt(const sse_f32pk& a)
	{
		return _mm_sqrt_ps(a.v);
	}

	LMAT_ENSURE_INLINE
	inline sse_f64pk sqrt(const sse_f64pk& a)
	{
		return _mm_sqrt_pd(a.v);
	}

	LMAT_ENSURE_INLINE
	inline sse_f32pk rcp(const sse_f32pk& a)
	{
		return _mm_div_ps(_mm_set1_ps(1.0f), a.v);
	}

	LMAT_ENSURE_INLINE
	inline sse_f32pk approx_rcp(const sse_f32pk& a)
	{
		return _mm_rcp_ps(a.v);
	}

	LMAT_ENSURE_INLINE
	inline sse_f64pk rcp(const sse_f64pk& a)
	{
		return _mm_div_pd(_mm_set1_pd(1.0), a.v);
	}

	LMAT_ENSURE_INLINE
	inline sse_f32pk rsqrt(const sse_f32pk& a)
	{
		return _mm_div_ps(_mm_set1_ps(1.0f), _mm_sqrt_ps(a.v));
	}

	LMAT_ENSURE_INLINE
	inline sse_f32pk approx_rsqrt(const sse_f32pk& a)
	{
		return _mm_rsqrt_ps(a.v);
	}

	LMAT_ENSURE_INLINE
	inline sse_f64pk rsqrt(const sse_f64pk& a)
	{
		return _mm_div_pd(_mm_set1_pd(1.0), _mm_sqrt_pd(a.v));
	}


	/********************************************
	 *
	 *  comparison operator
	 *
	 ********************************************/

	LMAT_ENSURE_INLINE
	inline sse_f32bpk operator == (const sse_f32pk& a, const sse_f32pk& b)
	{
		return _mm_cmpeq_ps(a, b);
	}

	LMAT_ENSURE_INLINE
	inline sse_f64bpk operator == (const sse_f64pk& a, const sse_f64pk& b)
	{
		return _mm_cmpeq_pd(a, b);
	}

	LMAT_ENSURE_INLINE
	inline sse_f32bpk operator != (const sse_f32pk& a, const sse_f32pk& b)
	{
		return _mm_cmpneq_ps(a, b);
	}

	LMAT_ENSURE_INLINE
	inline sse_f64bpk operator != (const sse_f64pk& a, const sse_f64pk& b)
	{
		return _mm_cmpneq_pd(a, b);
	}

	LMAT_ENSURE_INLINE
	inline sse_f32bpk operator > (const sse_f32pk& a, const sse_f32pk& b)
	{
		return _mm_cmpgt_ps(a, b);
	}

	LMAT_ENSURE_INLINE
	inline sse_f64bpk operator > (const sse_f64pk& a, const sse_f64pk& b)
	{
		return _mm_cmpgt_pd(a, b);
	}

	LMAT_ENSURE_INLINE
	inline sse_f32bpk operator >= (const sse_f32pk& a, const sse_f32pk& b)
	{
		return _mm_cmpge_ps(a, b);
	}

	LMAT_ENSURE_INLINE
	inline sse_f64bpk operator >= (const sse_f64pk& a, const sse_f64pk& b)
	{
		return _mm_cmpge_pd(a, b);
	}

	LMAT_ENSURE_INLINE
	inline sse_f32bpk operator < (const sse_f32pk& a, const sse_f32pk& b)
	{
		return _mm_cmplt_ps(a, b);
	}

	LMAT_ENSURE_INLINE
	inline sse_f64bpk operator < (const sse_f64pk& a, const sse_f64pk& b)
	{
		return _mm_cmplt_pd(a, b);
	}


	LMAT_ENSURE_INLINE
	inline sse_f32bpk operator <= (const sse_f32pk& a, const sse_f32pk& b)
	{
		return _mm_cmple_ps(a, b);
	}

	LMAT_ENSURE_INLINE
	inline sse_f64bpk operator <= (const sse_f64pk& a, const sse_f64pk& b)
	{
		return _mm_cmple_pd(a, b);
	}


	/********************************************
	 *
	 *  logical operations
	 *
	 ********************************************/

	LMAT_ENSURE_INLINE
	inline sse_f32bpk operator ~ (const sse_f32bpk& a)
	{
		return _mm_cmpeq_epi32(_mm_castps_si128(a.v), _mm_setzero_si128());
	}

	LMAT_ENSURE_INLINE
	inline sse_f64bpk operator ~ (const sse_f64bpk& a)
	{
		return _mm_cmpeq_epi64(_mm_castpd_si128(a.v), _mm_setzero_si128());
	}

	LMAT_ENSURE_INLINE
	inline sse_f32bpk operator & (const sse_f32bpk& a, const sse_f32bpk& b)
	{
		return _mm_and_ps(a.v, b.v);
	}

	LMAT_ENSURE_INLINE
	inline sse_f64bpk operator & (const sse_f64bpk& a, const sse_f64bpk& b)
	{
		return _mm_and_pd(a.v, b.v);
	}

	LMAT_ENSURE_INLINE
	inline sse_f32bpk operator | (const sse_f32bpk& a, const sse_f32bpk& b)
	{
		return _mm_or_ps(a.v, b.v);
	}

	LMAT_ENSURE_INLINE
	inline sse_f64bpk operator | (const sse_f64bpk& a, const sse_f64bpk& b)
	{
		return _mm_or_pd(a.v, b.v);
	}

	LMAT_ENSURE_INLINE
	inline sse_f32bpk operator == (const sse_f32bpk& a, const sse_f32bpk& b)
	{
		return _mm_cmpeq_epi32(_mm_castps_si128(a.v), _mm_castps_si128(b.v));
	}

	LMAT_ENSURE_INLINE
	inline sse_f64bpk operator == (const sse_f64bpk& a, const sse_f64bpk& b)
	{
		return _mm_cmpeq_epi64(_mm_castpd_si128(a.v), _mm_castpd_si128(b.v));
	}

	LMAT_ENSURE_INLINE
	inline sse_f32bpk operator != (const sse_f32bpk& a, const sse_f32bpk& b)
	{
		return ~(a == b);
	}

	LMAT_ENSURE_INLINE
	inline sse_f64bpk operator != (const sse_f64bpk& a, const sse_f64bpk& b)
	{
		return ~(a == b);
	}


	LMAT_ENSURE_INLINE
	sse_f32bpk& operator &= (sse_f32bpk& a, const sse_f32bpk& b)
	{
		a = a & b;
		return a;
	}

	LMAT_ENSURE_INLINE
	sse_f64bpk& operator &= (sse_f64bpk& a, const sse_f64bpk& b)
	{
		a = a & b;
		return a;
	}

	LMAT_ENSURE_INLINE
	inline sse_f32bpk& operator |= (sse_f32bpk& a, const sse_f32bpk& b)
	{
		a = a | b;
		return a;
	}

	LMAT_ENSURE_INLINE
	inline sse_f64bpk& operator |= (sse_f64bpk& a, const sse_f64bpk& b)
	{
		a = a | b;
		return a;
	}


	/********************************************
	 *
	 *  rounding
	 *
	 ********************************************/

	LMAT_ENSURE_INLINE
	inline sse_f32pk round(const sse_f32pk& a)
	{
#ifdef LMAT_HAS_SSE4_1
		return _mm_round_ps(a.v, 0);
#else
		return internal::round_sse2(a);
#endif
	}

	LMAT_ENSURE_INLINE
	inline sse_f64pk round(const sse_f64pk& a)
	{
#ifdef LMAT_HAS_SSE4_1
		return _mm_round_pd(a.v, 0);
#else
		return internal::round_sse2(a);
#endif
	}

	LMAT_ENSURE_INLINE
	inline sse_f32pk floor(const sse_f32pk& a)
	{
#ifdef LMAT_HAS_SSE4_1
		return _mm_round_ps(a.v, 1);
#else
		return internal::floor_sse2(a);
#endif
	}

	LMAT_ENSURE_INLINE
	inline sse_f64pk floor(const sse_f64pk& a)
	{
#ifdef LMAT_HAS_SSE4_1
		return _mm_round_pd(a.v, 1);
#else
		return internal::floor_sse2(a);
#endif
	}

	LMAT_ENSURE_INLINE
	inline sse_f32pk ceil(const sse_f32pk& a)
	{
#ifdef LMAT_HAS_SSE4_1
		return _mm_round_ps(a.v, 2);
#else
		return internal::ceil_sse2(a);
#endif
	}

	LMAT_ENSURE_INLINE
	inline sse_f64pk ceil(const sse_f64pk& a)
	{
#ifdef LMAT_HAS_SSE4_1
		return _mm_round_pd(a.v, 2);
#else
		return internal::ceil_sse2(a);
#endif
	}

	LMAT_ENSURE_INLINE
	inline sse_f32pk trunc(const sse_f32pk& a)
	{
#ifdef LMAT_HAS_SSE4_1
		return _mm_round_ps(a.v, 3);
#else
		return internal::trunc_sse2(a);
#endif
	}

	LMAT_ENSURE_INLINE
	inline sse_f64pk trunc(const sse_f64pk& a)
	{
#ifdef LMAT_HAS_SSE4_1
		return _mm_round_pd(a.v, 3);
#else
		return internal::trunc_sse2(a);
#endif
	}



	/********************************************
	 *
	 *  FP classification
	 *
	 ********************************************/

	LMAT_ENSURE_INLINE
	inline sse_f32bpk is_neg(const sse_f32pk& a)
	{
		typedef internal::num_fmt<float> fmt;

		__m128i msk = _mm_set1_epi32(fmt::sign_bit);
		__m128i u = _mm_and_si128(_mm_castps_si128(a.v), msk);
		return _mm_castsi128_ps(_mm_cmpeq_epi32(u, msk));
	}

	LMAT_ENSURE_INLINE
	inline sse_f64bpk is_neg(const sse_f64pk& a)
	{
		typedef internal::num_fmt<double> fmt;

		__m128i msk = _mm_set1_epi64x(fmt::sign_bit);
		__m128i u = _mm_and_si128(_mm_castpd_si128(a.v), msk);
		return _mm_castsi128_pd(_mm_cmpeq_epi64(u, msk));
	}


	LMAT_ENSURE_INLINE
	inline sse_f32bpk is_finite(const sse_f32pk& a)
	{
		typedef internal::num_fmt<float> fmt;

		__m128i e_msk = _mm_set1_epi32(fmt::exponent_bits);
		__m128i e = _mm_and_si128(_mm_castps_si128(a.v), e_msk);
		__m128i not_finite = _mm_castsi128_ps(_mm_cmpeq_epi32(e, e_msk));
		return _mm_castsi128_ps(
				_mm_cmpeq_epi32(not_finite, _mm_setzero_si128()));
	}

	LMAT_ENSURE_INLINE
	inline sse_f64bpk is_finite(const sse_f64pk& a)
	{
		typedef internal::num_fmt<double> fmt;

		__m128i e_msk = _mm_set1_epi64x(fmt::exponent_bits);
		__m128i e = _mm_and_si128(_mm_castpd_si128(a.v), e_msk);
		__m128i not_finite = _mm_castsi128_pd(_mm_cmpeq_epi64(e, e_msk));
		return _mm_castsi128_pd(
				_mm_cmpeq_epi64(not_finite, _mm_setzero_si128()));
	}


	LMAT_ENSURE_INLINE
	inline sse_f32bpk is_inf(const sse_f32pk& a)
	{
		typedef internal::num_fmt<float> fmt;

		__m128i aa = _mm_castps_si128(abs(a).v);
		__m128i val = _mm_set1_epi32(fmt::exponent_bits);
		return _mm_castsi128_ps(_mm_cmpeq_epi32(aa, val));
	}

	LMAT_ENSURE_INLINE
	inline sse_f64bpk is_inf(const sse_f64pk& a)
	{
		typedef internal::num_fmt<double> fmt;

		__m128i aa = _mm_castpd_si128(abs(a).v);
		__m128i val = _mm_set1_epi64x(fmt::exponent_bits);
		return _mm_castsi128_pd(_mm_cmpeq_epi64(aa, val));
	}


	LMAT_ENSURE_INLINE
	inline sse_f32bpk is_nan(const sse_f32pk& a)
	{
		typedef internal::num_fmt<float> fmt;

		__m128i ai = _mm_castps_si128(a.v);
		__m128i e_msk = _mm_set1_epi32(fmt::exponent_bits);
		__m128i s_msk = _mm_set1_epi32(fmt::mantissa_bits);

		__m128i e = _mm_and_si128(ai, e_msk);
		__m128i s = _mm_and_si128(ai, s_msk);

		return _mm_castsi128_ps(_mm_andnot_si128(
				_mm_cmpeq_epi32(s, _mm_setzero_si128()),
				_mm_cmpeq_epi32(e, e_msk)));
	}

	LMAT_ENSURE_INLINE
	inline sse_f64bpk is_nan(const sse_f64pk& a)
	{
		typedef internal::num_fmt<double> fmt;

		__m128i ai = _mm_castpd_si128(a.v);
		__m128i e_msk = _mm_set1_epi64x(fmt::exponent_bits);
		__m128i s_msk = _mm_set1_epi64x(fmt::mantissa_bits);

		__m128i e = _mm_and_si128(ai, e_msk);
		__m128i s = _mm_and_si128(ai, s_msk);

		return _mm_castsi128_pd(_mm_andnot_si128(
				_mm_cmpeq_epi64(s, _mm_setzero_si128()),
				_mm_cmpeq_epi64(e, e_msk)));
	}




} }

#endif /* SSE_ARITH_H_ */
