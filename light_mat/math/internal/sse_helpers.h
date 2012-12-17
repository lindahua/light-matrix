/**
 * @file sse_helpers.h
 *
 * @brief Auxiliary devices for SSE operations
 *
 * @author Dahua Lin
 */

#ifndef LIGHTMAT_SSE_HELPERS_H_
#define LIGHTMAT_SSE_HELPERS_H_

#include <light_mat/math/simd_base.h>
#include "numrepr_format.h"

namespace lmat { namespace math { namespace internal {

	// bitwise not

	LMAT_ENSURE_INLINE
	inline __m128i sse_bitwise_not(const __m128i& a)
	{
		return _mm_xor_si128(a, _mm_set1_epi32(-1));
	}


	// partial load

	LMAT_ENSURE_INLINE
	inline __m128 sse_loadpart_f32(unsigned int n, const float *p)
	{
		__m128 v;

        switch (n)
        {
        case 1:
        	v = _mm_load_ss(p);
        	break;
        case 2:
        	v = _mm_castpd_ps(_mm_load_sd((double*)p));
        	break;
        case 3:
        	v = _mm_movelh_ps(
        		_mm_castpd_ps(_mm_load_sd((double*)p)),
        		_mm_load_ss(p + 2));
        	break;
        case 4:
        	v = _mm_loadu_ps(p);
        	break;
        default:
        	v = _mm_setzero_ps();
        	break;
        }

        return v;
	}

	LMAT_ENSURE_INLINE
	inline __m128d sse_loadpart_f64(unsigned int n, const double *p)
	{
		__m128d v;

        switch (n)
        {
        case 1:
        	v = _mm_load_sd(p);
        	break;
        case 2:
        	v = _mm_loadu_pd(p);
        	break;
        default:
        	v = _mm_setzero_pd();
        	break;
        }

        return v;
	}


	// partial store

	LMAT_ENSURE_INLINE
	inline void sse_storepart_f32(unsigned int n, float *p, const __m128& v)
	{
        switch (n)
        {
        case 1:
            _mm_store_ss(p, v);
            break;
        case 2:
            _mm_store_sd((double*)p, _mm_castps_pd(v));
            break;
        case 3:
            _mm_store_sd((double*)p, _mm_castps_pd(v));
            _mm_store_ss(p + 2, _mm_movehl_ps(v, v));
            break;
        case 4:
            _mm_storeu_ps(p, v);
            break;
        }
	}


	LMAT_ENSURE_INLINE
	inline void sse_storepart_f64(unsigned int n, double *p, const __m128d& v)
	{
        switch (n)
        {
        case 1:
            _mm_store_sd(p, v);
            break;
        case 2:
            _mm_storeu_pd(p, v);
            break;
        }
	}


	// extract scalar

	LMAT_ENSURE_INLINE
	inline float sse_extract_f32(const __m128& v, unsigned int i)
	{
		float s;

    	switch (i)
    	{
    	case 0:
    		s = _mm_cvtss_f32(v);
    		break;
    	case 1:
    		s = _mm_cvtss_f32(
    				_mm_castsi128_ps(_mm_srli_si128(_mm_castps_si128(v), 4)));
    		break;
    	case 2:
    		s = _mm_cvtss_f32(
    				_mm_castsi128_ps(_mm_srli_si128(_mm_castps_si128(v), 8)));
    		break;
    	default:
    		s = _mm_cvtss_f32(
    				_mm_castsi128_ps(_mm_srli_si128(_mm_castps_si128(v), 12)));
    		break;
    	}

    	return s;
	}


	LMAT_ENSURE_INLINE
	inline double sse_extract_f64(const __m128d& v, unsigned int i)
	{
    	double s;

    	if (i == 0)
    		s = _mm_cvtsd_f64(v);
    	else
    		s = _mm_cvtsd_f64(
    				_mm_castsi128_pd(_mm_srli_si128(_mm_castpd_si128(v), 8)));

    	return s;
	}


	LMAT_ENSURE_INLINE
	inline int32_t sse_extract_i32(const __m128i& v, unsigned int i)
	{
		int32_t s;

    	switch (i)
    	{
    	case 0:
    		s = _mm_cvtsi128_si32(v);
    		break;
    	case 1:
    		s = _mm_cvtsi128_si32(_mm_srli_si128(v, 4));
    		break;
    	case 2:
    		s = _mm_cvtsi128_si32(_mm_srli_si128(v, 8));
    		break;
    	default:
    		s = _mm_cvtsi128_si32(_mm_srli_si128(v, 12));
    		break;
    	}

    	return s;
	}


	LMAT_ENSURE_INLINE
	inline int64_t sse_extract_i64(const __m128i& v, unsigned int i)
	{
		int64_t s;

    	if (i == 0)
    		s = _mm_cvtsi128_si64(v);
    	else
    		s = _mm_cvtsi128_si64(_mm_srli_si128(v, 8));

    	return s;
	}


	// classification

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



} } }


#endif /* SSE_HELPERS_H_ */
