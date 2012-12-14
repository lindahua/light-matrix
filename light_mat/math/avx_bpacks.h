/**
 * @file avx_bpacks.h
 *
 * AVX boolean packs
 * 
 * @author Dahua Lin 
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_AVX_BPACKS_H_
#define LIGHTMAT_AVX_BPACKS_H_

#include <light_mat/math/simd_base.h>
#include "internal/avx_helpers.h"

namespace lmat { namespace math {

	typedef simd_bpack<float, avx_t> avx_f32bpk;
	typedef simd_bpack<double, avx_t> avx_f64bpk;


	template<>
	struct simd_bpack<float, avx_t>
	{
		typedef int32_t bint_type;
		static const unsigned int pack_width = 8;

		union
		{
			__m256 v;
			LMAT_ALIGN_SSE int32_t e[8];
		};

		LMAT_ENSURE_INLINE
		unsigned int width() const
		{
			return pack_width;
		}

		// constructors

		LMAT_ENSURE_INLINE simd_bpack() { }

		LMAT_ENSURE_INLINE simd_bpack(const __m256& v_) : v(v_) { }

		LMAT_ENSURE_INLINE simd_bpack( bool b )
		{
			set(b);
		}

		LMAT_ENSURE_INLINE simd_bpack(
				bool b0, bool b1, bool b2, bool b3, bool b4, bool b5, bool b6, bool b7 )
		{
			set(b0, b1, b2, b3, b4, b5, b6, b7);
		}

		LMAT_ENSURE_INLINE
		static simd_bpack all_false()
		{
			return _mm256_setzero_ps();
		}

		LMAT_ENSURE_INLINE
		static simd_bpack all_true()
		{
			return _mm256_castsi256_ps(_mm256_set1_epi32(-1));
		}

		// converters

	    LMAT_ENSURE_INLINE
	    operator __m256() const
	    {
	    	return v;
	    }

	    // set values

	    LMAT_ENSURE_INLINE
	    void set(bool b)
		{
	    	if (b)
	    		v = _mm256_castsi256_ps(_mm256_set1_epi32(-1));
	    	else
	    		v = _mm256_setzero_ps();

		}

		LMAT_ENSURE_INLINE void set(bool b0, bool b1, bool b2, bool b3, bool b4, bool b5, bool b6, bool b7)
		{
			v = _mm256_castsi256_ps(_mm256_setr_epi32(
					-(int)b0, -(int)b1, -(int)b2, -(int)b3, -(int)b4, -(int)b5, -(int)b6, -(int)b7));
		}

		// extract

		LMAT_ENSURE_INLINE __m128 get_low() const
		{
			return _mm256_castps256_ps128(v);
		}

		LMAT_ENSURE_INLINE __m128 get_high() const
		{
			return _mm256_extractf128_ps(v, 1);
		}

	    LMAT_ENSURE_INLINE bool to_scalar() const
	    {
	    	return (bool)_mm_cvtsi128_si32(_mm_castps_si128(get_low()));
	    }

	    LMAT_ENSURE_INLINE bool extract(unsigned int i) const
	    {
	    	int bi;

	    	if (i < 4)
	    		bi = internal::sse_extract_i32(_mm_castps_si128(get_low()), i);
	    	else
	    		bi = internal::sse_extract_i32(_mm_castps_si128(get_high()), i & 3);

	    	return (bool)(bi);
	    }

	};



	template<>
	struct simd_bpack<double, avx_t>
	{
		typedef int32_t bint_type;
		static const unsigned int pack_width = 4;

		union
		{
			__m256d v;
			LMAT_ALIGN_SSE int64_t e[4];
		};

		LMAT_ENSURE_INLINE
		unsigned int width() const
		{
			return pack_width;
		}

		// constructors

		LMAT_ENSURE_INLINE simd_bpack() { }

		LMAT_ENSURE_INLINE simd_bpack(const __m256d& v_) : v(v_) { }

		LMAT_ENSURE_INLINE simd_bpack( bool b )
		{
			set(b);
		}

		LMAT_ENSURE_INLINE simd_bpack(bool b0, bool b1, bool b2, bool b3)
		{
			set(b0, b1, b2, b3);
		}

		LMAT_ENSURE_INLINE
		static simd_bpack all_false()
		{
			return _mm256_setzero_pd();
		}

		LMAT_ENSURE_INLINE
		static simd_bpack all_true()
		{
			return _mm256_castsi256_pd(_mm256_set1_epi32(-1));
		}

		// converters

	    LMAT_ENSURE_INLINE
	    operator __m256d() const
	    {
	    	return v;
	    }

	    // set values

	    LMAT_ENSURE_INLINE
	    void set(bool b)
		{
	    	if (b)
	    		v = _mm256_castsi256_pd(_mm256_set1_epi32(-1));
	    	else
	    		v = _mm256_setzero_pd();

		}

		LMAT_ENSURE_INLINE void set(bool b0, bool b1, bool b2, bool b3)
		{
			v = _mm256_castsi256_pd(_mm256_setr_epi32(
					-(int)b0, -(int)b0, -(int)b1, -(int)b1, -(int)b2, -(int)b2, -(int)b3, -(int)b3));
		}

		// extract

		LMAT_ENSURE_INLINE __m128d get_low() const
		{
			return _mm256_castpd256_pd128(v);
		}

		LMAT_ENSURE_INLINE __m128d get_high() const
		{
			return _mm256_extractf128_pd(v, 1);
		}

	    LMAT_ENSURE_INLINE bool to_scalar() const
	    {
	    	return (bool)_mm_cvtsi128_si64(_mm_castpd_si128(get_low()));
	    }

	    LMAT_ENSURE_INLINE bool extract(unsigned int i) const
	    {
	    	int64_t bi;

	    	if (i < 2)
	    		bi = internal::sse_extract_i64(_mm_castpd_si128(get_low()), i);
	    	else
	    		bi = internal::sse_extract_i64(_mm_castpd_si128(get_high()), i & 1);

	    	return (bool)(bi);
	    }

	};

} }


#endif 
