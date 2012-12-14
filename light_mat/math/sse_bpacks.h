/**
 * @file sse_bpacks.h
 *
 * @brief SSE boolean pack classes
 *
 * @author Dahua Lin
 */

#ifndef LIGHTMAT_SSE_BPACKS_H_
#define LIGHTMAT_SSE_BPACKS_H_

#include <light_mat/math/simd_base.h>

namespace lmat { namespace math {

	typedef simd_bpack<float, sse_t> sse_f32bpk;
	typedef simd_bpack<double, sse_t> sse_f64bpk;


	template<>
	struct simd_bpack<float, sse_t>
	{
		typedef int32_t bint_type;
		static const unsigned int pack_width = 4;

		union
		{
			__m128 v;
			LMAT_ALIGN_SSE int32_t e[4];
		};

		LMAT_ENSURE_INLINE
		unsigned int width() const
		{
			return pack_width;
		}

		// constructors

		LMAT_ENSURE_INLINE simd_bpack() { }

		LMAT_ENSURE_INLINE simd_bpack(const __m128& v_) : v(v_) { }

		LMAT_ENSURE_INLINE simd_bpack( bool b )
		{
			set(b);
		}

		LMAT_ENSURE_INLINE simd_bpack( bool b0, bool b1, bool b2, bool b3 )
		{
			set(b0, b1, b2, b3);
		}

		LMAT_ENSURE_INLINE
		static simd_bpack all_false()
		{
			return _mm_setzero_ps();
		}

		LMAT_ENSURE_INLINE
		static simd_bpack all_true()
		{
			return _mm_castsi128_ps(_mm_set1_epi32(-1));
		}

		// converters

	    LMAT_ENSURE_INLINE
	    operator __m128() const
	    {
	    	return v;
	    }

	    // set values

	    LMAT_ENSURE_INLINE
	    void set(bool b)
		{
	    	if (b)
	    		v = _mm_castsi128_ps(_mm_set1_epi32(-1));
	    	else
	    		v = _mm_setzero_ps();

		}

		LMAT_ENSURE_INLINE void set(bool b0, bool b1, bool b2, bool b3)
		{
			v = _mm_castsi128_ps(_mm_setr_epi32(-(int)b0, -(int)b1, -(int)b2, -(int)b3));
		}

		// extract

	    LMAT_ENSURE_INLINE bool to_scalar() const
	    {
	    	return (bool)_mm_extract_epi32(_mm_castps_si128(v), 0);
	    }

	    LMAT_ENSURE_INLINE bool extract(int i) const
	    {
	    	return (bool)_mm_extract_epi32(_mm_castps_si128(v), i);
	    }
	};


	template<>
	struct simd_bpack<double, sse_t>
	{
		typedef int64_t bint_type;
		static const unsigned int pack_width = 2;

		union
		{
			__m128d v;
			LMAT_ALIGN_SSE int64_t e[2];
		};

		LMAT_ENSURE_INLINE
		unsigned int width() const
		{
			return pack_width;
		}

		// constructors

		LMAT_ENSURE_INLINE simd_bpack() { }

		LMAT_ENSURE_INLINE simd_bpack(const __m128d& v_) : v(v_) { }

		LMAT_ENSURE_INLINE simd_bpack( bool b )
		{
			set(b);
		}

		LMAT_ENSURE_INLINE simd_bpack( bool b0, bool b1 )
		{
			set(b0, b1);
		}

		LMAT_ENSURE_INLINE
		static simd_bpack all_false()
		{
			return _mm_setzero_pd();
		}

		LMAT_ENSURE_INLINE
		static simd_bpack all_true()
		{
			return _mm_castsi128_pd(_mm_set1_epi32(-1));
		}

		// converters

	    LMAT_ENSURE_INLINE
	    operator __m128d() const
	    {
	    	return v;
	    }

	    // set values

	    LMAT_ENSURE_INLINE
	    void set(bool b)
		{
	    	if (b)
	    		v = _mm_castsi128_pd(_mm_set1_epi32(-1));
	    	else
	    		v = _mm_setzero_pd();
		}

		LMAT_ENSURE_INLINE void set(bool b0, bool b1)
		{
			v = _mm_castsi128_pd(_mm_setr_epi32(-(int)b0, -(int)b0, -(int)b1, -(int)b1));
		}

		// extract

	    LMAT_ENSURE_INLINE bool to_scalar() const
	    {
	    	return (bool)_mm_extract_epi64(_mm_castpd_si128(v), 0);
	    }

	    LMAT_ENSURE_INLINE bool extract(int i) const
	    {
	    	return (bool)_mm_extract_epi64(_mm_castpd_si128(v), i);
	    }
	};


} }

#endif /* SSE_BPACKS_H_ */
