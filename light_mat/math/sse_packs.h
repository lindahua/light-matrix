/**
 * @file sse_packs.h
 *
 * @brief SSE pack classes
 *
 * @author Dahua Lin
 */

#ifndef LIGHTMAT_SSE_PACKS_H_
#define LIGHTMAT_SSE_PACKS_H_

#include <light_mat/math/simd_base.h>
#include "internal/sse_helpers.h"

namespace lmat { namespace math {


	/********************************************
	 *
	 *  trait classes
	 *
	 ********************************************/

	LMAT_DEFINE_SIMD_TRAITS( sse_t, float,  4, 16 )
	LMAT_DEFINE_SIMD_TRAITS( sse_t, double, 2, 16 )


	/********************************************
	 *
	 *  pack classes
	 *
	 ********************************************/

	typedef simd_pack<float,  sse_t> sse_f32pk;
	typedef simd_pack<double, sse_t> sse_f64pk;


	template<>
	struct simd_pack<float, sse_t>
	{
		LMAT_DEFINE_FOR_SIMD_PACK( sse_t, float, 4 )

		union
		{
			__m128 v;
			LMAT_ALIGN_SSE float e[4];
		};

		LMAT_ENSURE_INLINE
		unsigned int width() const
		{
			return pack_width;
		}


		// constructors

		LMAT_ENSURE_INLINE simd_pack() { }

		LMAT_ENSURE_INLINE simd_pack(const __m128& v_) : v(v_) { }

		LMAT_ENSURE_INLINE simd_pack(const float& ev)
		{
			v = _mm_set1_ps(ev);
		}

		LMAT_ENSURE_INLINE simd_pack(
				const float& e0, const float& e1, const float& e2, const float& e3)
		{
			v = _mm_setr_ps(e0, e1, e2, e3);
		}

	    LMAT_ENSURE_INLINE
	    static simd_pack zeros()
	    {
	    	return _mm_setzero_ps();
	    }

	    LMAT_ENSURE_INLINE
	    static simd_pack ones()
	    {
	    	return _mm_set1_ps(1.0f);
	    }

	    LMAT_ENSURE_INLINE
	    static simd_pack inf()
	    {
	    	return _mm_castsi128_ps(_mm_set1_epi32((int)0x7f800000));
	    }

	    LMAT_ENSURE_INLINE
	    static simd_pack neg_inf()
	    {
	    	return _mm_castsi128_ps(_mm_set1_epi32((int)0xff800000));
	    }

	    LMAT_ENSURE_INLINE
	    static simd_pack nan()
	    {
	    	return _mm_set1_ps(std::numeric_limits<float>::quiet_NaN());
	    }


	    // converter

	    LMAT_ENSURE_INLINE
	    operator __m128() const
	    {
	    	return v;
	    }


		// set

		LMAT_ENSURE_INLINE void reset()
		{
			v = _mm_setzero_ps();
		}

		LMAT_ENSURE_INLINE void set(const float& ev)
		{
			v = _mm_set1_ps(ev);
		}

		LMAT_ENSURE_INLINE void set(
				const float& e0, const float& e1, const float& e2, const float& e3)
		{
			v = _mm_setr_ps(e0, e1, e2, e3);
		}


		// load

		LMAT_ENSURE_INLINE void load_u(const float *p)
		{
			v = _mm_loadu_ps(p);
		}

		LMAT_ENSURE_INLINE void load_a(const float *p)
		{
			v = _mm_load_ps(p);
		}

	    LMAT_ENSURE_INLINE void load_part(unsigned int n, float const * p)
	    {
	    	v = internal::sse_loadpart_f32(n, p);
	    }

	    // store

	    LMAT_ENSURE_INLINE void store_u(float *p) const
	    {
	    	_mm_storeu_ps(p, v);
	    }

	    LMAT_ENSURE_INLINE void store_a(float *p) const
	    {
	    	_mm_store_ps(p, v);
	    }

	    LMAT_ENSURE_INLINE void store_part(unsigned int n, float *p) const
	    {
	    	internal::sse_storepart_f32(n, p, v);
	    }


	    // extract

	    LMAT_ENSURE_INLINE float to_scalar() const
	    {
	    	return _mm_cvtss_f32(v);
	    }

	    LMAT_ENSURE_INLINE float extract(unsigned int i) const
	    {
	    	return internal::sse_extract_f32(v, i);
	    }


	}; // SSE f32 pack



	template<>
	struct simd_pack<double, sse_t>
	{
		LMAT_DEFINE_FOR_SIMD_PACK( sse_t, double, 2 )

		union
		{
			__m128d v;
			LMAT_ALIGN_SSE double e[2];
		};

		LMAT_ENSURE_INLINE
		unsigned int width() const
		{
			return pack_width;
		}

		// constructors

		LMAT_ENSURE_INLINE simd_pack() { }

		LMAT_ENSURE_INLINE simd_pack(const __m128d& v_) : v(v_) { }

		LMAT_ENSURE_INLINE simd_pack(const double& ev)
		{
			v = _mm_set1_pd(ev);
		}

		LMAT_ENSURE_INLINE simd_pack(const double& e0, const double& e1)
		{
			v = _mm_setr_pd(e0, e1);
		}

	    LMAT_ENSURE_INLINE
	    static simd_pack zeros()
	    {
	    	return _mm_setzero_pd();
	    }

	    LMAT_ENSURE_INLINE
	    static simd_pack ones()
	    {
	    	return _mm_set1_pd(1.0);
	    }

	    LMAT_ENSURE_INLINE
	    static simd_pack inf()
	    {
	    	return _mm_castsi128_pd(_mm_setr_epi32(0, (int)0x7ff00000, 0, (int)0x7ff00000));
	    }

	    LMAT_ENSURE_INLINE
	    static simd_pack neg_inf()
	    {
	    	return _mm_castsi128_pd(_mm_setr_epi32(0, (int)0xfff00000, 0, (int)0xfff00000));
	    }

	    LMAT_ENSURE_INLINE
	    static simd_pack nan()
	    {
	    	return _mm_set1_pd(std::numeric_limits<double>::quiet_NaN());
	    }

	    // converter

	    LMAT_ENSURE_INLINE
	    operator __m128d() const
	    {
	    	return v;
	    }


		// set

		LMAT_ENSURE_INLINE void reset()
		{
			v = _mm_setzero_pd();
		}

		LMAT_ENSURE_INLINE void set(const double& ev)
		{
			v = _mm_set1_pd(ev);
		}

		LMAT_ENSURE_INLINE void set(const double& e0, const double& e1)
		{
			v = _mm_setr_pd(e0, e1);
		}


		// load

		LMAT_ENSURE_INLINE void load_u(const double *p)
		{
			v = _mm_loadu_pd(p);
		}

		LMAT_ENSURE_INLINE void load_a(const double *p)
		{
			v = _mm_load_pd(p);
		}

	    LMAT_ENSURE_INLINE void load_part(unsigned int n, double const * p)
	    {
	    	v = internal::sse_loadpart_f64(n, p);
	    }

	    // store

	    LMAT_ENSURE_INLINE void store_u(double *p) const
	    {
	    	_mm_storeu_pd(p, v);
	    }

	    LMAT_ENSURE_INLINE void store_a(double *p) const
	    {
	    	_mm_store_pd(p, v);
	    }

	    LMAT_ENSURE_INLINE void store_part(unsigned int n, double *p) const
	    {
	    	internal::sse_storepart_f64(n, p, v);
	    }


	    // extract

	    LMAT_ENSURE_INLINE double to_scalar() const
	    {
	    	return _mm_cvtsd_f64(v);
	    }

	    LMAT_ENSURE_INLINE double extract(unsigned int i) const
	    {
	    	return internal::sse_extract_f64(v, i);
	    }


	}; // SSE f64 pack


} }

#endif /* SSE_PACKS_H_ */
