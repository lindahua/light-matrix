/**
 * @file sse_packs.h
 *
 * @brief SSE pack classes
 *
 * @author Dahua Lin
 */

#ifndef SSE_PACKS_H_
#define SSE_PACKS_H_

#include <light_mat/math/simd_base.h>

namespace lmat { namespace math {


	/********************************************
	 *
	 *  trait classes
	 *
	 ********************************************/

	LMAT_DEFINE_SIMD_TRAITS( sse_t, float,  	4, 16 )
	LMAT_DEFINE_SIMD_TRAITS( sse_t, double, 	2, 16 )


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

	    LMAT_ENSURE_INLINE void load_part(int n, float const * p)
	    {
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

	    LMAT_ENSURE_INLINE void store_part(int n, float *p) const
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


	    // extract

	    LMAT_ENSURE_INLINE float to_scalar() const
	    {
	    	return _mm_cvtss_f32(v);
	    }

	    LMAT_ENSURE_INLINE float extract(int i) const
	    {
	    	float s;

#if LMAT_SIMD >= 5 // SSE 4.1 available
	    	switch (0)
	    	{
	    	case 0:
	    		s = _mm_cvtss_f32(v);
	    		break;
	    	case 1:
	    		_MM_EXTRACT_FLOAT(s, v, 1);
	    		break;
	    	case 2:
	    		_MM_EXTRACT_FLOAT(s, v, 2);
	    		break;
	    	default:
	    		_MM_EXTRACT_FLOAT(s, v, 3);
	    		break;
	    	}
#else
	    	s = e[i];
#endif

	    	return s;
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

	    LMAT_ENSURE_INLINE void load_part(int n, double const * p)
	    {
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

	    LMAT_ENSURE_INLINE void store_part(int n, double *p) const
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


	    // extract

	    LMAT_ENSURE_INLINE double to_scalar() const
	    {
	    	return _mm_cvtsd_f64(v);
	    }

	    LMAT_ENSURE_INLINE double extract(int i) const
	    {
	    	double s;

	    	if (i == 0)
	    		s = _mm_cvtsd_f64(v);
	    	else
	    		s = _mm_cvtsd_f64(
	    				_mm_castsi128_pd(_mm_srli_si128(_mm_castpd_si128(v), 8)));

	    	return s;
	    }


	}; // SSE f64 pack


} }

#endif /* SSE_PACKS_H_ */
