/**
 * @file avx_base.h
 *
 * @brief The AVX pack classes
 *
 * @author Dahua Lin
 */

#ifndef LIGHTMAT_AVX_PACKS_H_
#define LIGHTMAT_AVX_PACKS_H_

#include <light_mat/math/simd_base.h>
#include "internal/avx_helpers.h"

#ifndef LMAT_HAS_AVX
#error Only include avx_packs.h when AVX is enabled.
#endif

namespace lmat { namespace math {


	/********************************************
	 *
	 *  trait classes
	 *
	 ********************************************/

	LMAT_DEFINE_SIMD_TRAITS( avx_t, float,  8, 32 )
	LMAT_DEFINE_SIMD_TRAITS( avx_t, double, 4, 32 )


	/********************************************
	 *
	 *  pack classes
	 *
	 ********************************************/

	typedef simd_pack<float,  avx_t> avx_f32pk;
	typedef simd_pack<double, avx_t> avx_f64pk;


	template<>
	struct simd_pack<float, avx_t>
	{
		LMAT_DEFINE_FOR_SIMD_PACK( avx_t, float, 8 )

		union
		{
			__m256 v;
			LMAT_ALIGN_AVX float e[8];
		};

		LMAT_ENSURE_INLINE
		unsigned int width() const
		{
			return pack_width;
		}


		// constructors

		LMAT_ENSURE_INLINE simd_pack() { }

		LMAT_ENSURE_INLINE simd_pack(const __m256& v_) : v(v_) { }

		LMAT_ENSURE_INLINE simd_pack(const float& ev)
		{
			v = _mm256_set1_ps(ev);
		}

		LMAT_ENSURE_INLINE simd_pack(
				const float& e0, const float& e1, const float& e2, const float& e3,
				const float& e4, const float& e5, const float& e6, const float& e7)
		{
			v = _mm256_setr_ps(e0, e1, e2, e3, e4, e5, e6, e7);
		}

	    LMAT_ENSURE_INLINE
	    static simd_pack zeros()
	    {
	    	return _mm256_setzero_ps();
	    }

	    LMAT_ENSURE_INLINE
	    static simd_pack ones()
	    {
	    	return _mm256_set1_ps(1.0f);
	    }

	    LMAT_ENSURE_INLINE
	    static simd_pack inf()
	    {
	    	return _mm256_castsi256_ps(_mm256_set1_epi32((int)0x7f800000));
	    }

	    LMAT_ENSURE_INLINE
	    static simd_pack neg_inf()
	    {
	    	return _mm256_castsi256_ps(_mm256_set1_epi32((int)0xff800000));
	    }

	    LMAT_ENSURE_INLINE
	    static simd_pack nan()
	    {
	    	return _mm256_set1_ps(std::numeric_limits<float>::quiet_NaN());
	    }


	    // converter

	    LMAT_ENSURE_INLINE
	    operator __m256() const
	    {
	    	return v;
	    }


		// set

		LMAT_ENSURE_INLINE void reset()
		{
			v = _mm256_setzero_ps();
		}

		LMAT_ENSURE_INLINE void set(const float& ev)
		{
			v = _mm256_set1_ps(ev);
		}

		LMAT_ENSURE_INLINE void set(
				const float& e0, const float& e1, const float& e2, const float& e3,
				const float& e4, const float& e5, const float& e6, const float& e7)
		{
			v = _mm256_setr_ps(e0, e1, e2, e3, e4, e5, e6, e7);
		}


		// load

		LMAT_ENSURE_INLINE void load_u(const float *p)
		{
			v = _mm256_loadu_ps(p);
		}

		LMAT_ENSURE_INLINE void load_a(const float *p)
		{
			v = _mm256_load_ps(p);
		}

	    LMAT_ENSURE_INLINE void load_part(unsigned int n, float const * p)
	    {
	    	v = _mm256_maskload_ps(p, internal::avx_part_mask_32(n));
	    }

	    // store

	    LMAT_ENSURE_INLINE void store_u(float *p) const
	    {
	    	_mm256_storeu_ps(p, v);
	    }

	    LMAT_ENSURE_INLINE void store_a(float *p) const
	    {
	    	_mm256_store_ps(p, v);
	    }

	    LMAT_ENSURE_INLINE void store_part(unsigned int n, float *p) const
	    {
	    	_mm256_maskstore_ps(p, internal::avx_part_mask_32(n), v);
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

	    LMAT_ENSURE_INLINE float to_scalar() const
	    {
	    	return _mm_cvtss_f32(get_low());
	    }

	    LMAT_ENSURE_INLINE float extract(unsigned int i) const
	    {
	    	return e[i];
	    }

	}; // AVX f32 pack


	template<>
	struct simd_pack<double, avx_t>
	{
		LMAT_DEFINE_FOR_SIMD_PACK( avx_t, double, 4 )

		union
		{
			__m256d v;
			LMAT_ALIGN_AVX double e[4];
		};

		LMAT_ENSURE_INLINE
		unsigned int width() const
		{
			return pack_width;
		}


		// constructors

		LMAT_ENSURE_INLINE simd_pack() { }

		LMAT_ENSURE_INLINE simd_pack(const __m256d& v_) : v(v_) { }

		LMAT_ENSURE_INLINE simd_pack(const double& ev)
		{
			v = _mm256_set1_pd(ev);
		}

		LMAT_ENSURE_INLINE simd_pack(
				const double& e0, const double& e1, const double& e2, const double& e3)
		{
			v = _mm256_setr_pd(e0, e1, e2, e3);
		}

	    LMAT_ENSURE_INLINE
	    static simd_pack zeros()
	    {
	    	return _mm256_setzero_pd();
	    }

	    LMAT_ENSURE_INLINE
	    static simd_pack ones()
	    {
	    	return _mm256_set1_pd(1.0);
	    }

	    LMAT_ENSURE_INLINE
	    static simd_pack inf()
	    {
	    	return _mm256_castsi256_pd(
	    			_mm256_set1_epi64x((int64_t)0x7ff0000000000000LL));
	    }

	    LMAT_ENSURE_INLINE
	    static simd_pack neg_inf()
	    {
	    	return _mm256_castsi256_pd(
	    			_mm256_set1_epi64x((int64_t)0xfff0000000000000LL));
	    }

	    LMAT_ENSURE_INLINE
	    static simd_pack nan()
	    {
	    	return _mm256_set1_pd(std::numeric_limits<double>::quiet_NaN());
	    }


	    // converters

	    LMAT_ENSURE_INLINE
	    operator __m256d() const
	    {
	    	return v;
	    }


		// set

		LMAT_ENSURE_INLINE void reset()
		{
			v = _mm256_setzero_pd();
		}

		LMAT_ENSURE_INLINE void set(const double& ev)
		{
			v = _mm256_set1_pd(ev);
		}

		LMAT_ENSURE_INLINE void set(
				const double& e0, const double& e1, const double& e2, const double& e3)
		{
			v = _mm256_setr_pd(e0, e1, e2, e3);
		}


		// load

		LMAT_ENSURE_INLINE void load_u(const double *p)
		{
			v = _mm256_loadu_pd(p);
		}

		LMAT_ENSURE_INLINE void load_a(const double *p)
		{
			v = _mm256_load_pd(p);
		}

	    LMAT_ENSURE_INLINE void load_part(unsigned int n, double const * p)
	    {
	    	v = _mm256_maskload_pd(p, internal::avx_part_mask_64(n));
	    }

	    // store

	    LMAT_ENSURE_INLINE void store_u(double *p) const
	    {
	    	_mm256_storeu_pd(p, v);
	    }

	    LMAT_ENSURE_INLINE void store_a(double *p) const
	    {
	    	_mm256_store_pd(p, v);
	    }

	    LMAT_ENSURE_INLINE void store_part(unsigned int n, double *p) const
	    {
	    	_mm256_maskstore_pd(p, internal::avx_part_mask_64(n), v);
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

	    LMAT_ENSURE_INLINE double to_scalar() const
	    {
	    	return _mm_cvtsd_f64(get_low());
	    }

	    LMAT_ENSURE_INLINE double extract(unsigned int i) const
	    {
	    	return e[i];
	    }


	}; // AVX f64 pack


} }


#endif
