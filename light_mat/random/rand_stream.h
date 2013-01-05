/**
 * @file rand_stream.h
 *
 * @brief The concept of random streams
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_RAND_STREAM_H_
#define LIGHTMAT_RAND_STREAM_H_

#include <light_mat/common/basic_defs.h>
#include <light_mat/math/simd_base.h>

namespace lmat { namespace random {

	/********************************************
	 *
	 *  The concept of random streams
	 *
	 ********************************************/

	template<class RStream>
	struct rand_stream_traits;

	template<class Derived>
	class IRandStream
	{
	public:
		LMAT_CRTP_REF
		typedef typename rand_stream_traits<Derived>::seed_type seed_type;

		LMAT_ENSURE_INLINE
		void set_seed(const seed_type& seed)
		{
			derived().set_seed(seed);
		}

		LMAT_ENSURE_INLINE
		size_t state_size() const  // in terms of bytes
		{
			return derived().state_size();
		}

		LMAT_ENSURE_INLINE
		uint32_t rand_u32()
		{
			return derived().rand_u32();
		}

		LMAT_ENSURE_INLINE
		uint64_t rand_u64()
		{
			return derived().rand_u64();
		}

		LMAT_ENSURE_INLINE
		void rand_seq(size_t nbytes, void *buf)
		{
			derived().rand_seq(nbytes, buf);
		}
	};


	/********************************************
	 *
	 *  pre-defined random stream classes
	 *
	 ********************************************/

	template<unsigned int MEXP> class sfmt_rand_stream;

	typedef sfmt_rand_stream<19937> default_rand_stream;


	/********************************************
	 *
	 *  random number conversion
	 *
	 ********************************************/

	// from random bits to real value in [1, 2)

	LMAT_ENSURE_INLINE
	inline float randbits_to_c1o2_f32(uint32_t u)
	{
		u = (u & 0x007fffff) | (0x3f800000);
		return float(reinterpret_cast<const float&>(u));
	}

	LMAT_ENSURE_INLINE
	inline double randbits_to_c1o2_f64(uint64_t u)
	{
		u = (u & 0x000fffffffffffffULL) | (0x3ff0000000000000);
		return double(reinterpret_cast<const double&>(u));
	}

	LMAT_ENSURE_INLINE
	inline __m128 randbits_to_c1o2_f32(const __m128i& u, sse_t)
	{
		return _mm_or_ps(
			_mm_set1_ps(1.0f),
			_mm_andnot_ps(_mm_set1_ps(-2.0f), _mm_castsi128_ps(u)));
	}

	LMAT_ENSURE_INLINE
	inline __m128d randbits_to_c1o2_f64(const __m128i& u, sse_t)
	{
		return _mm_or_pd(
			_mm_set1_pd(1.0),
			_mm_andnot_pd(_mm_set1_pd(-2.0), _mm_castsi128_pd(u)));
	}

#ifdef LMAT_HAS_AVX

	LMAT_ENSURE_INLINE
	inline __m256 randbits_to_c1o2_f32(const __m256i& u, avx_t)
	{
		return _mm256_or_ps(
			_mm256_set1_ps(1.0f),
			_mm256_andnot_ps(_mm256_set1_ps(-2.0f), _mm256_castsi256_ps(u)));
	}

	LMAT_ENSURE_INLINE
	inline __m256d randbits_to_c1o2_f64(const __m256i& u, avx_t)
	{
		return _mm256_or_pd(
			_mm256_set1_pd(1.0),
			_mm256_andnot_pd(_mm256_set1_pd(-2.0), _mm256_castsi256_pd(u)));
	}

#endif

} }

#endif /* RAND_STREAM_H_ */
