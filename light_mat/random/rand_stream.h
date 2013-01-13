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
#include <light_mat/simd/simd_base.h>

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


} }

#endif /* RAND_STREAM_H_ */
