/**
 * @file sfmt.h
 *
 * @brief SFMT random stream
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_SFMT_H_
#define LIGHTMAT_SFMT_H_

#include <light_mat/random/rand_stream.h>
#include <light_mat/random/stream_tracker.h>

#include "internal/sfmt_params.h"
#include "internal/rand_stream_internal.h"

#define LMAT_SFMT_IDXOF(i) i

namespace lmat { namespace random {

	/********************************************
	 *
	 *  sfmt_state
	 *
	 ********************************************/

	union sfmt_pack_t
	{
		__m128i si;
		uint64_t u64[2];
	    uint32_t u[4];
	};

	template<unsigned MEXP>
	class sfmt_state
	{
	public:
		typedef internal::sfmt_params<MEXP> param_t;

		LMAT_ENSURE_INLINE
		explicit sfmt_state(unsigned int seed=1234)
		{
			init_states(seed);
			param_mask = internal::sfmt_init_param_mask<MEXP>();
			pbase = state[0].u;
		}

		void init_states(uint32_t seed);

		void next()
		{
		    __m128i r1 = state[param_t::N - 2].si;
		    __m128i r2 = state[param_t::N - 1].si;

			unsigned int i;
		    for (i = 0; i < param_t::N - param_t::POS1; i++)
			{
				state[i].si = mm_recursion(state[i].si,
				     state[i + param_t::POS1].si, r1, r2, param_mask);
				r1 = r2;
				r2 = state[i].si;
		    }

		    for (; i < param_t::N; i++)
			{
				state[i].si = mm_recursion(state[i].si,
				     state[i + param_t::POS1 - param_t::N].si,
				     r1, r2, param_mask);
				r1 = r2;
				r2 = state[i].si;
		    }
		}

		const uint32_t* ptr_base() const
		{
			return pbase;
		}

		__m128i pack(size_t offset) const  // offset must be multiples of four
		{
			return _mm_load_si128(reinterpret_cast<const __m128i*>(pbase + offset));
		}

#ifdef LMAT_HAS_AVX
		__m256i avx_pack(size_t offset) const  // offset must be multiples of eight & at least eight u32 remain
		{
			return _mm256_load_si256(reinterpret_cast<const __m256i*>(pbase + offset));
		}
#endif

		uint64_t u64(size_t offset) const // offset must be multiples of two
		{
			return *(reinterpret_cast<const uint64_t*>(pbase + offset));
		}

		uint32_t u32(size_t offset) const
		{
			return pbase[offset];
		}

	private:

		LMAT_ENSURE_INLINE
		static __m128i mm_recursion(__m128i a, __m128i b,
						__m128i c, __m128i d, __m128i msk)
		{
		    __m128i v, x, y, z;

		    y = _mm_srli_epi32(b, param_t::SR1);
		    z = _mm_srli_si128(c, param_t::SR2);
		    v = _mm_slli_epi32(d, param_t::SL1);
		    z = _mm_xor_si128(z, a);
		    z = _mm_xor_si128(z, v);
		    x = _mm_slli_si128(a, param_t::SL2);
		    y = _mm_and_si128(y, msk);
		    z = _mm_xor_si128(z, x);
		    return _mm_xor_si128(z, y);
		}

	private:
		LMAT_ALIGN(64) sfmt_pack_t state[param_t::N];
		__m128i param_mask;
		uint32_t *pbase;
	};


	template<unsigned int MEXP>
	void sfmt_state<MEXP>::init_states(uint32_t seed)
	{
		const uint32_t parity[4] = {
			param_t::PARITY1, param_t::PARITY2,
			param_t::PARITY3, param_t::PARITY4};

	    uint32_t *psfmt32 = &state[0].u[0];

	    psfmt32[LMAT_SFMT_IDXOF(0)] = seed;
	    for (unsigned int i = 1; i < param_t::N32; i++)
		{
			psfmt32[LMAT_SFMT_IDXOF(i)] =
				1812433253UL * (psfmt32[LMAT_SFMT_IDXOF(i - 1)]
				^ (psfmt32[LMAT_SFMT_IDXOF(i - 1)] >> 30)) + i;
	    }

	    // period certification

	    int inner = 0;

	    for (unsigned int i = 0; i < 4; i++)
		{
			inner ^= psfmt32[LMAT_SFMT_IDXOF(i)] & parity[i];
		}

	    for (unsigned int i = 16; i > 0; i >>= 1)
			inner ^= inner >> i;

	    inner &= 1;

	    // check OK

		if (inner == 1) return;

	    // test and modify

	    for (unsigned int i = 0; i < 4; i++)
		{
			uint32_t work = 1;
			for (unsigned int j = 0; j < 32; j++)
			{
		    	if ((work & parity[i]) != 0)
				{
					psfmt32[LMAT_SFMT_IDXOF(i)] ^= work;
					return;
		    	}
		    	work = work << 1;
			}
	    }
	}


	/********************************************
	 *
	 *  sfmt_rand_stream
	 *
	 ********************************************/

	template<unsigned int MEXP>
	struct rand_stream_traits<sfmt_rand_stream<MEXP> >
	{
		typedef uint32_t seed_type;
	};


	template<unsigned int MEXP>
	class sfmt_rand_stream : public IRandStream<sfmt_rand_stream<MEXP> >
	{
		typedef internal::sfmt_params<MEXP> param_t;

	public:
		typedef uint32_t seed_type;

		LMAT_ENSURE_INLINE
		sfmt_rand_stream(uint32_t seed=1234)
		: m_intern(seed), m_tracker(param_t::N32) { }

		LMAT_ENSURE_INLINE
		unsigned int period_exponent() const
		{
			return MEXP;
		}

		LMAT_ENSURE_INLINE
		unsigned int internal_npacks() const
		{
			return param_t::N;
		}

		LMAT_ENSURE_INLINE
		unsigned int internal_nunits() const
		{
			return (unsigned int)m_tracker.length();
		}

	public:
		LMAT_ENSURE_INLINE
		void set_seed(const seed_type& seed)
		{
			m_intern.init_states(seed);
			m_tracker.set_end();
		}

		LMAT_ENSURE_INLINE
		size_t state_size() const  // in terms of bytes
		{
			return param_t::N * 16;
		}

		LMAT_ENSURE_INLINE uint32_t rand_u32()
		{
			check_end();
			uint32_t x = m_intern.u32(m_tracker.offset());
			m_tracker.forward(1);
			return x;
		}

		LMAT_ENSURE_INLINE uint64_t rand_u64()
		{
			m_tracker.to_boundary(bdtags::dbl());

			check_end();
			uint64_t x = m_intern.u64(m_tracker.offset());
			m_tracker.forward(2);
			return x;
		}

		LMAT_ENSURE_INLINE __m128i rand_pack()
		{
			m_tracker.to_boundary(bdtags::quad());

			check_end();
			__m128i u = m_intern.pack(m_tracker.offset());
			m_tracker.forward(4);
			return u;
		}

#ifdef LMAT_HAS_AVX
		LMAT_ENSURE_INLINE __m256i rand_avx_pack()
		{
			m_tracker.to_boundary(bdtags::oct());

			check_end();  // note: it was assured that the state has even number of 128-bit packs
			__m256i u = m_intern.avx_pack(m_tracker.offset());
			m_tracker.forward(8);
			return u;
		}
#endif

		LMAT_ENSURE_INLINE void rand_seq(size_t nbytes, void *buf)
		{
			internal::gen_rand_seq(m_intern, m_tracker, buf, nbytes);
		}

	private:
		LMAT_ENSURE_INLINE
		void check_end()
		{
			if (m_tracker.is_end())
			{
				m_intern.next();
				m_tracker.rewind();
			}
		}

	private:
		sfmt_state<MEXP> m_intern;
		stream_tracker<uint32_t> m_tracker;
	};


	// Typdefs

	typedef sfmt_rand_stream<607>    sfmt607_t;
	typedef sfmt_rand_stream<1279>   sfmt1279_t;
	typedef sfmt_rand_stream<2281>   sfmt2281_t;
	typedef sfmt_rand_stream<4253>   sfmt4253_t;
	typedef sfmt_rand_stream<11213>  sfmt11213_t;
	typedef sfmt_rand_stream<19937>  sfmt19937_t;
	typedef sfmt_rand_stream<44497>  sfmt44497_t;
	typedef sfmt_rand_stream<86243>  sfmt86243_t;
	typedef sfmt_rand_stream<132049> sfmt132049_t;
	typedef sfmt_rand_stream<216091> sfmt216091_t;

} }

#endif /* SFMT_H_ */
