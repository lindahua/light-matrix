/**
 * @file sfmt_params.h
 *
 * @brief SFMT PRNG parameters
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_SFMT_PARAMS_H_
#define LIGHTMAT_SFMT_PARAMS_H_

#include <light_mat/common/prim_types.h>

namespace lmat { namespace random { namespace internal {

	template<unsigned int MEXP>
	struct sfmt_params_base
	{
		static const unsigned int N = (MEXP / 128 + 1);
		static const unsigned int N32 = N * 4;
		static const unsigned int N64 = N * 2;
	};

	template<unsigned int MEXP> struct sfmt_params;

	template<unsigned int MEXP>
	LMAT_ENSURE_INLINE
	inline __m128i sfmt_init_param_mask()
	{
		typedef sfmt_params<MEXP> param_t;
		return _mm_setr_epi32(
			param_t::MSK1, param_t::MSK2,
			param_t::MSK3, param_t::MSK4);
	}

	template<>
	struct sfmt_params<607> : public sfmt_params_base<607>
	{
		static const char * id_str()
		{
			return "SFMT-607:2-15-3-13-3:fdff37ff-ef7f3f7d-ff777b7d-7ff7fb2f";
		}

		static const unsigned int POS1 = 2;
		static const unsigned int SL1 = 15;
		static const unsigned int SL2 = 3;
		static const unsigned int SR1 = 13;
		static const unsigned int SR2 = 3;
		static const unsigned int MSK1 = 0xfdff37ffU;
		static const unsigned int MSK2 = 0xef7f3f7dU;
		static const unsigned int MSK3 = 0xff777b7dU;
		static const unsigned int MSK4 = 0x7ff7fb2fU;
		static const unsigned int PARITY1 = 0x00000001U;
		static const unsigned int PARITY2 = 0x00000000U;
		static const unsigned int PARITY3 = 0x00000000U;
		static const unsigned int PARITY4 = 0x5986f054U;
	};

	template<>
	struct sfmt_params<19937> : public sfmt_params_base<19937>
	{
		static const char * id_str()
		{
			return "SFMT-19937:122-18-1-11-1:dfffffef-ddfecb7f-bffaffff-bffffff6";
		}

		static const unsigned int POS1 = 122;
		static const unsigned int SL1 = 18;
		static const unsigned int SL2 = 1;
		static const unsigned int SR1 = 11;
		static const unsigned int SR2 = 1;
		static const unsigned int MSK1 = 0xdfffffefU;
		static const unsigned int MSK2 = 0xddfecb7fU;
		static const unsigned int MSK3 = 0xbffaffffU;
		static const unsigned int MSK4 = 0xbffffff6U;
		static const unsigned int PARITY1 = 0x00000001U;
		static const unsigned int PARITY2 = 0x00000000U;
		static const unsigned int PARITY3 = 0x00000000U;
		static const unsigned int PARITY4 = 0x13c9e684U;
	};

	template<>
	struct sfmt_params<132049> : public sfmt_params_base<132049>
	{
		static const char * id_str()
		{
			return "SFMT-132049:110-19-1-21-1:ffffbb5f-fb6ebf95-fffefffa-cff77fff";
		}

		static const unsigned int POS1 = 110;
		static const unsigned int SL1 = 19;
		static const unsigned int SL2 = 1;
		static const unsigned int SR1 = 21;
		static const unsigned int SR2 = 1;
		static const unsigned int MSK1 = 0xffffbb5fU;
		static const unsigned int MSK2 = 0xfb6ebf95U;
		static const unsigned int MSK3 = 0xfffefffaU;
		static const unsigned int MSK4 = 0xcff77fffU;
		static const unsigned int PARITY1 = 0x00000001U;
		static const unsigned int PARITY2 = 0x00000000U;
		static const unsigned int PARITY3 = 0xcb520000U;
		static const unsigned int PARITY4 = 0xc7e91c7dU;
	};


} } }

#endif /* SFMT_PARAMS_H_ */
