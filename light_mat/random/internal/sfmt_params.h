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
			(int)param_t::MSK1, (int)param_t::MSK2,
			(int)param_t::MSK3, (int)param_t::MSK4);
	}

	// Note: we disable 607 and 216091, as their state contains odd number of packs
	// in other settings, the state contains even number of packs, which is good for AVX

/*
	template<>
	struct sfmt_params<607> : public sfmt_params_base<607>
	{
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
*/

	template<>
	struct sfmt_params<1279> : public sfmt_params_base<1279>
	{
		static const unsigned int POS1 = 7;
		static const unsigned int SL1 = 14;
		static const unsigned int SL2 = 3;
		static const unsigned int SR1 = 5;
		static const unsigned int SR2 = 1;
		static const unsigned int MSK1 = 0xf7fefffdU;
		static const unsigned int MSK2 = 0x7fefcfffU;
		static const unsigned int MSK3 = 0xaff3ef3fU;
		static const unsigned int MSK4 = 0xb5ffff7fU;
		static const unsigned int PARITY1 = 0x00000001U;
		static const unsigned int PARITY2 = 0x00000000U;
		static const unsigned int PARITY3 = 0x00000000U;
		static const unsigned int PARITY4 = 0x20000000U;
	};

	template<>
	struct sfmt_params<2281> : public sfmt_params_base<2281>
	{
		static const unsigned int POS1 = 12;
		static const unsigned int SL1 = 19;
		static const unsigned int SL2 = 1;
		static const unsigned int SR1 = 5;
		static const unsigned int SR2 = 1;
		static const unsigned int MSK1 = 0xbff7ffbfU;
		static const unsigned int MSK2 = 0xfdfffffeU;
		static const unsigned int MSK3 = 0xf7ffef7fU;
		static const unsigned int MSK4 = 0xf2f7cbbfU;
		static const unsigned int PARITY1 = 0x00000001U;
		static const unsigned int PARITY2 = 0x00000000U;
		static const unsigned int PARITY3 = 0x00000000U;
		static const unsigned int PARITY4 = 0x41dfa600U;
	};

	template<>
	struct sfmt_params<4253> : public sfmt_params_base<4253>
	{
		static const unsigned int POS1 = 17;
		static const unsigned int SL1 = 20;
		static const unsigned int SL2 = 1;
		static const unsigned int SR1 = 7;
		static const unsigned int SR2 = 1;
		static const unsigned int MSK1 = 0x9f7bffffU;
		static const unsigned int MSK2 = 0x9fffff5fU;
		static const unsigned int MSK3 = 0x3efffffbU;
		static const unsigned int MSK4 = 0xfffff7bbU;
		static const unsigned int PARITY1 = 0xa8000001U;
		static const unsigned int PARITY2 = 0xaf5390a3U;
		static const unsigned int PARITY3 = 0xb740b3f8U;
		static const unsigned int PARITY4 = 0x6c11486dU;
	};

	template<>
	struct sfmt_params<11213> : public sfmt_params_base<11213>
	{
		static const unsigned int POS1 = 68;
		static const unsigned int SL1 = 14;
		static const unsigned int SL2 = 3;
		static const unsigned int SR1 = 7;
		static const unsigned int SR2 = 3;
		static const unsigned int MSK1 = 0xeffff7fbU;
		static const unsigned int MSK2 = 0xffffffefU;
		static const unsigned int MSK3 = 0xdfdfbfffU;
		static const unsigned int MSK4 = 0x7fffdbfdU;
		static const unsigned int PARITY1 = 0x00000001U;
		static const unsigned int PARITY2 = 0x00000000U;
		static const unsigned int PARITY3 = 0xe8148000U;
		static const unsigned int PARITY4 = 0xd0c7afa3U;
	};

	template<>
	struct sfmt_params<19937> : public sfmt_params_base<19937>
	{
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
	struct sfmt_params<44497> : public sfmt_params_base<44497>
	{
		static const unsigned int POS1 = 330;
		static const unsigned int SL1 = 5;
		static const unsigned int SL2 = 3;
		static const unsigned int SR1 = 9;
		static const unsigned int SR2 = 3;
		static const unsigned int MSK1 = 0xeffffffbU;
		static const unsigned int MSK2 = 0xdfbebfffU;
		static const unsigned int MSK3 = 0xbfbf7befU;
		static const unsigned int MSK4 = 0x9ffd7bffU;
		static const unsigned int PARITY1 = 0x00000001U;
		static const unsigned int PARITY2 = 0x00000000U;
		static const unsigned int PARITY3 = 0xa3ac4000U;
		static const unsigned int PARITY4 = 0xecc1327aU;
	};

	template<>
	struct sfmt_params<86243> : public sfmt_params_base<86243>
	{
		static const unsigned int POS1 = 366;
		static const unsigned int SL1 = 6;
		static const unsigned int SL2 = 7;
		static const unsigned int SR1 = 19;
		static const unsigned int SR2 = 1;
		static const unsigned int MSK1 = 0xfdbffbffU;
		static const unsigned int MSK2 = 0xbff7ff3fU;
		static const unsigned int MSK3 = 0xfd77efffU;
		static const unsigned int MSK4 = 0xbf9ff3ffU;
		static const unsigned int PARITY1 = 0x00000001U;
		static const unsigned int PARITY2 = 0x00000000U;
		static const unsigned int PARITY3 = 0x00000000U;
		static const unsigned int PARITY4 = 0xe9528d85U;
	};

	template<>
	struct sfmt_params<132049> : public sfmt_params_base<132049>
	{
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

/*
	template<>
	struct sfmt_params<216091> : public sfmt_params_base<216091>
	{
		static const unsigned int POS1 = 627;
		static const unsigned int SL1 = 11;
		static const unsigned int SL2 = 3;
		static const unsigned int SR1 = 10;
		static const unsigned int SR2 = 1;
		static const unsigned int MSK1 = 0xbff7bff7U;
		static const unsigned int MSK2 = 0xbfffffffU;
		static const unsigned int MSK3 = 0xbffffa7fU;
		static const unsigned int MSK4 = 0xffddfbfbU;
		static const unsigned int PARITY1 = 0xf8000001U;
		static const unsigned int PARITY2 = 0x89e80709U;
		static const unsigned int PARITY3 = 0x3bd2b64bU;
		static const unsigned int PARITY4 = 0x0c64b1e4U;
	};
*/

} } }

#endif /* SFMT_PARAMS_H_ */
