/**
 * @file numrepr_format.h
 *
 * @brief Information about number representation formats
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_NUMREPR_FORMAT_H_
#define LIGHTMAT_NUMREPR_FORMAT_H_

#include <light_mat/common/prim_types.h>

namespace lmat { namespace internal {

	template<typename T> struct num_fmt;

	template<>
	struct num_fmt<float>
	{
		static const unsigned int nbytes = 4;
		static const unsigned int nbits = 32;

		typedef int32_t  sint_type;
		typedef uint32_t uint_type;

		static const unsigned int n_exponent_bits = 8;
		static const unsigned int n_mantissa_bits = 23;

		static const sint_type sign_bit      = (sint_type)0x80000000;
		static const sint_type exponent_bits = (sint_type)0x7f800000;
		static const sint_type mantissa_bits = (sint_type)0x007fffff;
	};


	template<>
	struct num_fmt<double>
	{
		static const unsigned int nbytes = 8;
		static const unsigned int nbits = 64;

		typedef int64_t  sint_type;
		typedef uint64_t uint_type;

		static const unsigned int n_exponent_bits = 11;
		static const unsigned int n_mantissa_bits = 52;

		static const sint_type sign_bit      = (sint_type)0x8000000000000000LL;
		static const sint_type exponent_bits = (sint_type)0x7ff0000000000000LL;
		static const sint_type mantissa_bits = (sint_type)0x000fffffffffffffLL;
	};


} }

#endif /* NUMREPR_FORMAT_H_ */
