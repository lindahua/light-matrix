/**
 * @file simd_debug.h
 *
 * @brief Devices for SIMD debugging
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_SIMD_DEBUG_H_
#define LIGHTMAT_SIMD_DEBUG_H_

#include <light_mat/simd/simd_base.h>
#include <cstdio>

namespace lmat { namespace math {

	template<typename T, typename Kind>
	inline void print_pack(const char *fmt, const simd_pack<T, Kind>& pk)
	{
		const unsigned int width = pk.width();

		std::printf("(");

		for (unsigned int i = 0; i < width - 1; ++i)
		{
			std::printf(fmt, pk[i]);
			std::printf(", ");
		}
		std::printf(fmt, pk[width-1]);

		std::printf(")");
	}


	template<typename T, typename Kind>
	inline void print_pack(const char *fmt, const simd_bpack<T, Kind>& pk)
	{
		const unsigned int width = pk.width();

		std::printf("(");

		for (unsigned int i = 0; i < width - 1; ++i)
		{
			std::printf(fmt, pk[i]);
			std::printf(", ");
		}
		std::printf(fmt, pk[width-1]);

		std::printf(")");
	}


	template<typename T, typename Kind>
	inline bool test_equal(const simd_pack<T, Kind>& pk, const T& v)
	{
		const unsigned int width = pk.width();

		for (unsigned int i = 0; i < width; ++i)
		{
			if (!(pk[i] == v)) return false;
		}

		return true;
	}


	template<typename T, typename Kind>
	inline bool test_equal(const simd_pack<T, Kind>& pk, const T *ref)
	{
		const unsigned int width = pk.width();

		for (unsigned int i = 0; i < width; ++i)
		{
			if (!(pk[i] == ref[i])) return false;
		}

		return true;
	}

	template<typename T, typename Kind>
	inline bool test_equal(const simd_bpack<T, Kind>& pk, const typename simd_bpack<T, Kind>::bint_type& v)
	{
		const unsigned int width = pk.width();

		for (unsigned int i = 0; i < width; ++i)
		{
			if (!(pk[i] == v)) return false;
		}

		return true;
	}

	template<typename T, typename Kind>
	inline bool test_equal(const simd_bpack<T, Kind>& pk, const typename simd_bpack<T, Kind>::bint_type* ref)
	{
		const unsigned int width = pk.width();

		for (unsigned int i = 0; i < width; ++i)
		{
			if (!(pk[i] == ref[i])) return false;
		}

		return true;
	}

} }

#endif /* SIMD_DEBUG_H_ */
