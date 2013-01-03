/**
 * @file int_div.h
 *
 * @brief Facilities for integer division
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_INT_DIV_H_
#define LIGHTMAT_INT_DIV_H_

#include <light_mat/common/prim_types.h>

#define LMAT_DEFINE_INT_DIV_C( L, D ) \
	template<> struct int_div<D> { \
		LMAT_ENSURE_INLINE \
		static size_t quo(size_t n) { return n >> L; } \
		LMAT_ENSURE_INLINE \
		static size_t rem(size_t n) { return n & (D-1); } \
		LMAT_ENSURE_INLINE \
		static size_t maj(size_t n) { return (n >> L) << L; } \
	};


namespace lmat
{

	template<unsigned int D> struct int_div;

	template<>
	struct int_div<1>
	{
		LMAT_ENSURE_INLINE
		static size_t quo(size_t n) { return n; }

		LMAT_ENSURE_INLINE
		static size_t rem(size_t n) { return 0; }

		LMAT_ENSURE_INLINE
		static size_t maj(size_t n) { return n; }
	};

	LMAT_DEFINE_INT_DIV_C(1, 2)
	LMAT_DEFINE_INT_DIV_C(2, 4)
	LMAT_DEFINE_INT_DIV_C(3, 8)
	LMAT_DEFINE_INT_DIV_C(4, 16)
	LMAT_DEFINE_INT_DIV_C(5, 32)
	LMAT_DEFINE_INT_DIV_C(6, 64)
	LMAT_DEFINE_INT_DIV_C(7, 128)
	LMAT_DEFINE_INT_DIV_C(8, 256)
	LMAT_DEFINE_INT_DIV_C(9, 512)

	LMAT_DEFINE_INT_DIV_C(10, 1024)


	namespace bdtags
	{
		struct dbl { };
		struct quad { };
		struct oct { };
		struct hex { };
	}
}

#endif /* INT_DIV_H_ */
