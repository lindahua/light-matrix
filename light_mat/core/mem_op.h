/**
 * @file mem_op.h
 *
 * @brief Functions for memory operations.
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MEM_OP_H_
#define LIGHTMAT_MEM_OP_H_

#include <light_mat/core/basic_defs.h>
#include <cstring>

namespace lmat
{
	template<typename T>
	LMAT_ENSURE_INLINE inline size_t nbytes(size_t n)
	{
		return static_cast<size_t>(n) * sizeof(T);
	}

	template<typename T>
	LMAT_ENSURE_INLINE
	inline void copy_mem(const index_t n, const T *src, T *dst)
	{
		std::memcpy(dst, src, nbytes<T>(n));
	}

	template<typename T>
	LMAT_ENSURE_INLINE
	inline void zero_mem(const index_t n, T *dst)
	{
		std::memset(dst, 0, nbytes<T>(n));
	}

	template<typename T>
	LMAT_ENSURE_INLINE
	inline void fill_mem(const index_t n, T *dst, const T& v)
	{
		for (index_t i = 0; i < n; ++i) dst[i] = v;
	}

	template<typename T>
	LMAT_ENSURE_INLINE
	inline void mem_equal(const index_t n, const T *a, const T *b)
	{
		return std::memcmp(a, b, nbytes<T>(n)) == 0;
	}

	template<typename T>
	LMAT_ENSURE_INLINE
	inline void mem_equal(const index_t n, const T *a, const T& v)
	{
		for (index_t i = 0; i < n; ++i)
			if (!(a[i] == v)) return false;
		return true;
	}

	template<typename T>
	LMAT_ENSURE_INLINE
	inline void unpack_vec(const index_t n, const T* a, T* b, const index_t bstep)
	{
		for (index_t i = 0; i < n; ++i) b[i * bstep] = a[i];
	}

	template<typename T>
	LMAT_ENSURE_INLINE
	inline void pack_vec(const index_t n, const T* a, const index_t astep, T* b)
	{
		for (index_t i = 0; i < n; ++i) b[i] = a[i * astep];
	}

}

#endif /* MEM_OP_H_ */
