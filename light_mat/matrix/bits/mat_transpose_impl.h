/**
 * @file mat_transpose_impl.h
 *
 * Implementation of 2D matrix transposition
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MAT_TRANSPOSE_IMPL_H_
#define LIGHTMAT_MAT_TRANSPOSE_IMPL_H_

#include <light_mat/core/mem_op.h>

namespace lmat { namespace detail {

	template<typename T>
	inline void direct_transpose(const index_t m, const index_t n,
			const T *src, const index_t ldim_s, T *dst, const index_t ldim_d)  // m x n --> n x m
	{
		if (m < n)
		{
			for (index_t i = 0; i < m; ++i)
			{
				pack_vec(n, src + i, ldim_s, dst + i * ldim_d);
			}
		}
		else
		{
			for (index_t j = 0; j < n; ++j)
			{
				unpack_vec(m, src + j * ldim_s, dst + j, ldim_d);
			}
		}
	}


	template<typename T>
	LMAT_ENSURE_INLINE
	inline void transpose(const index_t m, const index_t n,
			const T *src, const index_t ldim_s, T *dst, const index_t ldim_d)
	{
		direct_transpose(m, n, src, ldim_s, dst, ldim_d);
	}


} }

#endif /* MAT_TRANSPOSE_IMPL_H_ */
