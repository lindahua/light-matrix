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

#include <light_mat/common/memory.h>

namespace lmat { namespace internal {

	template<typename T>
	LMAT_ENSURE_INLINE
	inline void vec_transpose(index_t len, const T *src, index_t src_step, T *dst, index_t dst_step)
	{
		if (src_step == 1)
		{
			if (dst_step == 1)
				copy_vec(len, src, dst);
			else
				copy_vec(len, src, dst, dst_step);
		}
		else
		{
			if (dst_step == 1)
				copy_vec(len, src, src_step, dst);
			else
				copy_vec(len, src, src_step, dst, dst_step);
		}
	}

	template<typename T>
	inline void simple_transpose(index_t m, index_t n,
			const T *src, index_t src_rs, index_t src_cs,
			T *dst, index_t dst_rs, index_t dst_cs)  // m x n --> n x m
	{
		if (m < n)  // src rows ==> dst columns
		{
			if (dst_rs == 1)
			{
				for (index_t i = 0; i < m; ++i)
					copy_vec(n, src + i * src_rs, src_cs, dst + i * dst_cs);
			}
			else
			{
				for (index_t i = 0; i < m; ++i)
					copy_vec(n, src + i * src_rs, src_cs, dst + i * dst_cs, dst_rs);
			}
		}
		else  // src columns ==> dst rows
		{
			if (src_rs == 1)
			{
				for (index_t j = 0; j < n; ++j)
					copy_vec(m, src + j * src_cs, dst + j * dst_rs, dst_cs);
			}
			else
			{
				for (index_t j = 0; j < n; ++j)
					copy_vec(m, src + j * src_cs, src_rs, dst + j * dst_rs, dst_cs);
			}
		}
	}


	template<typename T>
	LMAT_ENSURE_INLINE
	inline void transpose(index_t m, index_t n,
			const T *src, index_t src_rs, index_t src_cs,
			T *dst, index_t dst_rs, index_t dst_cs)
	{
		simple_transpose(m, n, src, src_rs, src_cs, dst, dst_rs, dst_cs);
	}


} }

#endif /* MAT_TRANSPOSE_IMPL_H_ */
