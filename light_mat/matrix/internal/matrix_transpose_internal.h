/**
 * @file matrix_transpose_internal.h
 *
 * @brief Internal implementation of matrix transposition
 *
 * @author Dahua Lin
 */

#ifndef LIGHTMAT_MATRIX_TRANSPOSE_INTERNAL_H_
#define LIGHTMAT_MATRIX_TRANSPOSE_INTERNAL_H_

#include <light_mat/common/memory.h>
#include <light_mat/matrix/matrix_classes.h>

namespace lmat { namespace internal {

	template<typename T>
	inline void naive_transpose(index_t m, index_t n, const T *src, T *dst)
	{
		if (m == 1 || n == 1)
		{
			copy_vec(m * n, src, dst);
		}
		else if (m >= n) // cols --> rows
		{
			for (index_t j = 0; j < n; ++j)
			{
				copy_vec(m, src + j * m, step_ptr(dst + j, n));
			}
		}
		else // rows --> cols
		{
			for (index_t i = 0; i < m; ++i)
			{
				copy_vec(n, step_ptr(src + i, m), dst + i * n);
			}
		}
	}


	template<typename T>
	inline void naive_transpose(index_t m, index_t n,
			const T* src, index_t src_cs, T *dst, index_t dst_cs)
	{
		if (m == 1)  // 1 x n --> n x 1
		{
			if (src_cs == 1)
				copy_vec(n, src, dst);
			else
				copy_vec(n, step_ptr(src, src_cs), dst);
		}
		else if (n == 1) // m x 1 --> 1 x m
		{
			if (dst_cs == 1)
				copy_vec(m, src, dst);
			else
				copy_vec(m, src, step_ptr(dst, dst_cs));
		}
		else if (m >= n) // cols --> rows
		{
			for (index_t j = 0; j < n; ++j)
				copy_vec(m, src + j * src_cs, step_ptr(dst + j, dst_cs));
		}
		else // rows --> cols
		{
			for (index_t i = 0; i < m; ++i)
				copy_vec(n, step_ptr(src + i, src_cs), dst + i * dst_cs);
		}
	}


	template<typename T>
	inline void naive_transpose(index_t m, index_t n,
			const T* src, index_t src_rs, index_t src_cs,
			T *dst, index_t dst_rs, index_t dst_cs)
	{
		if (src_rs == 1 && dst_rs == 1)
		{
			naive_transpose(m, n, src, src_cs, dst, dst_cs);
		}
		else
		{
			if (m == 1) // 1 x n --> n x 1
			{
				copy_vec(n, step_ptr(src, src_cs), step_ptr(dst, dst_rs));
			}
			else if (n == 1) // m x 1 --> 1 x m
			{
				copy_vec(m, step_ptr(src, src_rs), step_ptr(dst, dst_cs));
			}
			else if (m >= n) // cols --> rows
			{
				for (index_t j = 0; j < n; ++j)
					copy_vec(m, step_ptr(src + j * src_cs, src_rs), step_ptr(dst + j * dst_rs, dst_cs));
			}
			else // rows --> cols
			{
				for (index_t i = 0; i < m; ++i)
					copy_vec(n, step_ptr(src + i * src_rs, src_cs), step_ptr(dst + i * dst_cs, dst_rs));
			}
		}
	}


	template<typename T, class SMat, class DMat>
	inline void direct_transpose(index_t m, index_t n, const IRegularMatrix<SMat, T>& smat, IRegularMatrix<DMat, T>& dmat)
	{
		if (meta::is_continuous<SMat>::value && meta::is_continuous<DMat>::value)
		{
			naive_transpose(m, n, smat.ptr_data(), dmat.ptr_data());
		}
		else if (meta::is_percol_continuous<SMat>::value && meta::is_percol_continuous<DMat>::value)
		{
			naive_transpose(m, n, smat.ptr_data(), smat.col_stride(), dmat.ptr_data(), dmat.col_stride());
		}
		else
		{
			naive_transpose(m, n,
					smat.ptr_data(), smat.row_stride(), smat.col_stride(),
					dmat.ptr_data(), dmat.row_stride(), dmat.col_stride());
		}
	}


} }

#endif /* MATRIX_TRANSPOSE_INTERNAL_H_ */
