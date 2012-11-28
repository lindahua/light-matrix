/**
 * @file matrix_transpose.h
 *
 * Matrix transposition
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_TRANSPOSE_H_
#define LIGHTMAT_MATRIX_TRANSPOSE_H_

#include <light_mat/matexpr/matrix_ewise_eval.h>
#include "bits/direct_transpose_impl.h"

namespace lmat
{

	/********************************************
	 *
	 *  Direct transpose
	 *
	 ********************************************/

	template<typename T, class SMat, class DMat>
	LMAT_ENSURE_INLINE
	inline void direct_transpose(const IDenseMatrix<SMat, T>& src, IDenseMatrix<DMat, T>& dst)
	{
		const index_t m = src.nrows();
		const index_t n = src.ncolumns();

		detail::transpose(m, n,
				src.ptr_data(), src.row_stride(), src.col_stride(),
				dst.ptr_data(), dst.row_stride(), dst.col_stride());
	}


	/********************************************
	 *
	 *  Transpose expression
	 *
	 ********************************************/

	template<typename Arg_HP, class Arg> class transpose_expr;


}

#endif /* MATRIX_TRANSPOSE_H_ */
