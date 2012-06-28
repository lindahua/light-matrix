/*
 * @file matrix_check.h
 *
 * Functions for matrix argument checking
 *
 * @author Dahua Lin
 */

#ifndef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_CHECK_H_
#define LIGHTMAT_MATRIX_CHECK_H_

#include <light_mat/matrix/matrix_properties.h>

namespace lmat
{
	template<class Mat1, typename T1, class Mat2, typename T2>
	LMAT_ENSURE_INLINE
	inline void check_same_size(
			const IMatrixXpr<Mat1, T1>& A,
			const IMatrixXpr<Mat2, T2>& B,
			const char *msg)
	{
		check_arg(has_same_size(A, B), msg);
	}

	template<class Mat, typename T>
	LMAT_ENSURE_INLINE
	inline void check_square(
			const IMatrixXpr<Mat, T>& A,
			const char *msg)
	{
		check_arg(is_square(A), msg);
	}

	template<class Mat1, typename T1, class Mat2, typename T2>
	LMAT_ENSURE_INLINE
	inline void check_same_innerdim(
			const IMatrixXpr<Mat1, T1>& A,
			const IMatrixXpr<Mat2, T2>& B,
			const char *msg)
	{
		check_arg(A.ncolumns() == B.nrows(), msg);
	}
}


#endif 
