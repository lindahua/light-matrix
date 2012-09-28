/*
 * @file matrix_properties.h
 *
 * Basic properties of matrices
 *
 * @author Dahua Lin
 */

#ifndef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_PROPERTIES_H_
#define LIGHTMAT_MATRIX_PROPERTIES_H_

#include <light_mat/matrix/matrix_concepts.h>

#ifdef LMAT_ENABLE_DIM_CHECKING
#define LMAT_CHECK_SAME_SHAPE(a, b) check_arg(has_same_size((a), (b)), "Inconsistent matrix size.");
#else
#define LMAT_CHECK_SAME_SHAPE(a, b)
#endif

namespace lmat
{
	template<class Mat, typename T>
	LMAT_ENSURE_INLINE
	inline bool is_empty(const IMatrixXpr<Mat, T>& X)
	{
		return (has_dynamic_nrows<Mat>::value && X.nrows() == 0) ||
			(has_dynamic_ncols<Mat>::value && X.ncolumns() == 0);
	}

	template<class Mat, typename T>
	LMAT_ENSURE_INLINE
	inline bool is_column(const IMatrixXpr<Mat, T>& X)
	{
		return ct_is_col<Mat>::value || X.ncolumns() == 1;
	}

	template<class Mat, typename T>
	LMAT_ENSURE_INLINE
	inline bool is_row(const IMatrixXpr<Mat, T>& X)
	{
		return ct_is_row<Mat>::value || X.nrows() == 1;
	}

	template<class Mat, typename T>
	LMAT_ENSURE_INLINE
	inline bool is_scalar(const IMatrixXpr<Mat, T>& X)
	{
		return is_column(X) && is_row(X);
	}

	template<class Mat, typename T>
	LMAT_ENSURE_INLINE
	inline bool is_vector(const IMatrixXpr<Mat, T>& X)
	{
		return is_column(X) || is_row(X);
	}

	template<class Mat, typename T>
	LMAT_ENSURE_INLINE
	inline bool is_square(const IMatrixXpr<Mat, T>& X)
	{
		return (ct_rows<Mat>::value == ct_cols<Mat>::value) ||
				(X.nrows() == X.ncolumns());
	}

	template<class Mat1, typename T1, class Mat2, typename T2>
	LMAT_ENSURE_INLINE
	inline bool has_same_nrows(const IMatrixXpr<Mat1, T1>& A, const IMatrixXpr<Mat2, T2>& B)
	{
		return ct_has_same_nrows<Mat1, Mat2>::value || A.nrows() == B.nrows();
	}

	template<class Mat1, typename T1, class Mat2, typename T2>
	LMAT_ENSURE_INLINE
	inline bool has_same_ncolumns(const IMatrixXpr<Mat1, T1>& A, const IMatrixXpr<Mat2, T2>& B)
	{
		return ct_has_same_ncols<Mat1, Mat2>::value || A.ncolumns() == B.ncolumns();
	}

	template<class Mat1, typename T1, class Mat2, typename T2>
	LMAT_ENSURE_INLINE
	inline bool has_same_size(const IMatrixXpr<Mat1, T1>& A, const IMatrixXpr<Mat2, T2>& B)
	{
		return has_same_nrows(A, B) && has_same_ncolumns(A, B);
	}

	template<class Mat, typename T>
	LMAT_ENSURE_INLINE
	inline bool has_continuous_layout(const IDenseMatrix<Mat, T>& A)
	{
		return ct_has_continuous_layout<Mat>::value || A.lead_dim() == A.nrows();
	}
}

#endif 
