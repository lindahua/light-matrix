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
#define LMAT_CHECK_DIMS( cond ) check_arg( cond , "Inconsistent matrix dimensions");
#else
#define LMAT_CHECK_DIMS(n1, n2)
#endif

namespace lmat
{
	template<class Mat, typename T>
	LMAT_ENSURE_INLINE
	inline bool is_empty(const IMatrixXpr<Mat, T>& X)
	{
		return (meta::has_dynamic_nrows<Mat>::value && X.nrows() == 0) ||
			(meta::has_dynamic_ncols<Mat>::value && X.ncolumns() == 0);
	}

	template<class Mat, typename T>
	LMAT_ENSURE_INLINE
	inline bool is_column(const IMatrixXpr<Mat, T>& X)
	{
		return meta::is_col<Mat>::value || X.ncolumns() == 1;
	}

	template<class Mat, typename T>
	LMAT_ENSURE_INLINE
	inline bool is_row(const IMatrixXpr<Mat, T>& X)
	{
		return meta::is_row<Mat>::value || X.nrows() == 1;
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
		return X.nrows() == X.ncolumns();
	}

	template<class Mat1, typename T1, class Mat2, typename T2>
	LMAT_ENSURE_INLINE
	inline bool have_same_nrows(const IMatrixXpr<Mat1, T1>& A, const IMatrixXpr<Mat2, T2>& B)
	{
		return all_equal(A.nrows(), B.nrows());
	}

	template<class Mat1, typename T1, class Mat2, typename T2, class Mat3, typename T3>
	LMAT_ENSURE_INLINE
	inline bool have_same_nrows(
			const IMatrixXpr<Mat1, T1>& A,
			const IMatrixXpr<Mat2, T2>& B,
			const IMatrixXpr<Mat3, T3>& C)
	{
		return all_equal(A.nrows(), B.nrows(), C.nrows());
	}

	template<class Mat1, typename T1, class Mat2, typename T2, class Mat3, typename T3, class Mat4, typename T4>
	LMAT_ENSURE_INLINE
	inline bool have_same_nrows(
			const IMatrixXpr<Mat1, T1>& A,
			const IMatrixXpr<Mat2, T2>& B,
			const IMatrixXpr<Mat3, T3>& C,
			const IMatrixXpr<Mat4, T4>& D)
	{
		return all_equal(A.nrows(), B.nrows(), C.nrows(), D.nrows());
	}


	template<class Mat1, typename T1, class Mat2, typename T2>
	LMAT_ENSURE_INLINE
	inline bool have_same_ncols(const IMatrixXpr<Mat1, T1>& A, const IMatrixXpr<Mat2, T2>& B)
	{
		return all_equal(A.ncolumns(), B.ncolumns());
	}

	template<class Mat1, typename T1, class Mat2, typename T2, class Mat3, typename T3>
	LMAT_ENSURE_INLINE
	inline bool have_same_ncols(
			const IMatrixXpr<Mat1, T1>& A,
			const IMatrixXpr<Mat2, T2>& B,
			const IMatrixXpr<Mat3, T3>& C)
	{
		return all_equal(A.ncolumns(), B.ncolumns(), C.ncolumns());
	}

	template<class Mat1, typename T1, class Mat2, typename T2, class Mat3, typename T3, class Mat4, typename T4>
	LMAT_ENSURE_INLINE
	inline bool have_same_ncols(
			const IMatrixXpr<Mat1, T1>& A,
			const IMatrixXpr<Mat2, T2>& B,
			const IMatrixXpr<Mat3, T3>& C,
			const IMatrixXpr<Mat4, T4>& D)
	{
		return all_equal(A.ncolumns(), B.ncolumns(), C.ncolumns(), D.ncolumns());
	}


	template<class Mat1, typename T1, class Mat2, typename T2>
	LMAT_ENSURE_INLINE
	inline bool have_same_shape(const IMatrixXpr<Mat1, T1>& A, const IMatrixXpr<Mat2, T2>& B)
	{
		return have_same_nrows(A, B) && have_same_ncols(A, B);
	}

	template<class Mat1, typename T1, class Mat2, typename T2, class Mat3, typename T3>
	LMAT_ENSURE_INLINE
	inline bool have_same_shape(
			const IMatrixXpr<Mat1, T1>& A,
			const IMatrixXpr<Mat2, T2>& B,
			const IMatrixXpr<Mat3, T3>& C)
	{
		return have_same_nrows(A, B, C) && have_same_ncols(A, B, C);
	}

	template<class Mat1, typename T1, class Mat2, typename T2, class Mat3, typename T3, class Mat4, typename T4>
	LMAT_ENSURE_INLINE
	inline bool have_same_shape(
			const IMatrixXpr<Mat1, T1>& A,
			const IMatrixXpr<Mat2, T2>& B,
			const IMatrixXpr<Mat3, T3>& C,
			const IMatrixXpr<Mat4, T4>& D)
	{
		return have_same_nrows(A, B, C, D) && have_same_ncols(A, B, C, D);
	}


	// common shape

	// common nrows

	template<class Mat1, typename T1, class Mat2, typename T2>
	LMAT_ENSURE_INLINE
	inline index_t common_nrows(const IMatrixXpr<Mat1, T1>& A, const IMatrixXpr<Mat2, T2>& B)
	{
		LMAT_CHECK_DIMS( have_same_nrows(A, B) )
		return A.nrows();
	}

	template<typename T1, class Mat2, typename T2>
	LMAT_ENSURE_INLINE
	inline index_t common_nrows(const scalar_expr<T1>& A, const IMatrixXpr<Mat2, T2>& B)
	{
		return B.nrows();
	}

	template<class Mat1, typename T1, typename T2>
	LMAT_ENSURE_INLINE
	inline index_t common_nrows(const IMatrixXpr<Mat1, T1>& A, const scalar_expr<T2>& B)
	{
		return A.nrows();
	}

	template<class Mat1, typename T1, class Mat2, typename T2, class Mat3, typename T3>
	LMAT_ENSURE_INLINE
	inline index_t common_nrows(
			const IMatrixXpr<Mat1, T1>& A,
			const IMatrixXpr<Mat2, T2>& B,
			const IMatrixXpr<Mat3, T3>& C)
	{
		LMAT_CHECK_DIMS( have_same_nrows(A, B, C) )
		return A.nrows();
	}

	template<typename T1, class Mat2, typename T2, class Mat3, typename T3>
	LMAT_ENSURE_INLINE
	inline index_t common_nrows(
			const scalar_expr<T1>& A,
			const IMatrixXpr<Mat2, T2>& B,
			const IMatrixXpr<Mat3, T3>& C)
	{
		LMAT_CHECK_DIMS( have_same_nrows(B, C) )
		return B.nrows();
	}

	template<class Mat1, typename T1, typename T2, class Mat3, typename T3>
	LMAT_ENSURE_INLINE
	inline index_t common_nrows(
			const IMatrixXpr<Mat1, T1>& A,
			const scalar_expr<T2>& B,
			const IMatrixXpr<Mat3, T3>& C)
	{
		LMAT_CHECK_DIMS( have_same_nrows(A, C) )
		return A.nrows();
	}

	template<class Mat1, typename T1, class Mat2, typename T2, typename T3>
	LMAT_ENSURE_INLINE
	inline index_t common_nrows(
			const IMatrixXpr<Mat1, T1>& A,
			const IMatrixXpr<Mat2, T2>& B,
			const scalar_expr<T3>& C)
	{
		LMAT_CHECK_DIMS( have_same_nrows(A, B) )
		return A.nrows();
	}


	template<typename T1, typename T2, class Mat3, typename T3>
	LMAT_ENSURE_INLINE
	inline index_t common_nrows(
			const scalar_expr<T1>& A,
			const scalar_expr<T2>& B,
			const IMatrixXpr<Mat3, T3>& C)
	{
		return C.nrows();
	}

	template<typename T1, class Mat2, typename T2, typename T3>
	LMAT_ENSURE_INLINE
	inline index_t common_nrows(
			const scalar_expr<T1>& A,
			const IMatrixXpr<Mat2, T2>& B,
			const scalar_expr<T3>& C)
	{
		LMAT_CHECK_DIMS( have_same_nrows(A, B) )
		return B.nrows();
	}

	template<class Mat1, typename T1, typename T2, typename T3>
	LMAT_ENSURE_INLINE
	inline index_t common_nrows(
			const IMatrixXpr<Mat1, T1>& A,
			const scalar_expr<T2>& B,
			const scalar_expr<T3>& C)
	{
		LMAT_CHECK_DIMS( have_same_nrows(A, C) )
		return A.nrows();
	}


	// common ncols

	template<class Mat1, typename T1, class Mat2, typename T2>
	LMAT_ENSURE_INLINE
	inline index_t common_ncols(const IMatrixXpr<Mat1, T1>& A, const IMatrixXpr<Mat2, T2>& B)
	{
		LMAT_CHECK_DIMS( have_same_ncols(A, B) )
		return A.ncolumns();
	}

	template<typename T1, class Mat2, typename T2>
	LMAT_ENSURE_INLINE
	inline index_t common_ncols(const scalar_expr<T1>& A, const IMatrixXpr<Mat2, T2>& B)
	{
		return B.ncolumns();
	}

	template<class Mat1, typename T1, typename T2>
	LMAT_ENSURE_INLINE
	inline index_t common_ncols(const IMatrixXpr<Mat1, T1>& A, const scalar_expr<T2>& B)
	{
		return A.ncolumns();
	}

	template<class Mat1, typename T1, class Mat2, typename T2, class Mat3, typename T3>
	LMAT_ENSURE_INLINE
	inline index_t common_ncols(
			const IMatrixXpr<Mat1, T1>& A,
			const IMatrixXpr<Mat2, T2>& B,
			const IMatrixXpr<Mat3, T3>& C)
	{
		LMAT_CHECK_DIMS( have_same_ncols(A, B, C) )
		return A.ncolumns();
	}

	template<typename T1, class Mat2, typename T2, class Mat3, typename T3>
	LMAT_ENSURE_INLINE
	inline index_t common_ncols(
			const scalar_expr<T1>& A,
			const IMatrixXpr<Mat2, T2>& B,
			const IMatrixXpr<Mat3, T3>& C)
	{
		LMAT_CHECK_DIMS( have_same_ncols(B, C) )
		return B.ncolumns();
	}

	template<class Mat1, typename T1, typename T2, class Mat3, typename T3>
	LMAT_ENSURE_INLINE
	inline index_t common_ncols(
			const IMatrixXpr<Mat1, T1>& A,
			const scalar_expr<T2>& B,
			const IMatrixXpr<Mat3, T3>& C)
	{
		LMAT_CHECK_DIMS( have_same_ncols(A, C) )
		return A.ncolumns();
	}

	template<class Mat1, typename T1, class Mat2, typename T2, typename T3>
	LMAT_ENSURE_INLINE
	inline index_t common_ncols(
			const IMatrixXpr<Mat1, T1>& A,
			const IMatrixXpr<Mat2, T2>& B,
			const scalar_expr<T3>& C)
	{
		LMAT_CHECK_DIMS( have_same_ncols(A, B) )
		return A.ncolumns();
	}


	template<typename T1, typename T2, class Mat3, typename T3>
	LMAT_ENSURE_INLINE
	inline index_t common_ncols(
			const scalar_expr<T1>& A,
			const scalar_expr<T2>& B,
			const IMatrixXpr<Mat3, T3>& C)
	{
		return C.ncolumns();
	}

	template<typename T1, class Mat2, typename T2, typename T3>
	LMAT_ENSURE_INLINE
	inline index_t common_ncols(
			const scalar_expr<T1>& A,
			const IMatrixXpr<Mat2, T2>& B,
			const scalar_expr<T3>& C)
	{
		LMAT_CHECK_DIMS( have_same_ncols(A, B) )
		return B.ncolumns();
	}

	template<class Mat1, typename T1, typename T2, typename T3>
	LMAT_ENSURE_INLINE
	inline index_t common_ncols(
			const IMatrixXpr<Mat1, T1>& A,
			const scalar_expr<T2>& B,
			const scalar_expr<T3>& C)
	{
		LMAT_CHECK_DIMS( have_same_ncols(A, C) )
		return A.ncolumns();
	}

}

#endif 
