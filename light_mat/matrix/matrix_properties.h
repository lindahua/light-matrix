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
#define LMAT_CHECK_SAME_DIM(n1, n2) check_arg(n1 == n2, "Inconsistent matrix dimensions");
#define LMAT_CHECK_SAME_SHAPE(a, b) check_arg(has_same_size((a), (b)), "Inconsistent matrix shape.");
#else
#define LMAT_CHECK_SAME_DIM(n1, n2)
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


	// common shape

	namespace internal
	{
		template<int D>
		struct get_common_dim_helper
		{
			LMAT_ENSURE_INLINE
			static index_t get(index_t, index_t) { return D; }

			LMAT_ENSURE_INLINE
			static index_t get(index_t, index_t, index_t) { return D; }
		};

		template<>
		struct get_common_dim_helper<0>
		{
			LMAT_ENSURE_INLINE
			static index_t get(index_t d1, index_t d2)
			{
				LMAT_CHECK_SAME_DIM(d1, d2)
				return d1;
			}

			LMAT_ENSURE_INLINE
			static index_t get(index_t d1, index_t d2, index_t d3)
			{
				LMAT_CHECK_SAME_DIM(d1, d2)
				LMAT_CHECK_SAME_DIM(d2, d3)
				return d1;
			}
		};
	}


	template<int D1, int D2>
	LMAT_ENSURE_INLINE
	inline index_t get_common_dim(const index_t& d1, const index_t& d2)
	{
		return internal::get_common_dim_helper<common_ctdim_c<D1, D2>::value>::get(d1, d2);
	}

	template<int D1, int D2, int D3>
	LMAT_ENSURE_INLINE
	inline index_t get_common_dim(const index_t& d1, const index_t& d2, const index_t& d3)
	{
		return internal::get_common_dim_helper<common_ctdim_c<D1, D2, D3>::value>::get(d1, d2, d3);
	}


	// common nrows

	template<class Mat1, typename T1, class Mat2, typename T2>
	LMAT_ENSURE_INLINE
	inline index_t get_common_nrows(
			const IMatrixXpr<Mat1, T1>& A,
			const IMatrixXpr<Mat2, T2>& B)
	{
		return get_common_dim<ct_rows<Mat1>::value, ct_rows<Mat2>::value>(A.nrows(), B.nrows());
	}

	template<typename T1, class Mat2, typename T2>
	LMAT_ENSURE_INLINE
	inline index_t get_common_nrows(
			const scalar_expr<T1>& A,
			const IMatrixXpr<Mat2, T2>& B)
	{
		return B.nrows();
	}

	template<typename T1, class Mat1, typename T2>
	LMAT_ENSURE_INLINE
	inline index_t get_common_nrows(
			const IMatrixXpr<Mat1, T1>& A,
			const scalar_expr<T2>& B)
	{
		return B.nrows();
	}

	template<class Mat1, typename T1, class Mat2, typename T2, class Mat3, typename T3>
	LMAT_ENSURE_INLINE
	inline index_t get_common_nrows(
			const IMatrixXpr<Mat1, T1>& A,
			const IMatrixXpr<Mat2, T2>& B,
			const IMatrixXpr<Mat3, T3>& C)
	{
		return get_common_dim<ct_rows<Mat1>::value, ct_rows<Mat2>::value, ct_rows<Mat3>::value>(
				A.nrows(), B.nrows(), C.nrows());
	}


	// common ncolumns

	template<class Mat1, typename T1, class Mat2, typename T2>
	LMAT_ENSURE_INLINE
	inline index_t get_common_ncolumns(
			const IMatrixXpr<Mat1, T1>& A,
			const IMatrixXpr<Mat2, T2>& B)
	{
		return get_common_dim<ct_cols<Mat1>::value, ct_cols<Mat2>::value>(A.ncolumns(), B.ncolumns());
	}

	template<typename T1, class Mat2, typename T2>
	LMAT_ENSURE_INLINE
	inline index_t get_common_ncolumns(
			const scalar_expr<T1>& A,
			const IMatrixXpr<Mat2, T2>& B)
	{
		return B.ncolumns();
	}

	template<typename T1, class Mat1, typename T2>
	LMAT_ENSURE_INLINE
	inline index_t get_common_ncolumns(
			const IMatrixXpr<Mat1, T1>& A,
			const scalar_expr<T2>& B)
	{
		return B.ncolumns();
	}

	template<class Mat1, typename T1, class Mat2, typename T2, class Mat3, typename T3>
	LMAT_ENSURE_INLINE
	inline index_t get_common_ncolumns(
			const IMatrixXpr<Mat1, T1>& A,
			const IMatrixXpr<Mat2, T2>& B,
			const IMatrixXpr<Mat3, T3>& C)
	{
		return get_common_dim<ct_cols<Mat1>::value, ct_cols<Mat2>::value, ct_cols<Mat3>::value>(
				A.ncolumns(), B.ncolumns(), C.ncolumns());
	}

}

#endif 
