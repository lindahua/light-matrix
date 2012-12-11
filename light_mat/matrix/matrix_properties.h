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

	template<typename... Mats>
	LMAT_ENSURE_INLINE
	inline bool have_same_nrows(const Mats&... mats)
	{
		return args_equal(mats.nrows()...);
	}


	template<typename... Mats>
	LMAT_ENSURE_INLINE
	inline bool have_same_ncols(const Mats&... mats)
	{
		return args_equal(mats.ncolumns()...);
	}


	template<typename... Mats>
	LMAT_ENSURE_INLINE
	inline bool have_same_shape(const Mats&... mats)
	{
		return have_same_nrows(mats...) && have_same_ncols(mats...);
	}

	// common shape

	template<typename Mat0, typename... Mats>
	LMAT_ENSURE_INLINE
	inline index_t common_nrows(const Mat0& mat0, const Mats&... mats)
	{
		LMAT_CHECK_DIMS( have_same_nrows(mat0, mats...) );
		return mat0.nrows();
	}

	template<typename Mat0, typename... Mats>
	LMAT_ENSURE_INLINE
	inline index_t common_ncols(const Mat0& mat0, const Mats&... mats)
	{
		LMAT_CHECK_DIMS( have_same_ncols(mat0, mats...) );
		return mat0.ncolumns();
	}

	template<typename Mat0, typename... Mats>
	LMAT_ENSURE_INLINE
	inline matrix_shape<
		meta::common_nrows<Mat0, Mats...>::value,
		meta::common_ncols<Mat0, Mats...>::value>
	common_shape(const Mat0& mat0, const Mats&... mats)
	{
		typedef matrix_shape<
			meta::common_nrows<Mat0, Mats...>::value,
			meta::common_ncols<Mat0, Mats...>::value> shape_t;

		return shape_t(common_nrows(mat0, mats...), common_ncols(mat0, mats...));
	}

}

#endif 
