/**
 * @file matrix_asvec.h
 *
 * @brief Matrix as-vector views
 *
 * @author Dahua Lin
 */

#ifndef LIGHTMAT_MATRIX_ASVEC_H_
#define LIGHTMAT_MATRIX_ASVEC_H_

#include <light_mat/matrix/matrix_concepts.h>
#include "internal/matrix_asvec_internal.h"

namespace lmat
{
	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline typename internal::as_col_map<Mat>::const_type
	as_col(const IRegularMatrix<Mat, T>& mat)
	{
		return internal::as_col_map<Mat>::get(mat.derived());
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline typename internal::as_col_map<Mat>::_type
	as_col(IRegularMatrix<Mat, T>& mat)
	{
		return internal::as_col_map<Mat>::get(mat.derived());
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline typename internal::as_row_map<Mat>::const_type
	as_row(const IRegularMatrix<Mat, T>& mat)
	{
		return internal::as_row_map<Mat>::get(mat.derived());
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline typename internal::as_row_map<Mat>::_type
	as_row(IRegularMatrix<Mat, T>& mat)
	{
		return internal::as_row_map<Mat>::get(mat.derived());
	}

}

#endif /* MATRIX_ASVEC_H_ */
