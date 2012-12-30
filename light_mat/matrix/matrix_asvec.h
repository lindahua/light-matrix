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
#include <vector>

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


	template<typename T, typename Alloc>
	LMAT_ENSURE_INLINE
	inline cref_col<T> as_col(const std::vector<T, Alloc>& vec)
	{
		const T* p = vec.data();
		return cref_col<T>(p, static_cast<index_t>(vec.size()));
	}

	template<typename T, typename Alloc>
	LMAT_ENSURE_INLINE
	inline ref_col<T> as_col(std::vector<T, Alloc>& vec)
	{
		T* p = vec.data();
		return ref_col<T>(p, static_cast<index_t>(vec.size()));
	}

	template<typename T, typename Alloc>
	LMAT_ENSURE_INLINE
	inline cref_row<T> as_row(const std::vector<T, Alloc>& vec)
	{
		const T* p = vec.data();
		return cref_row<T>(p, static_cast<index_t>(vec.size()));
	}

	template<typename T, typename Alloc>
	LMAT_ENSURE_INLINE
	inline ref_row<T> as_row(std::vector<T, Alloc>& vec)
	{
		T* p = vec.data();
		return ref_row<T>(p, static_cast<index_t>(vec.size()));
	}

}

#endif /* MATRIX_ASVEC_H_ */
