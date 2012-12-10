/**
 * @file matrix_fill.h
 *
 * Filling elements to matrices
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_FILL_H_
#define LIGHTMAT_MATRIX_FILL_H_

#include <light_mat/matrix/matrix_properties.h>
#include "internal/matrix_fill_internal.h"

namespace lmat
{
	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline void zero(IRegularMatrix<Mat, T>& dst)
	{
		internal::zero(dst, internal::get_fill_scheme(dst));
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline void fill(IRegularMatrix<Mat, T>& dst, const T& val)
	{
		internal::fill(val, dst, internal::get_fill_scheme(dst));
	}

	template<typename T, class DMat>
	LMAT_ENSURE_INLINE
	inline DMat& operator << (IRegularMatrix<DMat, T>& dmat, const T& v)
	{
		fill(dmat, v);
		return dmat.derived();
	}

}

#endif /* MATRIX_FILL_H_ */
