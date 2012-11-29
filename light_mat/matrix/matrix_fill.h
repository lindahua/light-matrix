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
#include "bits/matrix_fill_internal.h"

namespace lmat
{
	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline void zero(IDenseMatrix<Mat, T>& dst)
	{
		typedef typename internal::mat_filler_map<Mat>::type filler_t;
		filler_t::zero(dst.derived());
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline void fill(IDenseMatrix<Mat, T>& dst, const T& val)
	{
		typedef typename internal::mat_filler_map<Mat>::type filler_t;
		filler_t::fill(val, dst.derived());
	}

	template<typename T, class DMat>
	LMAT_ENSURE_INLINE
	inline DMat& operator << (IDenseMatrix<DMat, T>& dmat, const T& v)
	{
		fill(dmat.derived(), v);
		return dmat.derived();
	}

}

#endif /* MATRIX_FILL_H_ */
