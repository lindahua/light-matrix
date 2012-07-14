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
	void zero(IDenseMatrix<Mat, T>& dst)
	{
		typedef typename detail::mat_filler<T,
				ct_rows<Mat>::value,
				ct_cols<Mat>::value>::type filler_t;

		if (has_continuous_layout(dst))
		{
			filler_t::zero(dst.nrows(), dst.ncolumns(), dst.ptr_data());
		}
		else
		{
			filler_t::zero(dst.nrows(), dst.ncolumns(),
					dst.ptr_data(), dst.lead_dim());
		}
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	void fill(IDenseMatrix<Mat, T>& dst, const T& val)
	{
		typedef typename detail::mat_filler<T,
				ct_rows<Mat>::value,
				ct_cols<Mat>::value>::type filler_t;

		if (has_continuous_layout(dst))
		{
			filler_t::fill(val, dst.nrows(), dst.ncolumns(), dst.ptr_data());
		}
		else
		{
			filler_t::fill(val, dst.nrows(), dst.ncolumns(),
					dst.ptr_data(), dst.lead_dim());
		}
	}

}

#endif /* MATRIX_FILL_H_ */
