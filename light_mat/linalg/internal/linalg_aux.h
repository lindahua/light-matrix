/**
 * @file linalg_aux.h
 *
 * Auxiliary function to support linear algebra routines
 * 
 * @author Dahua Lin 
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_LINALG_AUX_H_
#define LIGHTMAT_LINALG_AUX_H_

#include <light_mat/matrix/matrix_properties.h>

namespace lmat { namespace internal {

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline index_t get_vector_intv(const IRegularMatrix<Mat, T>& mat)
	{
		index_t intv;

		if (meta::is_continuous<Mat>::value)
		{
			intv = 1;
		}
		else if (is_column(mat))
		{
			intv = mat.row_stride();
		}
		else if (is_row(mat))
		{
			intv = mat.col_stride();
		}
		else
		{
			throw invalid_argument("The input is not a proper vector.");
		}

		return intv;
	}

} }

#endif 
