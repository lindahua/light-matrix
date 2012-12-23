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

#include <light_mat/linalg/linalg_fwd.h>

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


	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline index_t scale_columns(index_t m, index_t n, index_t lda, double *A, double *s)
	{
		for (index_t j = 0; j < n; ++j)
		{
			double *a = A + j * lda;
			for (index_t i = 0; i < m; ++i) a[i] *= s[i];
		}
	}

} }

#endif 
