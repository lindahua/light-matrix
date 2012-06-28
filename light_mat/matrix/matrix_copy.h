/*
 * @file matrix_copy.h
 *
 * Functions to copy matrices
 *
 * @author Dahua Lin
 */

#ifndef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_COPY_H_
#define LIGHTMAT_MATRIX_COPY_H_

#include <light_mat/matrix/matrix_concepts.h>
#include "bits/matrix_copy_internal.h"

namespace lmat
{
	template<typename T, class LMat, class RMat>
	LMAT_ENSURE_INLINE
	void copy(const IDenseMatrix<LMat, T>& src, IDenseMatrix<RMat, T>& dst)
	{
		check_same_size(src, dst, "copy: inconsistent sizes of src and dst.");

		const int M = binary_ct_rows<LMat, RMat>::value;
		const int N = binary_ct_cols<LMat, RMat>::value;
		typedef detail::mat_copier<T, M, N> copier_t;

		copier_t::_copy2d(src.nrows(), src.ncolumns(),
				src.ptr_data(), src.lead_dim(),
				dst.ptr_data(), dst.lead_dim());
	}
}

#endif 
