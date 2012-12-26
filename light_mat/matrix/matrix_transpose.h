/**
 * @file matrix_transpose.h
 *
 * @brief Matrix transposition
 *
 * @author Dahua Lin
 */

#ifndef LIGHTMAT_MATRIX_TRANSPOSE_H_
#define LIGHTMAT_MATRIX_TRANSPOSE_H_

#include "internal/matrix_transpose_internal.h"

namespace lmat
{
	template<typename T, class SMat, class DMat>
	LMAT_ENSURE_INLINE
	inline void transpose(const IRegularMatrix<SMat, T>& smat, IRegularMatrix<DMat, T>& dmat)
	{
		index_t m = smat.nrows();
		index_t n = smat.ncolumns();

		LMAT_CHECK_DIMS( dmat.nrows() == n && dmat.ncolumns() == m );

		internal::direct_transpose(m, n, smat, dmat);
	}

}

#endif /* MAT_TRANSPOSE_H_ */
