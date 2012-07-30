/**
 * @file matrix_blasl2.h
 *
 * BLAS Level 2 on matrices
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_BLASL2_H_
#define LIGHTMAT_MATRIX_BLASL2_H_

#include "bits/matrix_blasl2_internal.h"

namespace lmat
{
	template<class MatA, class VecX, class VecY>
	LMAT_ENSURE_INLINE
	void gemv(
			const IMatrixXpr<MatA, double>& A,
			const IMatrixXpr<VecX, double>& x,
			IDenseMatrix<VecY, double>& y)
	{
		typedef typename lmat::detail::gemv_n_internal_map<MatA, VecX, VecY>::type intern_t;
		intern_t::eval_n(A, x, y);
	}



}

#endif /* MATRIX_BLASL2_H_ */
