/**
 * @file matrix_blasl1.h
 *
 * BLAS Level 1 on matrices
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_BLASL1_H_
#define LIGHTMAT_MATRIX_BLASL1_H_

#include "bits/matrix_blasl1_internal.h"

namespace lmat { namespace blas {

	template<class Mat>
	LMAT_ENSURE_INLINE
	float asum(const IMatrixXpr<Mat, float>& x)
	{
		typedef typename lmat::detail::unary_blasl1_internal_map<Mat>::type intern_t;
		intern_t::asum(x);
	}

	template<class Mat>
	LMAT_ENSURE_INLINE
	double asum(const IMatrixXpr<Mat, double>& x)
	{
		typedef typename lmat::detail::unary_blasl1_internal_map<Mat>::type intern_t;
		intern_t::asum(x);
	}



} }

#endif /* MATRIX_BLASL1_H_ */
