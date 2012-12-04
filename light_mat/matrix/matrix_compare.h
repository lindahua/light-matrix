/**
 * @file matrix_compare.h
 *
 * Functions for matrix comparison
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_COMPARE_H_
#define LIGHTMAT_MATRIX_COMPARE_H_

#include <light_mat/matrix/matrix_properties.h>
#include "bits/matrix_compare_internal.h"
#include <cmath>

namespace lmat
{
	template<typename T, class LMat, class RMat>
	LMAT_ENSURE_INLINE
	inline bool is_equal(const IDenseMatrix<LMat, T>& a, const IDenseMatrix<RMat, T>& b)
	{
		typedef typename internal::mat_comparer_map<LMat, RMat>::type comparer_t;

		return have_same_shape(a, b) && comparer_t::all_equal(a.derived(), b.derived());
	}

	template<typename T, class LMat, class RMat>
	LMAT_ENSURE_INLINE
	inline bool is_approx(const IDenseMatrix<LMat, T>& a, const IDenseMatrix<RMat, T>& b, const T& tol)
	{
		typedef typename internal::mat_approx_comparer_map<LMat, RMat>::type comparer_t;

		return have_same_shape(a, b) && comparer_t::all_approx(a.derived(), b.derived(), tol);
	}

}

#endif
