/**
 * @file mat_compare.h
 *
 * @brief Matrix comparison
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MAT_COMPARE_H_
#define LIGHTMAT_MAT_COMPARE_H_

#include <light_mat/mateval/mat_allany.h>
#include <light_mat/matexpr/mat_arith.h>

namespace lmat
{

	template<typename T, class A, class B>
	inline bool is_equal(const IEWiseMatrix<A, T>& a, const IEWiseMatrix<B, T>& b)
	{
		return have_same_shape(a, b) && all(a == b);
	}

	template<typename T, class A, class B>
	inline bool is_approx(const IEWiseMatrix<A, T>& a, const IEWiseMatrix<B, T>& b, const T& tol)
	{
		return have_same_shape(a, b) && all(abs(a - b) < tol);
	}

}


#endif /* MAT_COMPARE_H_ */
