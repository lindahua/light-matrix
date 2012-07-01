/**
 * @file matrix_arith.h
 *
 * Element-wise arithmetics on matrices
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_ARITH_H_
#define LIGHTMAT_MATRIX_ARITH_H_

#include <light_mat/matrix/matrix_ewise_eval.h>
#include <light_mat/math/arith_functors.h>

namespace lmat
{
	template<typename T, class LMat, class RMat>
	LMAT_ENSURE_INLINE
	inline binary_ewise_expr<add_op<T>, LMat, RMat>
	operator + (const IMatrixXpr<LMat, T>& A, const IMatrixXpr<RMat, T>& B)
	{
		return ewise_expr(add_op<T>(), A, B);
	}

	template<typename T, class LMat>
	LMAT_ENSURE_INLINE
	inline binary_fix2_ewise_expr<add_op<T>, LMat>
	operator + (const IMatrixXpr<LMat, T>& A, const T& b)
	{
		typedef typename const_mat_same_form<LMat>::type cst_t;
		return ewise_expr(add_op<T>(), A, b);
	}

	template<typename T, class RMat>
	LMAT_ENSURE_INLINE
	inline binary_fix1_ewise_expr<add_op<T>, RMat>
	operator + (const T& a, const IMatrixXpr<RMat, T>& B)
	{
		typedef typename const_mat_same_form<RMat>::type cst_t;
		return ewise_expr(add_op<T>(), a, B);
	}

}



#endif /* MATRIX_ARITH_H_ */
