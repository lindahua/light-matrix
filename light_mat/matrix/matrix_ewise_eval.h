/**
 * @file matrix_ewise_eval.h
 *
 * Evaluation of element-wise matrix expressions
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_EWISE_EVAL_H_
#define LIGHTMAT_MATRIX_EWISE_EVAL_H_

#include <light_mat/matrix/matrix_ewise_expr.h>
#include <light_mat/matrix/vec_evaluator_concepts.h>
#include <light_mat/matrix/const_matrix.h>

#include <light_mat/math/arith_functors.h>

namespace lmat
{
	template<typename T, class LMat, class RMat>
	LMAT_ENSURE_INLINE
	simple_binary_ewise_expr<T, LMat, RMat>
	operator + (const IMatrixXpr<LMat, T>& A, const IMatrixXpr<RMat, T>& B)
	{
		return simple_ewise_expr(add_op<T>(), A, B);
	}

	template<typename T, class LMat>
	LMAT_ENSURE_INLINE
	simple_binary_ewise_expr<T, LMat, typename const_mat_same_form<LMat>::type>
	operator + (const IMatrixXpr<LMat, T>& A, const T& b)
	{
		return simple_ewise_expr(add_op<T>(), A, const_mat_same_form<LMat>::get(b));
	}

	template<typename T, class RMat>
	LMAT_ENSURE_INLINE
	simple_binary_ewise_expr<T, typename const_mat_same_form<RMat>::type, RMat>
	operator + (const T& a, const IMatrixXpr<RMat, T>& B)
	{
		return simple_ewise_expr(add_op<T>(), const_mat_same_form<RMat>::get(a), B);
	}

}

#endif /* MATRIX_EWISE_EVAL_H_ */
