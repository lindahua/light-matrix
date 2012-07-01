/*
 * @file matrix_elogical.h
 *
 * Element-wise boolean operations
 *
 * @author Dahua Lin
 */

#ifndef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_BOOLOP_H_
#define LIGHTMAT_MATRIX_BOOLOP_H_

#include <light_mat/matrix/matrix_ewise_eval.h>
#include <light_mat/math/bool_functors.h>

namespace lmat
{
	// not

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline unary_ewise_expr<not_op, Mat>
	operator ! (const IMatrixXpr<Mat, bool>& A)
	{
		return ewise_expr(and_op(), A);
	}


	// and

	template<typename T, class LMat, class RMat>
	LMAT_ENSURE_INLINE
	inline binary_ewise_expr<and_op, LMat, RMat>
	operator && (const IMatrixXpr<LMat, bool>& A, const IMatrixXpr<RMat, bool>& B)
	{
		return ewise_expr(and_op(), A, B);
	}

	template<typename T, class LMat>
	LMAT_ENSURE_INLINE
	inline binary_fix2_ewise_expr<and_op, LMat>
	operator && (const IMatrixXpr<LMat, bool>& A, bool b)
	{
		return ewise_expr(and_op(), A, b);
	}

	template<class RMat>
	LMAT_ENSURE_INLINE
	inline binary_fix1_ewise_expr<and_op, RMat>
	operator && (bool a, const IMatrixXpr<RMat, bool>& B)
	{
		return ewise_expr(and_op(), a, B);
	}


	// or

	template<typename T, class LMat, class RMat>
	LMAT_ENSURE_INLINE
	inline binary_ewise_expr<or_op, LMat, RMat>
	operator || (const IMatrixXpr<LMat, bool>& A, const IMatrixXpr<RMat, bool>& B)
	{
		return ewise_expr(or_op(), A, B);
	}

	template<typename T, class LMat>
	LMAT_ENSURE_INLINE
	inline binary_fix2_ewise_expr<or_op, LMat>
	operator || (const IMatrixXpr<LMat, bool>& A, bool b)
	{
		return ewise_expr(or_op(), A, b);
	}

	template<class RMat>
	LMAT_ENSURE_INLINE
	inline binary_fix1_ewise_expr<or_op, RMat>
	operator || (bool a, const IMatrixXpr<RMat, bool>& B)
	{
		return ewise_expr(or_op(), a, B);
	}

}

#endif 
