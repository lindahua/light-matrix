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

	template<class Mat>
	LMAT_ENSURE_INLINE
	inline typename unary_ewise_expr_map<not_op, Mat>::type
	operator ! (const IMatrixXpr<Mat, bool>& A)
	{
		return ewise(not_op(), A.derived());
	}

	// and

	template<class LMat, class RMat>
	LMAT_ENSURE_INLINE
	inline typename binary_ewise_expr_map<and_op, LMat, RMat>::type
	operator && (const IMatrixXpr<LMat, bool>& A, const IMatrixXpr<RMat, bool>& B)
	{
		return ewise(and_op(), A.derived(), B.derived());
	}

	// or

	template<class LMat, class RMat>
	LMAT_ENSURE_INLINE
	inline typename binary_ewise_expr_map<or_op, LMat, RMat>::type
	operator || (const IMatrixXpr<LMat, bool>& A, const IMatrixXpr<RMat, bool>& B)
	{
		return ewise(or_op(), A.derived(), B.derived());
	}

}

#endif 
