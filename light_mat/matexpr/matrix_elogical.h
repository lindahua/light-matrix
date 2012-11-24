/*
 * @file matrix_elogical.h
 *
 * Element-wise logical operations
 *
 * @author Dahua Lin
 */

#ifndef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_ELOGICAL_H_
#define LIGHTMAT_MATRIX_ELOGICAL_H_

#include <light_mat/matrix/matrix_ewise_expr.h>
#include <light_mat/math/mask_functors.h>

namespace lmat
{
	// not

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline typename unary_expr_map<
		ewise_t<mask_not_op<T> >,
		ref_arg_t, Mat >::type
	operator ~ (const IMatrixXpr<Mat, mask_t<T> >& A)
	{
		return ewise(mask_not_op<T>(), A.derived());
	}

	// and

	template<typename T, class LMat, class RMat>
	LMAT_ENSURE_INLINE
	inline typename binary_expr_map<
		ewise_t<mask_and_op<T> >,
		ref_arg_t, LMat,
		ref_arg_t, RMat >::type
	operator & (const IMatrixXpr<LMat, mask_t<T> >& A, const IMatrixXpr<RMat, mask_t<T> >& B)
	{
		return ewise(mask_and_op<T>(), A.derived(), B.derived());
	}

	// or

	template<typename T, class LMat, class RMat>
	LMAT_ENSURE_INLINE
	inline typename binary_expr_map<
		ewise_t<mask_or_op<T> >,
		ref_arg_t, LMat,
		ref_arg_t, RMat >::type
	operator | (const IMatrixXpr<LMat, mask_t<T> >& A, const IMatrixXpr<RMat, mask_t<T> >& B)
	{
		return ewise(mask_or_op<T>(), A.derived(), B.derived());
	}

	// xor

	template<typename T, class LMat, class RMat>
	LMAT_ENSURE_INLINE
	inline typename binary_expr_map<
		ewise_t<mask_xor_op<T> >,
		ref_arg_t, LMat,
		ref_arg_t, RMat >::type
	operator ^ (const IMatrixXpr<LMat, mask_t<T> >& A, const IMatrixXpr<RMat, mask_t<T> >& B)
	{
		return ewise(mask_xor_op<T>(), A.derived(), B.derived());
	}

}

#endif 
