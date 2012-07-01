/*
 * @file matrix_ecomp.h
 *
 * Element-wise comparison operations
 *
 * @author Dahua Lin
 */

#ifndef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_ECOMP_H_
#define LIGHTMAT_MATRIX_ECOMP_H_

#include <light_mat/matrix/matrix_ewise_eval.h>
#include <light_mat/math/comp_functors.h>

namespace lmat
{
	// eq

	template<typename T, class LMat, class RMat>
	LMAT_ENSURE_INLINE
	inline binary_ewise_expr<eq_op<T>, LMat, RMat>
	operator == (const IMatrixXpr<LMat, T>& A, const IMatrixXpr<RMat, T>& B)
	{
		return ewise_expr(eq_op<T>(), A, B);
	}

	template<typename T, class LMat>
	LMAT_ENSURE_INLINE
	inline binary_fix2_ewise_expr<eq_op<T>, LMat>
	operator ==(const IMatrixXpr<LMat, T>& A, const T& b)
	{
		typedef typename const_mat_same_form<LMat>::type cst_t;
		return ewise_expr(eq_op<T>(), A, b);
	}

	template<typename T, class RMat>
	LMAT_ENSURE_INLINE
	inline binary_fix1_ewise_expr<eq_op<T>, RMat>
	operator == (const T& a, const IMatrixXpr<RMat, T>& B)
	{
		typedef typename const_mat_same_form<RMat>::type cst_t;
		return ewise_expr(eq_op<T>(), a, B);
	}

	// ne

	template<typename T, class LMat, class RMat>
	LMAT_ENSURE_INLINE
	inline binary_ewise_expr<ne_op<T>, LMat, RMat>
	operator != (const IMatrixXpr<LMat, T>& A, const IMatrixXpr<RMat, T>& B)
	{
		return ewise_expr(ne_op<T>(), A, B);
	}

	template<typename T, class LMat>
	LMAT_ENSURE_INLINE
	inline binary_fix2_ewise_expr<ne_op<T>, LMat>
	operator !=(const IMatrixXpr<LMat, T>& A, const T& b)
	{
		typedef typename const_mat_same_form<LMat>::type cst_t;
		return ewise_expr(ne_op<T>(), A, b);
	}

	template<typename T, class RMat>
	LMAT_ENSURE_INLINE
	inline binary_fix1_ewise_expr<ne_op<T>, RMat>
	operator != (const T& a, const IMatrixXpr<RMat, T>& B)
	{
		typedef typename const_mat_same_form<RMat>::type cst_t;
		return ewise_expr(ne_op<T>(), a, B);
	}


	// le

	template<typename T, class LMat, class RMat>
	LMAT_ENSURE_INLINE
	inline binary_ewise_expr<le_op<T>, LMat, RMat>
	operator <= (const IMatrixXpr<LMat, T>& A, const IMatrixXpr<RMat, T>& B)
	{
		return ewise_expr(le_op<T>(), A, B);
	}

	template<typename T, class LMat>
	LMAT_ENSURE_INLINE
	inline binary_fix2_ewise_expr<le_op<T>, LMat>
	operator <=(const IMatrixXpr<LMat, T>& A, const T& b)
	{
		typedef typename const_mat_same_form<LMat>::type cst_t;
		return ewise_expr(le_op<T>(), A, b);
	}

	template<typename T, class RMat>
	LMAT_ENSURE_INLINE
	inline binary_fix1_ewise_expr<le_op<T>, RMat>
	operator <= (const T& a, const IMatrixXpr<RMat, T>& B)
	{
		typedef typename const_mat_same_form<RMat>::type cst_t;
		return ewise_expr(le_op<T>(), a, B);
	}

	// lt

	template<typename T, class LMat, class RMat>
	LMAT_ENSURE_INLINE
	inline binary_ewise_expr<lt_op<T>, LMat, RMat>
	operator < (const IMatrixXpr<LMat, T>& A, const IMatrixXpr<RMat, T>& B)
	{
		return ewise_expr(lt_op<T>(), A, B);
	}

	template<typename T, class LMat>
	LMAT_ENSURE_INLINE
	inline binary_fix2_ewise_expr<lt_op<T>, LMat>
	operator < (const IMatrixXpr<LMat, T>& A, const T& b)
	{
		typedef typename const_mat_same_form<LMat>::type cst_t;
		return ewise_expr(lt_op<T>(), A, b);
	}

	template<typename T, class RMat>
	LMAT_ENSURE_INLINE
	inline binary_fix1_ewise_expr<lt_op<T>, RMat>
	operator < (const T& a, const IMatrixXpr<RMat, T>& B)
	{
		typedef typename const_mat_same_form<RMat>::type cst_t;
		return ewise_expr(lt_op<T>(), a, B);
	}



	// ge

	template<typename T, class LMat, class RMat>
	LMAT_ENSURE_INLINE
	inline binary_ewise_expr<ge_op<T>, LMat, RMat>
	operator >= (const IMatrixXpr<LMat, T>& A, const IMatrixXpr<RMat, T>& B)
	{
		return ewise_expr(ge_op<T>(), A, B);
	}

	template<typename T, class LMat>
	LMAT_ENSURE_INLINE
	inline binary_fix2_ewise_expr<ge_op<T>, LMat>
	operator >=(const IMatrixXpr<LMat, T>& A, const T& b)
	{
		typedef typename const_mat_same_form<LMat>::type cst_t;
		return ewise_expr(ge_op<T>(), A, b);
	}

	template<typename T, class RMat>
	LMAT_ENSURE_INLINE
	inline binary_fix1_ewise_expr<ge_op<T>, RMat>
	operator >= (const T& a, const IMatrixXpr<RMat, T>& B)
	{
		typedef typename const_mat_same_form<RMat>::type cst_t;
		return ewise_expr(ge_op<T>(), a, B);
	}

	// gt

	template<typename T, class LMat, class RMat>
	LMAT_ENSURE_INLINE
	inline binary_ewise_expr<gt_op<T>, LMat, RMat>
	operator > (const IMatrixXpr<LMat, T>& A, const IMatrixXpr<RMat, T>& B)
	{
		return ewise_expr(gt_op<T>(), A, B);
	}

	template<typename T, class LMat>
	LMAT_ENSURE_INLINE
	inline binary_fix2_ewise_expr<gt_op<T>, LMat>
	operator >(const IMatrixXpr<LMat, T>& A, const T& b)
	{
		typedef typename const_mat_same_form<LMat>::type cst_t;
		return ewise_expr(gt_op<T>(), A, b);
	}

	template<typename T, class RMat>
	LMAT_ENSURE_INLINE
	inline binary_fix1_ewise_expr<gt_op<T>, RMat>
	operator > (const T& a, const IMatrixXpr<RMat, T>& B)
	{
		typedef typename const_mat_same_form<RMat>::type cst_t;
		return ewise_expr(gt_op<T>(), a, B);
	}
}

#endif 
