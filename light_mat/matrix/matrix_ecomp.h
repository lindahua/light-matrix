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
	inline typename binary_ewise_expr_map<eq_op<T>, LMat, RMat>::type
	operator == (const IMatrixXpr<LMat, T>& A, const IMatrixXpr<RMat, T>& B)
	{
		return ewise(eq_op<T>(), A.derived(), B.derived());
	}

	template<typename T, class LMat>
	LMAT_ENSURE_INLINE
	inline typename binary_fix2_ewise_expr_map<eq_op<T>, LMat>::type
	operator ==(const IMatrixXpr<LMat, T>& A, const T& b)
	{
		typedef typename const_mat_same_form<LMat>::type cst_t;
		return ewise(eq_op<T>(), A.derived(), b);
	}

	template<typename T, class RMat>
	LMAT_ENSURE_INLINE
	inline typename binary_fix1_ewise_expr_map<eq_op<T>, RMat>::type
	operator == (const T& a, const IMatrixXpr<RMat, T>& B)
	{
		typedef typename const_mat_same_form<RMat>::type cst_t;
		return ewise(eq_op<T>(), a, B.derived());
	}

	// ne

	template<typename T, class LMat, class RMat>
	LMAT_ENSURE_INLINE
	inline typename binary_ewise_expr_map<ne_op<T>, LMat, RMat>::type
	operator != (const IMatrixXpr<LMat, T>& A, const IMatrixXpr<RMat, T>& B)
	{
		return ewise(ne_op<T>(), A.derived(), B.derived());
	}

	template<typename T, class LMat>
	LMAT_ENSURE_INLINE
	inline typename binary_fix2_ewise_expr_map<ne_op<T>, LMat>::type
	operator !=(const IMatrixXpr<LMat, T>& A, const T& b)
	{
		typedef typename const_mat_same_form<LMat>::type cst_t;
		return ewise(ne_op<T>(), A.derived(), b);
	}

	template<typename T, class RMat>
	LMAT_ENSURE_INLINE
	inline typename binary_fix1_ewise_expr_map<ne_op<T>, RMat>::type
	operator != (const T& a, const IMatrixXpr<RMat, T>& B)
	{
		typedef typename const_mat_same_form<RMat>::type cst_t;
		return ewise(ne_op<T>(), a, B.derived());
	}


	// le

	template<typename T, class LMat, class RMat>
	LMAT_ENSURE_INLINE
	inline typename binary_ewise_expr_map<le_op<T>, LMat, RMat>::type
	operator <= (const IMatrixXpr<LMat, T>& A, const IMatrixXpr<RMat, T>& B)
	{
		return ewise(le_op<T>(), A.derived(), B.derived());
	}

	template<typename T, class LMat>
	LMAT_ENSURE_INLINE
	inline typename binary_fix2_ewise_expr_map<le_op<T>, LMat>::type
	operator <=(const IMatrixXpr<LMat, T>& A, const T& b)
	{
		typedef typename const_mat_same_form<LMat>::type cst_t;
		return ewise(le_op<T>(), A.derived(), b);
	}

	template<typename T, class RMat>
	LMAT_ENSURE_INLINE
	inline typename binary_fix1_ewise_expr_map<le_op<T>, RMat>::type
	operator <= (const T& a, const IMatrixXpr<RMat, T>& B)
	{
		typedef typename const_mat_same_form<RMat>::type cst_t;
		return ewise(le_op<T>(), a, B.derived());
	}

	// lt

	template<typename T, class LMat, class RMat>
	LMAT_ENSURE_INLINE
	inline typename binary_ewise_expr_map<lt_op<T>, LMat, RMat>::type
	operator < (const IMatrixXpr<LMat, T>& A, const IMatrixXpr<RMat, T>& B)
	{
		return ewise(lt_op<T>(), A.derived(), B.derived());
	}

	template<typename T, class LMat>
	LMAT_ENSURE_INLINE
	inline typename binary_fix2_ewise_expr_map<lt_op<T>, LMat>::type
	operator < (const IMatrixXpr<LMat, T>& A, const T& b)
	{
		typedef typename const_mat_same_form<LMat>::type cst_t;
		return ewise(lt_op<T>(), A.derived(), b);
	}

	template<typename T, class RMat>
	LMAT_ENSURE_INLINE
	inline typename binary_fix1_ewise_expr_map<lt_op<T>, RMat>::type
	operator < (const T& a, const IMatrixXpr<RMat, T>& B)
	{
		typedef typename const_mat_same_form<RMat>::type cst_t;
		return ewise(lt_op<T>(), a, B.derived());
	}



	// ge

	template<typename T, class LMat, class RMat>
	LMAT_ENSURE_INLINE
	inline typename binary_ewise_expr_map<ge_op<T>, LMat, RMat>::type
	operator >= (const IMatrixXpr<LMat, T>& A, const IMatrixXpr<RMat, T>& B)
	{
		return ewise(ge_op<T>(), A.derived(), B.derived());
	}

	template<typename T, class LMat>
	LMAT_ENSURE_INLINE
	inline typename binary_fix2_ewise_expr_map<ge_op<T>, LMat>::type
	operator >=(const IMatrixXpr<LMat, T>& A, const T& b)
	{
		typedef typename const_mat_same_form<LMat>::type cst_t;
		return ewise(ge_op<T>(), A.derived(), b);
	}

	template<typename T, class RMat>
	LMAT_ENSURE_INLINE
	inline typename binary_fix1_ewise_expr_map<ge_op<T>, RMat>::type
	operator >= (const T& a, const IMatrixXpr<RMat, T>& B)
	{
		typedef typename const_mat_same_form<RMat>::type cst_t;
		return ewise(ge_op<T>(), a, B.derived());
	}

	// gt

	template<typename T, class LMat, class RMat>
	LMAT_ENSURE_INLINE
	inline typename binary_ewise_expr_map<gt_op<T>, LMat, RMat>::type
	operator > (const IMatrixXpr<LMat, T>& A, const IMatrixXpr<RMat, T>& B)
	{
		return ewise(gt_op<T>(), A.derived(), B.derived());
	}

	template<typename T, class LMat>
	LMAT_ENSURE_INLINE
	inline typename binary_fix2_ewise_expr_map<gt_op<T>, LMat>::type
	operator >(const IMatrixXpr<LMat, T>& A, const T& b)
	{
		typedef typename const_mat_same_form<LMat>::type cst_t;
		return ewise(gt_op<T>(), A.derived(), b);
	}

	template<typename T, class RMat>
	LMAT_ENSURE_INLINE
	inline typename binary_fix1_ewise_expr_map<gt_op<T>, RMat>::type
	operator > (const T& a, const IMatrixXpr<RMat, T>& B)
	{
		typedef typename const_mat_same_form<RMat>::type cst_t;
		return ewise(gt_op<T>(), a, B.derived());
	}
}

#endif 
