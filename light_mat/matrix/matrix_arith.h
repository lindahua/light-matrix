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
	/********************************************
	 *
	 *  Add, Sub, Mul, and Div
	 *
	 ********************************************/

	// Addition

	template<typename T, class LMat, class RMat>
	LMAT_ENSURE_INLINE
	inline typename binary_ewise_expr_map<add_op<T>, LMat, RMat>::type
	operator + (const IMatrixXpr<LMat, T>& A, const IMatrixXpr<RMat, T>& B)
	{
		return ewise(add_op<T>(), A.derived(), B.derived());
	}

	template<typename T, class LMat>
	LMAT_ENSURE_INLINE
	inline typename binary_fix2_ewise_expr_map<add_op<T>, LMat>::type
	operator + (const IMatrixXpr<LMat, T>& A, const T& b)
	{
		return ewise(add_op<T>(), A.derived(), b);
	}

	template<typename T, class RMat>
	LMAT_ENSURE_INLINE
	inline typename binary_fix1_ewise_expr_map<add_op<T>, RMat>::type
	operator + (const T& a, const IMatrixXpr<RMat, T>& B)
	{
		return ewise(add_op<T>(), a, B.derived().derived());
	}

	template<typename T, class LMat, class RMat>
	LMAT_ENSURE_INLINE
	inline LMat& operator += (IDenseMatrix<LMat, T>& A, const IMatrixXpr<RMat, T>& B)
	{
		A.derived() = A.derived() + B.derived();
		return A.derived();
	}

	template<typename T, class LMat>
	LMAT_ENSURE_INLINE
	inline LMat& operator += (IDenseMatrix<LMat, T>& A, const T& b)
	{
		A.derived() = A.derived() + b;
		return A.derived();
	}


	// Subtraction

	template<typename T, class LMat, class RMat>
	LMAT_ENSURE_INLINE
	inline typename binary_ewise_expr_map<sub_op<T>, LMat, RMat>::type
	operator - (const IMatrixXpr<LMat, T>& A, const IMatrixXpr<RMat, T>& B)
	{
		return ewise(sub_op<T>(), A.derived(), B.derived());
	}

	template<typename T, class LMat>
	LMAT_ENSURE_INLINE
	inline typename binary_fix2_ewise_expr_map<sub_op<T>, LMat>::type
	operator - (const IMatrixXpr<LMat, T>& A, const T& b)
	{
		return ewise(sub_op<T>(), A.derived(), b);
	}

	template<typename T, class RMat>
	LMAT_ENSURE_INLINE
	inline typename binary_fix1_ewise_expr_map<sub_op<T>, RMat>::type
	operator - (const T& a, const IMatrixXpr<RMat, T>& B)
	{
		return ewise(sub_op<T>(), a, B.derived());
	}

	template<typename T, class LMat, class RMat>
	LMAT_ENSURE_INLINE
	inline LMat& operator -= (IDenseMatrix<LMat, T>& A, const IMatrixXpr<RMat, T>& B)
	{
		A.derived() = A.derived() - B.derived();
		return A.derived();
	}

	template<typename T, class LMat>
	LMAT_ENSURE_INLINE
	inline LMat& operator -= (IDenseMatrix<LMat, T>& A, const T& b)
	{
		A.derived() = A.derived() - b;
		return A.derived();
	}


	// Multiplication

	template<typename T, class LMat, class RMat>
	LMAT_ENSURE_INLINE
	inline typename binary_ewise_expr_map<mul_op<T>, LMat, RMat>::type
	operator * (const IMatrixXpr<LMat, T>& A, const IMatrixXpr<RMat, T>& B)
	{
		return ewise(mul_op<T>(), A.derived(), B.derived());
	}

	template<typename T, class LMat>
	LMAT_ENSURE_INLINE
	inline typename binary_fix2_ewise_expr_map<mul_op<T>, LMat>::type
	operator * (const IMatrixXpr<LMat, T>& A, const T& b)
	{
		return ewise(mul_op<T>(), A.derived(), b);
	}

	template<typename T, class RMat>
	LMAT_ENSURE_INLINE
	inline typename binary_fix1_ewise_expr_map<mul_op<T>, RMat>::type
	operator * (const T& a, const IMatrixXpr<RMat, T>& B)
	{
		return ewise(mul_op<T>(), a, B.derived());
	}

	template<typename T, class LMat, class RMat>
	LMAT_ENSURE_INLINE
	inline LMat& operator *= (IDenseMatrix<LMat, T>& A, const IMatrixXpr<RMat, T>& B)
	{
		A.derived() = A.derived() * B.derived();
		return A.derived();
	}

	template<typename T, class LMat>
	LMAT_ENSURE_INLINE
	inline LMat& operator *= (IDenseMatrix<LMat, T>& A, const T& b)
	{
		A.derived() = A.derived() * b;
		return A.derived();
	}

	// Division

	template<typename T, class LMat, class RMat>
	LMAT_ENSURE_INLINE
	inline typename binary_ewise_expr_map<div_op<T>, LMat, RMat>::type
	operator / (const IMatrixXpr<LMat, T>& A, const IMatrixXpr<RMat, T>& B)
	{
		return ewise(div_op<T>(), A.derived(), B.derived());
	}

	template<typename T, class LMat>
	LMAT_ENSURE_INLINE
	inline typename binary_fix2_ewise_expr_map<div_op<T>, LMat>::type
	operator / (const IMatrixXpr<LMat, T>& A, const T& b)
	{
		return ewise(div_op<T>(), A.derived(), b);
	}

	template<typename T, class RMat>
	LMAT_ENSURE_INLINE
	inline typename binary_fix1_ewise_expr_map<div_op<T>, RMat>::type
	operator / (const T& a, const IMatrixXpr<RMat, T>& B)
	{
		return ewise(div_op<T>(), a, B.derived());
	}

	template<typename T, class LMat, class RMat>
	LMAT_ENSURE_INLINE
	inline LMat& operator /= (IDenseMatrix<LMat, T>& A, const IMatrixXpr<RMat, T>& B)
	{
		A.derived() = A.derived() / B.derived();
		return A.derived();
	}

	template<typename T, class LMat>
	LMAT_ENSURE_INLINE
	inline LMat& operator /= (IDenseMatrix<LMat, T>& A, const T& b)
	{
		A.derived() = A.derived() / b;
		return A.derived();
	}


	/********************************************
	 *
	 *  Other arithmetic functions
	 *
	 ********************************************/

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline typename unary_ewise_expr_map<neg_op<T>, Mat>::type
	operator - (const IMatrixXpr<Mat, T>& A)
	{
		return ewise(neg_op<T>(), A.derived());
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline typename unary_ewise_expr_map<abs_op<T>, Mat>::type
	abs(const IMatrixXpr<Mat, T>& A)
	{
		return ewise(abs_op<T>(), A.derived());
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline typename unary_ewise_expr_map<sqr_op<T>, Mat>::type
	sqr(const IMatrixXpr<Mat, T>& A)
	{
		return ewise(sqr_op<T>(), A.derived());
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline typename unary_ewise_expr_map<sqrt_op<T>, Mat>::type
	sqrt(const IMatrixXpr<Mat, T>& A)
	{
		return ewise(sqrt_op<T>(), A.derived());
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline typename unary_ewise_expr_map<rcp_op<T>, Mat>::type
	rcp(const IMatrixXpr<Mat, T>& A)
	{
		return ewise(rcp_op<T>(), A.derived());
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline typename unary_ewise_expr_map<rsqrt_op<T>, Mat>::type
	rsqrt(const IMatrixXpr<Mat, T>& A)
	{
		return ewise(rsqrt_op<T>(), A.derived());
	}

}



#endif /* MATRIX_ARITH_H_ */
