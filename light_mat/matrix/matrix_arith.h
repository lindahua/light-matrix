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
	inline binary_ewise_expr<add_op<T>, LMat, RMat>
	operator + (const IMatrixXpr<LMat, T>& A, const IMatrixXpr<RMat, T>& B)
	{
		return ewise(add_op<T>(), A.derived(), B.derived());
	}

	template<typename T, class LMat>
	LMAT_ENSURE_INLINE
	inline binary_fix2_ewise_expr<add_op<T>, LMat>
	operator + (const IMatrixXpr<LMat, T>& A, const T& b)
	{
		return ewise(add_op<T>(), A.derived(), b);
	}

	template<typename T, class RMat>
	LMAT_ENSURE_INLINE
	inline binary_fix1_ewise_expr<add_op<T>, RMat>
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
	inline binary_ewise_expr<sub_op<T>, LMat, RMat>
	operator - (const IMatrixXpr<LMat, T>& A, const IMatrixXpr<RMat, T>& B)
	{
		return ewise(sub_op<T>(), A.derived(), B.derived());
	}

	template<typename T, class LMat>
	LMAT_ENSURE_INLINE
	inline binary_fix2_ewise_expr<sub_op<T>, LMat>
	operator - (const IMatrixXpr<LMat, T>& A, const T& b)
	{
		return ewise(sub_op<T>(), A.derived(), b);
	}

	template<typename T, class RMat>
	LMAT_ENSURE_INLINE
	inline binary_fix1_ewise_expr<sub_op<T>, RMat>
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
	inline binary_ewise_expr<mul_op<T>, LMat, RMat>
	operator * (const IMatrixXpr<LMat, T>& A, const IMatrixXpr<RMat, T>& B)
	{
		return ewise(mul_op<T>(), A.derived(), B.derived());
	}

	template<typename T, class LMat>
	LMAT_ENSURE_INLINE
	inline binary_fix2_ewise_expr<mul_op<T>, LMat>
	operator * (const IMatrixXpr<LMat, T>& A, const T& b)
	{
		return ewise(mul_op<T>(), A.derived(), b);
	}

	template<typename T, class RMat>
	LMAT_ENSURE_INLINE
	inline binary_fix1_ewise_expr<mul_op<T>, RMat>
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
	inline binary_ewise_expr<div_op<T>, LMat, RMat>
	operator / (const IMatrixXpr<LMat, T>& A, const IMatrixXpr<RMat, T>& B)
	{
		return ewise(div_op<T>(), A.derived(), B.derived());
	}

	template<typename T, class LMat>
	LMAT_ENSURE_INLINE
	inline binary_fix2_ewise_expr<div_op<T>, LMat>
	operator / (const IMatrixXpr<LMat, T>& A, const T& b)
	{
		return ewise(div_op<T>(), A.derived(), b);
	}

	template<typename T, class RMat>
	LMAT_ENSURE_INLINE
	inline binary_fix1_ewise_expr<div_op<T>, RMat>
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
	inline unary_ewise_expr<neg_op<T>, Mat>
	operator - (const IMatrixXpr<Mat, T>& A)
	{
		return ewise(neg_op<T>(), A.derived());
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline unary_ewise_expr<abs_op<T>, Mat>
	abs(const IMatrixXpr<Mat, T>& A)
	{
		return ewise(abs_op<T>(), A.derived());
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline unary_ewise_expr<sqr_op<T>, Mat>
	sqr(const IMatrixXpr<Mat, T>& A)
	{
		return ewise(sqr_op<T>(), A.derived());
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline unary_ewise_expr<sqrt_op<T>, Mat>
	sqrt(const IMatrixXpr<Mat, T>& A)
	{
		return ewise(sqrt_op<T>(), A.derived());
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline unary_ewise_expr<rcp_op<T>, Mat>
	rcp(const IMatrixXpr<Mat, T>& A)
	{
		return ewise(rcp_op<T>(), A.derived());
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline unary_ewise_expr<rsqrt_op<T>, Mat>
	rsqrt(const IMatrixXpr<Mat, T>& A)
	{
		return ewise(rsqrt_op<T>(), A.derived());
	}

}



#endif /* MATRIX_ARITH_H_ */
