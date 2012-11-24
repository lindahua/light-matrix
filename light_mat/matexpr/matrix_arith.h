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

#include <light_mat/matrix/matrix_ewise_expr.h>
#include <light_mat/math/arith_functors.h>

namespace lmat
{

	/********************************************
	 *
	 *  Specific Expressions
	 *
	 ********************************************/

	LMAT_DEFINE_BINARY_MATFUNCTION_EX( operator +, add_op )
	LMAT_DEFINE_BINARY_MATFUNCTION_EX( operator -, sub_op )
	LMAT_DEFINE_BINARY_MATFUNCTION_EX( operator *, mul_op )
	LMAT_DEFINE_BINARY_MATFUNCTION_EX( operator /, div_op )

	LMAT_DEFINE_UNARY_MATFUNCTION( operator -, neg_op )

	LMAT_DEFINE_UNARY_MATFUNCTION( abs, abs_fun )
	LMAT_DEFINE_UNARY_MATFUNCTION( sqr, sqr_fun )
	LMAT_DEFINE_UNARY_MATFUNCTION( sqrt, sqrt_fun )
	LMAT_DEFINE_UNARY_MATFUNCTION( rcp, rcp_fun )
	LMAT_DEFINE_UNARY_MATFUNCTION( rsqrt, rsqrt_fun )

	LMAT_DEFINE_BINARY_MATFUNCTION_EX( max, max_fun )
	LMAT_DEFINE_BINARY_MATFUNCTION_EX( min, min_fun )

	/********************************************
	 *
	 *  In-place Expressions
	 *
	 ********************************************/

	template<typename T, class LArg, class RArg>
	LArg& operator += (IDenseMatrix<LArg, T>& A, const IMatrixXpr<RArg, T>& B)
	{
		A.derived() = A.derived() + B.derived();
		return A.derived();
	}

	template<typename T, class LArg>
	LArg& operator += (IDenseMatrix<LArg, T>& A, const T& b)
	{
		A.derived() = A.derived() + b;
		return A.derived();
	}

	template<typename T, class LArg, class RArg>
	LArg& operator -= (IDenseMatrix<LArg, T>& A, const IMatrixXpr<RArg, T>& B)
	{
		A.derived() = A.derived() - B.derived();
		return A.derived();
	}

	template<typename T, class LArg>
	LArg& operator -= (IDenseMatrix<LArg, T>& A, const T& b)
	{
		A.derived() = A.derived() - b;
		return A.derived();
	}

	template<typename T, class LArg, class RArg>
	LArg& operator *= (IDenseMatrix<LArg, T>& A, const IMatrixXpr<RArg, T>& B)
	{
		A.derived() = A.derived() * B.derived();
		return A.derived();
	}

	template<typename T, class LArg>
	LArg& operator *= (IDenseMatrix<LArg, T>& A, const T& b)
	{
		A.derived() = A.derived() * b;
		return A.derived();
	}

	template<typename T, class LArg, class RArg>
	LArg& operator /= (IDenseMatrix<LArg, T>& A, const IMatrixXpr<RArg, T>& B)
	{
		A.derived() = A.derived() / B.derived();
		return A.derived();
	}

	template<typename T, class LArg>
	LArg& operator /= (IDenseMatrix<LArg, T>& A, const T& b)
	{
		A.derived() = A.derived() / b;
		return A.derived();
	}


}



#endif /* MATRIX_ARITH_H_ */
