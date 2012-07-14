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
	 *  Expression Type mapping
	 *
	 ********************************************/

	LMAT_DECLARE_BINARY_TYPE_MAP_EX( add, add_op )
	LMAT_DECLARE_BINARY_TYPE_MAP_EX( sub, sub_op )
	LMAT_DECLARE_BINARY_TYPE_MAP_EX( mul, mul_op )
	LMAT_DECLARE_BINARY_TYPE_MAP_EX( div, div_op )

	LMAT_DECLARE_UNARY_TYPE_MAP( neg, neg_op )
	LMAT_DECLARE_UNARY_TYPE_MAP( abs, abs_fun )
	LMAT_DECLARE_UNARY_TYPE_MAP( sqr, sqr_fun )
	LMAT_DECLARE_UNARY_TYPE_MAP( sqrt, sqrt_fun )
	LMAT_DECLARE_UNARY_TYPE_MAP( rcp, rcp_fun )
	LMAT_DECLARE_UNARY_TYPE_MAP( rsqrt, rsqrt_fun )

	LMAT_DECLARE_BINARY_TYPE_MAP_EX( max, max_fun )
	LMAT_DECLARE_BINARY_TYPE_MAP_EX( min, min_fun )


	/********************************************
	 *
	 *  Specific Expressions
	 *
	 ********************************************/

	LMAT_DEFINE_BINARY_MATFUNCTION_EX( add, operator +, add_op )
	LMAT_DEFINE_BINARY_MATFUNCTION_EX( sub, operator -, sub_op )
	LMAT_DEFINE_BINARY_MATFUNCTION_EX( mul, operator *, mul_op )
	LMAT_DEFINE_BINARY_MATFUNCTION_EX( div, operator /, div_op )

	LMAT_DEFINE_UNARY_MATFUNCTION( neg, operator -, neg_op )

	LMAT_DEFINE_UNARY_MATFUNCTION( abs, abs, abs_fun )
	LMAT_DEFINE_UNARY_MATFUNCTION( sqr, sqr, sqr_fun )
	LMAT_DEFINE_UNARY_MATFUNCTION( sqrt, sqrt, sqrt_fun )
	LMAT_DEFINE_UNARY_MATFUNCTION( rcp, rcp, rcp_fun )
	LMAT_DEFINE_UNARY_MATFUNCTION( rsqrt, rsqrt, rsqrt_fun )

	LMAT_DEFINE_BINARY_MATFUNCTION_EX( max, max, max_fun )
	LMAT_DEFINE_BINARY_MATFUNCTION_EX( min, min, min_fun )

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
