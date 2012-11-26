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

#include <light_mat/matexpr/matrix_ewise_eval.h>
#include <light_mat/math/arith_functors.h>

namespace lmat
{

	/********************************************
	 *
	 *  Expressions
	 *
	 ********************************************/

	LMAT_DEFINE_BINARY_MATFUNCTION( operator +, add_t )
	LMAT_DEFINE_BINARY_MATFUNCTION( operator -, subtract_t )
	LMAT_DEFINE_BINARY_MATFUNCTION( operator *, multiply_t )
	LMAT_DEFINE_BINARY_MATFUNCTION( operator /, divide_t )

	LMAT_DEFINE_UNARY_MATFUNCTION( operator -, negate_t )

	LMAT_DEFINE_UNARY_MATFUNCTION( abs, abs_t )

	LMAT_DEFINE_BINARY_MATFUNCTION( (max), max_t )
	LMAT_DEFINE_BINARY_MATFUNCTION( (min), min_t )

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

