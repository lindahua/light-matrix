/**
 * @file full_reduce.h
 *
 * Full matrix reduction expression
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_REDUCE_H_
#define LIGHTMAT_MATRIX_REDUCE_H_

#include "bits/full_reduce_internal.h"

#define LMAT_DEFINE_UNARY_FULL_REDUCE_FUNCTION( Tag, FunName ) \
	template<typename T, class Xpr> \
	LMAT_ENSURE_INLINE \
	inline typename reduction_result<Tag, T>::type \
	FunName(const IMatrixXpr<Xpr, T>& a) { return reduce(Tag(), a); }

#define LMAT_DEFINE_BINARY_FULL_REDUCE_FUNCTION( Tag, FunName ) \
	template<typename T, class Xpr1, class Xpr2> \
	LMAT_ENSURE_INLINE \
	inline typename reduction_result<Tag, T, T>::type \
	FunName(const IMatrixXpr<Xpr1, T>& a, const IMatrixXpr<Xpr2, T>& b) { return reduce(Tag(), a, b); }

#define LMAT_DEFINE_UNARY_FULL_REDUCE_FUNCTION_( Name ) \
	LMAT_DEFINE_UNARY_FULL_REDUCE_FUNCTION( Name##_t, Name )

#define LMAT_DEFINE_BINARY_FULL_REDUCE_FUNCTION_( Name ) \
	LMAT_DEFINE_BINARY_FULL_REDUCE_FUNCTION( Name##_t, Name )


namespace lmat
{
	/********************************************
	 *
	 *  generic reduction functions
	 *
	 ********************************************/

	template<typename Tag, typename T, class Xpr, typename Acc, typename Ker>
	inline typename reduction_result<Tag, T>::type
	_reduce(Tag tag, const IMatrixXpr<Xpr, T>& a, Acc, Ker)
	{
		return internal::_reduce(tag, a, Acc(), Ker());
	}

	template<typename Tag, typename T1, class Xpr1, typename T2, class Xpr2, typename Acc, typename Ker>
	inline typename reduction_result<Tag, T1, T2>::type
	_reduce(Tag tag, const IMatrixXpr<Xpr1, T1>& a, const IMatrixXpr<Xpr2, T2>& b, Acc, Ker)
	{
		return internal::_reduce(tag, a, b, Acc(), Ker());
	}


	template<typename Tag, typename T, class Xpr>
	LMAT_ENSURE_INLINE
	inline typename reduction_result<Tag, T>::type
	reduce(Tag tag, const IMatrixXpr<Xpr, T>& a)
	{
		const index_t m = a.nrows();

		int macc_linear_cost = macc_cost<Xpr, linear_macc, scalar_kernel_t>::value;
		int macc_percol_cost = macc_cost<Xpr, percol_macc, scalar_kernel_t>::value;

		if (m <= MACC_SHORT_PERCOL_COST)
			macc_percol_cost += MACC_SHORT_PERCOL_COST;

		if (macc_linear_cost <= macc_percol_cost)
		{
			return _reduce(tag, a, linear_macc(), scalar_kernel_t());
		}
		else
		{
			return _reduce(tag, a, percol_macc(), scalar_kernel_t());
		}
	}


	template<typename Tag, typename T, class Xpr1, class Xpr2>
	LMAT_ENSURE_INLINE
	inline typename reduction_result<Tag, T, T>::type
	reduce(Tag tag, const IMatrixXpr<Xpr1, T>& a, const IMatrixXpr<Xpr2, T>& b)
	{
		const index_t m = a.nrows();

		int macc_linear_cost =
				macc_cost<Xpr1, linear_macc, scalar_kernel_t>::value +
				macc_cost<Xpr2, linear_macc, scalar_kernel_t>::value;

		int macc_percol_cost =
				macc_cost<Xpr1, percol_macc, scalar_kernel_t>::value +
				macc_cost<Xpr2, percol_macc, scalar_kernel_t>::value;

		if (m <= MACC_SHORT_PERCOL_COST)
			macc_percol_cost += MACC_SHORT_PERCOL_COST;

		if (macc_linear_cost <= macc_percol_cost)
		{
			return _reduce(tag, a, b, linear_macc(), scalar_kernel_t());
		}
		else
		{
			return _reduce(tag, a, b, percol_macc(), scalar_kernel_t());
		}
	}


	/********************************************
	 *
	 *  specific reduction
	 *
	 ********************************************/

	LMAT_DEFINE_UNARY_FULL_REDUCE_FUNCTION_( sum )
	LMAT_DEFINE_UNARY_FULL_REDUCE_FUNCTION_( maximum )
	LMAT_DEFINE_UNARY_FULL_REDUCE_FUNCTION_( minimum )
	LMAT_DEFINE_UNARY_FULL_REDUCE_FUNCTION_( mean )

	LMAT_DEFINE_UNARY_FULL_REDUCE_FUNCTION_( L1norm )
	LMAT_DEFINE_UNARY_FULL_REDUCE_FUNCTION_( sqL2norm )
	LMAT_DEFINE_UNARY_FULL_REDUCE_FUNCTION_( L2norm )
	LMAT_DEFINE_UNARY_FULL_REDUCE_FUNCTION_( Linfnorm )

	LMAT_DEFINE_UNARY_FULL_REDUCE_FUNCTION_( logsum )
	LMAT_DEFINE_UNARY_FULL_REDUCE_FUNCTION_( entropy )

	LMAT_DEFINE_BINARY_FULL_REDUCE_FUNCTION_( dot )
	LMAT_DEFINE_BINARY_FULL_REDUCE_FUNCTION_( nrmdot )
}

#endif /* MATRIX_REDUC_EXPR_H_ */


