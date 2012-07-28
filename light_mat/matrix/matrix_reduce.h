/**
 * @file matrix_reduce.h
 *
 * Matrix reduction expression
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_REDUCE_H_
#define LIGHTMAT_MATRIX_REDUCE_H_

#include <light_mat/math/reduction_functors.h>
#include <light_mat/matrix/matrix_arith.h>
#include <light_mat/matrix/matrix_ewise_eval.h>

#include "bits/matrix_reduce_internal.h"

namespace lmat
{

	/********************************************
	 *
	 *  generic reduction functions
	 *
	 ********************************************/

	template<class Fun, typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline typename Fun::result_type
	linear_scalar_reduce(const Fun& fun, const IMatrixXpr<Mat, T>& X)
	{
		return detail::full_reduce_linear_internal::evaluate(fun, X.derived());
	}

	template<class Fun, typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline typename Fun::result_type
	percol_scalar_reduce(const Fun& fun, const IMatrixXpr<Mat, T>& X)
	{
		return detail::full_reduce_percol_internal::evaluate(fun, X.derived());
	}

	template<class Fun, typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline typename Fun::result_type
	reduce(const Fun& fun, const IMatrixXpr<Mat, T>& X)
	{
		typedef typename detail::reduce_internal_map<Mat>::type internal_t;
		return internal_t::evaluate(fun, X.derived());
	}


	/********************************************
	 *
	 *  specific reduction functions
	 *
	 ********************************************/

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	T sum(const IMatrixXpr<Mat, T>& X)
	{
		return reduce(sum_fun<T>(), X);
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	T mean(const IMatrixXpr<Mat, T>& X)
	{
		return sum(X.derived()) / T(X.nelems());
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	T prod(const IMatrixXpr<Mat, T>& X)
	{
		return reduce(prod_fun<T>(), X);
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	T maximum(const IMatrixXpr<Mat, T>& X)
	{
		return reduce(maximum_fun<T>(), X);
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	T minimum(const IMatrixXpr<Mat, T>& X)
	{
		return reduce(minimum_fun<T>(), X);
	}

	template<typename T, class LMat, class RMat>
	LMAT_ENSURE_INLINE
	T dot(const IMatrixXpr<LMat, T>& X, const IMatrixXpr<RMat, T>& Y)
	{
		return sum(X.derived() * Y.derived());
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	T L1norm(const IMatrixXpr<Mat, T>& X)
	{
		return sum(abs(X));
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	T sqL2norm(const IMatrixXpr<Mat, T>& X)
	{
		return sum(sqr(X));
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	T L2norm(const IMatrixXpr<Mat, T>& X)
	{
		return math::sqrt(sqL2norm(X));
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	T Linfnorm(const IMatrixXpr<Mat, T>& X)
	{
		return X.nelems() > 0 ? maximum(abs(X)) : T(0);
	}

	template<typename T, class LMat, class RMat>
	LMAT_ENSURE_INLINE
	T nrmdot(const IMatrixXpr<LMat, T>& X, const IMatrixXpr<RMat, T>& Y)
	{
		return dot(X, Y) / (L2norm(X) * L2norm(Y));
	}
}

#endif /* MATRIX_REDUC_EXPR_H_ */


