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

	template<class Fun, class Expr, typename Means>
	LMAT_ENSURE_INLINE
	inline static typename Fun::result_type reduce(const Fun& fun,
			const IMatrixXpr<Expr, typename Fun::arg_type>& expr,
			vector_eval_policy<as_linear_vec, Means>)
	{
		typename vector_eval<Expr, as_linear_vec, Means>::evaluator_type evaluator(expr.derived());
		return detail::single_vec_reduce<ct_size<Expr>::value, Means>::eval(
				fun, expr.nelems(), evaluator);
	}


	template<class Fun, class Expr, typename Means>
	inline static typename Fun::result_type reduce(const Fun& fun,
			const IMatrixXpr<Expr, typename Fun::arg_type>& expr,
			vector_eval_policy<per_column, Means>)
	{
		typedef typename vector_eval<Expr, per_column, Means>::evaluator_type eval_t;
		eval_t evaluator(expr.derived());
		const index_t m = expr.nrows();
		const index_t n = expr.ncolumns();

		typedef typename Fun::result_type RT;

		RT r = detail::single_vec_reduce<ct_rows<Expr>::value, Means>::eval(
				fun, m, evaluator, evaluator.col_state(0));

		for (index_t j = 1; j < n; ++j)
		{
			typename percol_eval_state<eval_t>::type s = evaluator.col_state(j);

			RT rj = detail::single_vec_reduce<ct_rows<Expr>::value, Means>::eval(
					fun, m, evaluator, s);

			r = fun(r, rj);
		}

		return r;
	}

	template<class Fun, typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline typename Fun::result_type
	linear_by_scalars_reduce(const Fun& fun, const IMatrixXpr<Mat, T>& X)
	{
		return reduce(fun, X.derived(), vector_eval_policy<as_linear_vec, by_scalars>());
	}

	template<class Fun, typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline typename Fun::result_type
	percol_by_scalars_reduce(const Fun& fun, const IMatrixXpr<Mat, T>& X)
	{
		return reduce(fun, X.derived(), vector_eval_policy<per_column, by_scalars>());
	}


	template<class Fun, typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline typename Fun::result_type
	reduce(const Fun& fun, const IMatrixXpr<Mat, T>& X)
	{
		return reduce(fun, X.derived(), typename vector_eval_default_policy<Mat>::type());
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


