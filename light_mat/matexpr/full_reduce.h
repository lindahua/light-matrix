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

	template<class Fun, class Expr, typename KerCate>
	LMAT_ENSURE_INLINE
	inline static typename Fun::result_type reduce(const Fun& fun,
			const IMatrixXpr<Expr, typename Fun::arg_type>& expr,
			matrix_visit_policy<linear_vis, KerCate>)
	{
		typedef matrix_visit_setting<linear_vis, KerCate,
				ct_rows<Expr>::value,
				ct_cols<Expr>::value> setting_t;
		typename matrix_vismap<Expr, setting_t>::type visitor(expr.derived());

		return detail::single_vec_reduce<ct_size<Expr>::value, KerCate>::eval(
				fun, expr.nelems(), visitor);
	}


	template<class Fun, class Expr, typename KerCate>
	inline static typename Fun::result_type reduce(const Fun& fun,
			const IMatrixXpr<Expr, typename Fun::arg_type>& expr,
			matrix_visit_policy<percol_vis, KerCate>)
	{
		typedef matrix_visit_setting<percol_vis, KerCate,
				ct_rows<Expr>::value,
				ct_cols<Expr>::value> setting_t;
		typedef typename matrix_vismap<Expr, setting_t>::type visitor_t;
		visitor_t visitor(expr.derived());

		const index_t m = expr.nrows();
		const index_t n = expr.ncolumns();

		typedef typename Fun::result_type RT;

		RT r = detail::single_vec_reduce<ct_rows<Expr>::value, KerCate>::eval(
				fun, m, visitor, visitor.col_state(0));

		for (index_t j = 1; j < n; ++j)
		{
			typename matrix_visitor_state<visitor_t>::type s = visitor.col_state(j);

			RT rj = detail::single_vec_reduce<ct_rows<Expr>::value, KerCate>::eval(
					fun, m, visitor, s);

			r = fun(r, rj);
		}

		return r;
	}

	template<class Fun, typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline typename Fun::result_type
	linear_by_scalars_reduce(const Fun& fun, const IMatrixXpr<Mat, T>& X)
	{
		return reduce(fun, X.derived(), matrix_visit_policy<linear_vis, scalar_kernel_t>());
	}

	template<class Fun, typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline typename Fun::result_type
	percol_by_scalars_reduce(const Fun& fun, const IMatrixXpr<Mat, T>& X)
	{
		return reduce(fun, X.derived(), matrix_visit_policy<percol_vis, scalar_kernel_t>());
	}


	template<class Mat>
	struct default_matrix_reduce_policy
	{
		typedef dense_matrix<
				typename matrix_traits<Mat>::value_type,
				ct_rows<Mat>::value,
				ct_cols<Mat>::value> Dst;

		typedef typename default_matrix_visit_policy<Mat, Dst>::type type;
	};


	template<class Fun, typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline typename Fun::result_type
	reduce(const Fun& fun, const IMatrixXpr<Mat, T>& X)
	{
		return reduce(fun, X.derived(),
				typename default_matrix_reduce_policy<Mat>::type());
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


