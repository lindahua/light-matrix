/**
 * @file matrix_veval.h
 *
 * Generic vector-based matrix evaluation
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_VEVAL_H_
#define LIGHTMAT_MATRIX_VEVAL_H_

#include <light_mat/matrix/vec_evaluators.h>

namespace lmat
{

	const int VEC_EVAL_CACHE_COST = 1000;
	const int SHORTVEC_LENGTH_THRESHOLD = 4;
	const int SHORTVEC_PERCOL_COST = 200;

	/********************************************
	 *
	 *  Evaluator dispatch
	 *
	 ********************************************/

	template<class Expr>
	struct linear_eval
	{
		static const bool can_direct = is_dense_mat<Expr>::value && ct_has_continuous_layout<Expr>::value;

		typedef typename matrix_traits<Expr>::value_type T;

		typedef typename
				if_c<can_direct,
					continuous_linear_evaluator<T>,
					cached_linear_evaluator<T> >::type evaluator_type;

		static const int cost = can_direct ? 0 : VEC_EVAL_CACHE_COST;
	};

	template<class Expr>
	struct percol_eval
	{
		static const bool can_direct = is_dense_mat<Expr>::value;
		static const bool has_short_col = ct_rows<Expr>::value < SHORTVEC_LENGTH_THRESHOLD;

		typedef typename matrix_traits<Expr>::value_type T;

		typedef typename
				if_c<can_direct,
					dense_percol_evaluator<T>,
					cached_percol_evaluator<T> >::type evaluator_type;

		static const int normal_cost = can_direct ? 0 : VEC_EVAL_CACHE_COST;
		static const int shortv_cost = SHORTVEC_PERCOL_COST + normal_cost;
	};


	template<class Expr>
	struct linear_eval<embed_mat<Expr> > { };

	template<class Expr>
	struct percol_eval<embed_mat<Expr> > { };


	template<typename T, int CTRows, int CTCols>
	struct linear_eval<const_matrix<T, CTRows, CTCols> >
	{
		typedef const_linear_evaluator<T> evaluator_type;

		static const int cost = 0;
	};

	template<typename T, int CTRows, int CTCols>
	struct percol_eval<const_matrix<T, CTRows, CTCols> >
	{
		typedef const_percol_evaluator<T> evaluator_type;

		static const int normal_cost = 0;
		static const int shortv_cost = 0;
	};


	/********************************************
	 *
	 *  Evaluation functions
	 *
	 ********************************************/

	template<typename T, class Expr, class Dst>
	LMAT_ENSURE_INLINE
	void linear_scalar_evaluate(const IMatrixXpr<Expr, T>& expr, IDenseMatrix<Dst, T>& dst)
	{
		linear_scalar_evalctx<Expr, Dst>::evaluate(expr.derived(), dst.derived());
	}

	template<typename T, class Expr, class Dst>
	LMAT_ENSURE_INLINE
	void percol_scalar_evaluate(const IMatrixXpr<Expr, T>& expr, IDenseMatrix<Dst, T>& dst)
	{
		percol_scalar_evalctx<Expr, Dst>::evaluate(expr.derived(), dst.derived());
	}

}

#endif
