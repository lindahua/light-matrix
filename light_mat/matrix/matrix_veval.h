/**
 * @file matrix_eval.h
 *
 * Generic matrix evaluation
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_GENERIC_MATRIX_EVAL_H_
#define LIGHTMAT_GENERIC_MATRIX_EVAL_H_

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

		LMAT_ENSURE_INLINE
		static int cost_of(const Expr& )
		{
			return cost;
		}
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

		static const int cost = can_direct ?
				(has_short_col ? SHORTVEC_PERCOL_COST : 0) :
				(has_short_col ? (SHORTVEC_PERCOL_COST + VEC_EVAL_CACHE_COST) : VEC_EVAL_CACHE_COST);


		LMAT_ENSURE_INLINE
		static int cost_of(const Expr& )
		{
			return cost;
		}
	};


	template<class Expr>
	struct linear_eval<embed_mat<Expr> > { };

	template<class Expr>
	struct percol_eval<embed_mat<Expr> > { };


	template<typename T, int CTRows, int CTCols>
	struct linear_eval<const_matrix<T, CTRows, CTCols> >
	{
		typedef const_linear_evaluator<T> evaluator_type;

		LMAT_ENSURE_INLINE
		static int cost_of(const const_matrix<T, CTRows, CTCols>& )
		{
			return 0;
		}
	};

	template<typename T, int CTRows, int CTCols>
	struct percol_eval<const_matrix<T, CTRows, CTCols> >
	{
		typedef const_percol_evaluator<T> evaluator_type;

		LMAT_ENSURE_INLINE
		static int cost_of(const const_matrix<T, CTRows, CTCols>& )
		{
			return 0;
		}
	};


	/********************************************
	 *
	 *  Generic evaluation
	 *
	 ********************************************/

	namespace detail
	{

		template<typename T, int CTSize>
		struct ewise_linear_eval_internal
		{
			template<class Expr, class Dst>
			LMAT_ENSURE_INLINE
			static void evaluate(const Expr& expr, Dst& dst)
			{
				typename linear_eval<Expr>::evaluator_type evaluator(expr);
				linear_eval_context<T, CTSize> ctx(dst);
				ctx.eval_by_scalars(evaluator);
			}
		};

		template<typename T, int CTRows, int CTCols>
		struct ewise_percol_eval_internal
		{
			template<class Expr, class Dst>
			LMAT_ENSURE_INLINE
			static void evaluate(const Expr& expr, Dst& dst)
			{
				typename percol_eval<Expr>::evaluator_type evaluator(expr);
				percol_eval_context<T, CTRows, CTCols> ctx(dst);
				ctx.eval_by_scalars(evaluator);
			}
		};

		template<typename T, int CTRows, int CTCols>
		struct ewise_eval_by_scalars_internal
		{
			template<class Expr, class Dst>
			inline static void evaluate(const Expr& expr, Dst& dst)
			{
				if (!ct_has_continuous_layout<Dst>::value ||
						percol_eval<Expr>::cost_of(expr) < linear_eval<Expr>::cost_of(expr))
				{
					ewise_percol_eval_internal<T, CTRows, CTCols>::evaluate(expr, dst);
				}
				else
				{
					ewise_linear_eval_internal<T, CTRows * CTCols>::evaluate(expr, dst);
				}
			}
		};
	}


	template<typename S, typename T, class SExpr, class DMat>
	LMAT_ENSURE_INLINE
	inline void evaluate_by_scalars(const IMatrixXpr<SExpr, S>& expr, IDenseMatrix<DMat, T>& dst)
	{
		detail::ewise_eval_by_scalars_internal<T,
			binary_ct_rows<SExpr, DMat>::value,
			binary_ct_cols<SExpr, DMat>::value>::evaluate(expr.derived(), dst.derived());
	}

	template<typename T, class SExpr, class DMat>
	LMAT_ENSURE_INLINE
	inline void evaluate_to(const IMatrixXpr<SExpr, T>& expr, IDenseMatrix<DMat, T>& dst)
	{
		evaluate_by_scalars(expr, dst);
	}



	template<typename T, class Expr, class DMat>
	LMAT_ENSURE_INLINE
	inline void evaluate_to(const embed_mat<Expr>& expr, IDenseMatrix<DMat, T>& dst);
	// ensure that evaluate_to is never invoked on embed_mat<Expr>

}

#endif
