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

#ifndef LIGHTMAT_MATRIX_EVAL_H_
#define LIGHTMAT_MATRIX_EVAL_H_

#include <light_mat/matrix/matrix_vec_evaluators.h>
#include <light_mat/matrix/matrix_eval_cost.h>

namespace lmat
{
	namespace detail
	{
		template<typename T, int CTRows>
		struct coleval_internal
		{
			template<class Expr, class Dst>
			LMAT_ENSURE_INLINE
			static void evaluate(const Expr& expr, Dst& dst)
			{
				typename linear_vector_evaluator<Expr>::type evaluator(expr);
				linear_eval_context<T, CTRows> ctx(dst);
				ctx.eval_by_scalars(evaluator);
			}
		};


		template<typename T, int CTSize>
		struct linear_mateval_internal
		{
			template<class Expr, class Dst>
			LMAT_ENSURE_INLINE
			static void evaluate(const Expr& expr, Dst& dst)
			{
				typename linear_vector_evaluator<Expr>::type evaluator(expr);
				linear_eval_context<T, CTSize> ctx(dst);
				ctx.eval_by_scalars(evaluator);
			}
		};

		template<typename T, int CTRows, int CTCols>
		struct percol_mateval_internal
		{
			template<class Expr, class Dst>
			LMAT_ENSURE_INLINE
			static void evaluate(const Expr& expr, Dst& dst)
			{
				typename percol_vector_evaluator<Expr>::type evaluator(expr);
				percol_eval_context<T, CTRows, CTCols> ctx(dst);
				ctx.eval_by_scalars(evaluator);
			}
		};

		template<typename T, int CTRows, int CTCols>
		struct mateval_internal
		{
			template<class Expr, class Dst>
			LMAT_ENSURE_INLINE
			static void evaluate(const Expr& expr, Dst& dst)
			{
				if (!has_continuous_layout<Dst>::value ||
						percol_eval_cost<Expr>::of(expr) < linear_eval_cost<Expr>::of(expr))
				{
					percol_mateval_internal<T, CTRows, CTCols>::evaluate(expr, dst);
				}
				else
				{
					linear_mateval_internal<T, CTRows * CTCols>::evaluate(expr, dst);
				}
			}
		};
	}


	template<typename T, class SExpr, class DMat>
	inline void evaluate_by_scalars(const IMatrixXpr<SExpr, T>& expr, IDenseMatrix<DMat, T>& dst)
	{
		if (is_column(expr))
		{
			detail::coleval_internal<T,
				binary_ct_rows<SExpr, DMat>::value>::evaluate(expr.derived(), dst.derived());
		}
		else
		{
			detail::mateval_internal<T,
				binary_ct_rows<SExpr, DMat>::value,
				binary_ct_cols<SExpr, DMat>::value>::evaluate(expr.derived(), dst.derived());
		}
	}

	template<typename T, class SExpr, class DMat>
	LMAT_ENSURE_INLINE
	inline void evaluate_to(const IMatrixXpr<SExpr, T>& expr, IDenseMatrix<DMat, T>& dst)
	{
		evaluate_by_scalars(expr, dst);
	}

}

#endif /* MATRIX_EVAL_H_ */
