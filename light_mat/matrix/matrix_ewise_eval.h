/**
 * @file matrix_ewise_eval.h
 *
 * Evaluation of element-wise matrix expressions
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_EWISE_EVAL_H_
#define LIGHTMAT_MATRIX_EWISE_EVAL_H_

#include <light_mat/matrix/matrix_ewise_expr.h>
#include <light_mat/matrix/matrix_eval.h>

namespace lmat
{
	template<class Fun, class Arg>
	class simple_unary_ewise_linear_evaluator
	: public ILinearVectorEvaluator<
	  	  simple_unary_ewise_linear_evaluator<Fun, Arg>,
	  	  typename Fun::result_type>
	{
		simple_unary_ewise_linear_evaluator(const simple_unary_ewise_expr<Fun, Arg>& expr)
		: m_fun(), m_arg_eval(expr)
		{
		}

		LMAT_ENSURE_INLINE T get_value(const index_t i) const
		{
			return m_fun(m_arg_eval.get_value(i));
		}

	private:
		Fun m_fun;
		typedef typename linear_vector_evaluator<Arg>::type arg_eval_t;
		arg_eval_t m_arg_eval;
	};

	template<class Fun, class Arg>
	class simple_unary_ewise_percol_evaluator
	: public ILinearVectorEvaluator<
	  	  simple_unary_ewise_percol_evaluator<Fun, Arg>,
	  	  typename Fun::result_type>
	{
		simple_unary_ewise_percol_evaluator(const simple_unary_ewise_expr<Fun, Arg>& expr)
		: m_fun(), m_arg_eval(expr)
		{
		}

		LMAT_ENSURE_INLINE T get_value(const index_t i) const
		{
			return m_fun(m_arg_eval.get_value(i));
		}

	private:
		Fun m_fun;
		typedef typename percol_vector_evaluator<Arg>::type arg_eval_t;
		arg_eval_t m_arg_eval;
	};


}

#endif /* MATRIX_EWISE_EVAL_H_ */
