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
#include <light_mat/matrix/matrix_vector_eval.h>

namespace lmat
{

	/********************************************
	 *
	 *  Vector evaluator classes
	 *
	 ********************************************/

	template<class Fun, typename Arg_HP, class Arg>
	class unary_ewise_linear_evaluator
	: public ILinearVectorEvaluator<
	  	  unary_ewise_linear_evaluator<Fun, Arg_HP, Arg>,
	  	  typename Fun::result_type>
	{
	public:
		typedef unary_ewise_expr<Fun, Arg_HP, Arg> expr_t;
		typedef typename expr_t::arg_type arg_type;

		LMAT_ENSURE_INLINE
		unary_ewise_linear_evaluator(const expr_t& expr)
		: m_fun(expr.fun()), m_arg_eval(expr.arg())
		{
		}

		LMAT_ENSURE_INLINE typename Fun::result_type get_value(const index_t i) const
		{
			return m_fun(m_arg_eval.get_value(i));
		}

	private:
		Fun m_fun;
		typedef typename vector_eval<arg_type, as_linear_vec, by_scalars>::evaluator_type arg_eval_t;
		arg_eval_t m_arg_eval;
	};

	template<class Fun, typename Arg_HP, class Arg>
	class unary_ewise_percol_evaluator
	: public IPerColVectorEvaluator<
	  	  unary_ewise_percol_evaluator<Fun, Arg_HP, Arg>,
	  	  typename Fun::result_type>
	{
	public:
		typedef unary_ewise_expr<Fun, Arg_HP, Arg> expr_t;
		typedef typename expr_t::arg_type arg_type;

		LMAT_ENSURE_INLINE
		unary_ewise_percol_evaluator(const expr_t& expr)
		: m_fun(expr.fun()), m_arg_eval(expr.arg())
		{
		}

		LMAT_ENSURE_INLINE typename Fun::result_type get_value(const index_t i) const
		{
			return m_fun(m_arg_eval.get_value(i));
		}

		LMAT_ENSURE_INLINE void next_column()
		{
			m_arg_eval.next_column();
		}

	private:
		Fun m_fun;
		typedef typename vector_eval<arg_type, per_column, by_scalars>::evaluator_type arg_eval_t;
		arg_eval_t m_arg_eval;
	};


	template<class Fun, typename Arg1_HP, class Arg1, typename Arg2_HP, class Arg2>
	class binary_ewise_linear_evaluator
	: public ILinearVectorEvaluator<
	  	  binary_ewise_linear_evaluator<Fun, Arg1_HP, Arg1, Arg2_HP, Arg2>,
	  	  typename Fun::result_type>
	{
	public:
		typedef binary_ewise_expr<Fun, Arg1_HP, Arg1, Arg2_HP, Arg2> expr_t;
		typedef typename expr_t::arg1_type arg1_type;
		typedef typename expr_t::arg2_type arg2_type;

		LMAT_ENSURE_INLINE
		binary_ewise_linear_evaluator(const expr_t& expr)
		: m_fun(expr.fun()), m_arg1_eval(expr.first_arg()), m_arg2_eval(expr.second_arg())
		{
		}

		LMAT_ENSURE_INLINE typename Fun::result_type get_value(const index_t i) const
		{
			return m_fun(m_arg1_eval.get_value(i), m_arg2_eval.get_value(i));
		}

	private:
		Fun m_fun;
		typedef typename vector_eval<arg1_type, as_linear_vec, by_scalars>::evaluator_type arg1_eval_t;
		typedef typename vector_eval<arg2_type, as_linear_vec, by_scalars>::evaluator_type arg2_eval_t;

		arg1_eval_t m_arg1_eval;
		arg2_eval_t m_arg2_eval;
	};


	template<class Fun, typename Arg1_HP, class Arg1, typename Arg2_HP, class Arg2>
	class binary_ewise_percol_evaluator
	: public IPerColVectorEvaluator<
	  	  binary_ewise_percol_evaluator<Fun, Arg1_HP, Arg1, Arg2_HP, Arg2>,
	  	  typename Fun::result_type>
	{
	public:
		typedef binary_ewise_expr<Fun, Arg1_HP, Arg1, Arg2_HP, Arg2> expr_t;
		typedef typename expr_t::arg1_type arg1_type;
		typedef typename expr_t::arg2_type arg2_type;

		LMAT_ENSURE_INLINE
		binary_ewise_percol_evaluator(const expr_t& expr)
		: m_fun(expr.fun()), m_arg1_eval(expr.first_arg()), m_arg2_eval(expr.second_arg())
		{
		}

		LMAT_ENSURE_INLINE typename Fun::result_type get_value(const index_t i) const
		{
			return m_fun(m_arg1_eval.get_value(i), m_arg2_eval.get_value(i));
		}

		LMAT_ENSURE_INLINE void next_column()
		{
			m_arg1_eval.next_column();
			m_arg2_eval.next_column();
		}

	private:
		Fun m_fun;
		typedef typename vector_eval<arg1_type, per_column, by_scalars>::evaluator_type arg1_eval_t;
		typedef typename vector_eval<arg2_type, per_column, by_scalars>::evaluator_type arg2_eval_t;

		arg1_eval_t m_arg1_eval;
		arg2_eval_t m_arg2_eval;
	};



	/********************************************
	 *
	 *  Evaluator dispatch
	 *
	 ********************************************/

	template<typename Fun, typename Arg_HP, class Arg, typename Means>
	struct vector_eval<unary_ewise_expr<Fun, Arg_HP, Arg>, as_linear_vec, Means>
	{
		typedef unary_ewise_linear_evaluator<Fun, Arg_HP, Arg> evaluator_type;

		typedef unary_ewise_expr<Fun, Arg_HP, Arg> expr_t;
		typedef typename expr_t::arg_type arg_type;
		static const int cost = vector_eval<arg_type, as_linear_vec, Means>::cost;
	};

	template<typename Fun, typename Arg_HP, class Arg, typename Means>
	struct vector_eval<unary_ewise_expr<Fun, Arg_HP, Arg>, per_column, Means>
	{
		typedef unary_ewise_percol_evaluator<Fun, Arg_HP, Arg> evaluator_type;

		typedef unary_ewise_expr<Fun, Arg_HP, Arg> expr_t;
		typedef typename expr_t::arg_type arg_type;
		static const int normal_cost = vector_eval<arg_type, per_column, Means>::normal_cost;
		static const int shortv_cost = vector_eval<arg_type, per_column, Means>::shortv_cost;

		static const bool has_short_col = ct_rows<arg_type>::value < SHORTVEC_LENGTH_THRESHOLD;
		static const int cost = has_short_col ? shortv_cost : normal_cost;
	};

	template<typename Fun, typename Arg1_HP, class Arg1, typename Arg2_HP, class Arg2, typename Means>
	struct vector_eval<binary_ewise_expr<Fun, Arg1_HP, Arg1, Arg2_HP, Arg2>, as_linear_vec, Means>
	{
		typedef binary_ewise_linear_evaluator<Fun, Arg1_HP, Arg1, Arg2_HP, Arg2> evaluator_type;

		typedef binary_ewise_expr<Fun, Arg1_HP, Arg1, Arg2_HP, Arg2> expr_t;
		typedef typename expr_t::arg1_type arg1_type;
		typedef typename expr_t::arg2_type arg2_type;

		static const int cost =
				vector_eval<arg1_type, as_linear_vec, Means>::cost +
				vector_eval<arg2_type, as_linear_vec, Means>::cost;
	};

	template<typename Fun, typename Arg1_HP, class Arg1, typename Arg2_HP, class Arg2, typename Means>
	struct vector_eval<binary_ewise_expr<Fun, Arg1_HP, Arg1, Arg2_HP, Arg2>, per_column, Means>
	{
		typedef binary_ewise_percol_evaluator<Fun, Arg1_HP, Arg1, Arg2_HP, Arg2> evaluator_type;

		typedef binary_ewise_expr<Fun, Arg1_HP, Arg1, Arg2_HP, Arg2> expr_t;
		typedef typename expr_t::arg1_type arg1_type;
		typedef typename expr_t::arg2_type arg2_type;

		static const int normal_cost =
				vector_eval<arg1_type, per_column, Means>::normal_cost +
				vector_eval<arg2_type, per_column, Means>::normal_cost;

		static const int shortv_cost =
				vector_eval<arg1_type, per_column, Means>::shortv_cost +
				vector_eval<arg2_type, per_column, Means>::shortv_cost;

		static const bool has_short_col = binary_ct_rows<arg1_type, arg2_type>::value < SHORTVEC_LENGTH_THRESHOLD;
		static const int cost = has_short_col ? shortv_cost : normal_cost;
	};


	template<typename Fun, typename Arg_HP, class Arg, class Dst>
	struct default_matrix_eval_policy<unary_ewise_expr<Fun, Arg_HP, Arg>, Dst>
	{
		typedef typename vector_eval_default_policy<
				unary_ewise_expr<Fun, Arg_HP, Arg> >::type type;
	};

	template<typename Fun, typename Arg1_HP, class Arg1, typename Arg2_HP, class Arg2, class Dst>
	struct default_matrix_eval_policy<binary_ewise_expr<Fun, Arg1_HP, Arg1, Arg2_HP, Arg2>, Dst>
	{
		typedef typename vector_eval_default_policy<
				binary_ewise_expr<Fun, Arg1_HP, Arg1, Arg2_HP, Arg2> >::type type;
	};

}

#endif /* MATRIX_EWISE_EVAL_H_ */
