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

	template<class Fun, class Arg_Holder>
	class unary_ewise_linear_evaluator
	: public ILinearVectorEvaluator<
	  	  unary_ewise_linear_evaluator<Fun, Arg_Holder>,
	  	  typename Fun::result_type>
	{
	public:
		typedef typename Arg_Holder::arg_type arg_type;

		LMAT_ENSURE_INLINE
		unary_ewise_linear_evaluator(const unary_ewise_expr<Fun, Arg_Holder>& expr)
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

	template<class Fun, class Arg_Holder>
	class unary_ewise_percol_evaluator
	: public IPerColVectorEvaluator<
	  	  unary_ewise_percol_evaluator<Fun, Arg_Holder>,
	  	  typename Fun::result_type>
	{
	public:
		typedef typename Arg_Holder::arg_type arg_type;

		LMAT_ENSURE_INLINE
		unary_ewise_percol_evaluator(const unary_ewise_expr<Fun, Arg_Holder>& expr)
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


	template<class Fun, class Arg1_Holder, class Arg2_Holder>
	class binary_ewise_linear_evaluator
	: public ILinearVectorEvaluator<
	  	  binary_ewise_linear_evaluator<Fun, Arg1_Holder, Arg2_Holder>,
	  	  typename Fun::result_type>
	{
	public:
		typedef typename Arg1_Holder::arg_type arg1_type;
		typedef typename Arg2_Holder::arg_type arg2_type;

		LMAT_ENSURE_INLINE
		binary_ewise_linear_evaluator(const binary_ewise_expr<Fun, Arg1_Holder, Arg2_Holder>& expr)
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


	template<class Fun, class Arg1_Holder, class Arg2_Holder>
	class binary_ewise_percol_evaluator
	: public IPerColVectorEvaluator<
	  	  binary_ewise_percol_evaluator<Fun, Arg1_Holder, Arg2_Holder>,
	  	  typename Fun::result_type>
	{
	public:
		typedef typename Arg1_Holder::arg_type arg1_type;
		typedef typename Arg2_Holder::arg_type arg2_type;

		LMAT_ENSURE_INLINE
		binary_ewise_percol_evaluator(const binary_ewise_expr<Fun, Arg1_Holder, Arg2_Holder>& expr)
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

	template<typename Fun, class Arg_Holder, typename Means>
	struct vector_eval<unary_ewise_expr<Fun, Arg_Holder>, as_linear_vec, Means>
	{
		typedef unary_ewise_linear_evaluator<Fun, Arg_Holder> evaluator_type;

		typedef typename Arg_Holder::arg_type arg_type;
		static const int cost = vector_eval<arg_type, as_linear_vec, Means>::cost;
	};

	template<typename Fun, class Arg_Holder, typename Means>
	struct vector_eval<unary_ewise_expr<Fun, Arg_Holder>, per_column, Means>
	{
		typedef unary_ewise_percol_evaluator<Fun, Arg_Holder> evaluator_type;

		typedef typename Arg_Holder::arg_type arg_type;
		static const int normal_cost = vector_eval<arg_type, per_column, Means>::normal_cost;
		static const int shortv_cost = vector_eval<arg_type, per_column, Means>::shortv_cost;

		static const bool has_short_col = ct_rows<arg_type>::value < SHORTVEC_LENGTH_THRESHOLD;
		static const int cost = has_short_col ? shortv_cost : normal_cost;
	};

	template<typename Fun, class Arg1_Holder, class Arg2_Holder, typename Means>
	struct vector_eval<binary_ewise_expr<Fun, Arg1_Holder, Arg2_Holder>, as_linear_vec, Means>
	{
		typedef binary_ewise_linear_evaluator<Fun, Arg1_Holder, Arg2_Holder> evaluator_type;

		typedef typename Arg1_Holder::arg_type arg1_type;
		typedef typename Arg2_Holder::arg_type arg2_type;

		static const int cost =
				vector_eval<arg1_type, as_linear_vec, Means>::cost +
				vector_eval<arg2_type, as_linear_vec, Means>::cost;
	};

	template<typename Fun, class Arg1_Holder, class Arg2_Holder, typename Means>
	struct vector_eval<binary_ewise_expr<Fun, Arg1_Holder, Arg2_Holder>, per_column, Means>
	{
		typedef binary_ewise_percol_evaluator<Fun, Arg1_Holder, Arg2_Holder> evaluator_type;

		typedef typename Arg1_Holder::arg_type arg1_type;
		typedef typename Arg2_Holder::arg_type arg2_type;

		static const int normal_cost =
				vector_eval<arg1_type, per_column, Means>::normal_cost +
				vector_eval<arg2_type, per_column, Means>::normal_cost;

		static const int shortv_cost =
				vector_eval<arg1_type, per_column, Means>::shortv_cost +
				vector_eval<arg2_type, per_column, Means>::shortv_cost;

		static const bool has_short_col = binary_ct_rows<arg1_type, arg2_type>::value < SHORTVEC_LENGTH_THRESHOLD;
		static const int cost = has_short_col ? shortv_cost : normal_cost;
	};


	template<typename Fun, class Arg_Holder, class Dst>
	struct default_matrix_eval_policy<unary_ewise_expr<Fun, Arg_Holder>, Dst>
	{
		typedef typename vector_eval_default_policy<
				unary_ewise_expr<Fun, Arg_Holder> >::type type;
	};

	template<typename Fun, class Arg1_Holder, class Arg2_Holder, class Dst>
	struct default_matrix_eval_policy<binary_ewise_expr<Fun, Arg1_Holder, Arg2_Holder>, Dst>
	{
		typedef typename vector_eval_default_policy<
				binary_ewise_expr<Fun, Arg1_Holder, Arg2_Holder> >::type type;
	};

}

#endif /* MATRIX_EWISE_EVAL_H_ */
