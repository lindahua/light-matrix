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
#include <light_mat/matrix/matrix_visit_eval.h>
#include <utility>

namespace lmat
{

	// forward

	template<typename VisSetting, class Fun, class Arg>
	class unary_ewise_mvisitor;

	template<typename VisSetting, class Fun, class Arg1, class Arg2>
	class binary_ewise_mvisitor;


	/********************************************
	 *
	 *  State maps
	 *
	 ********************************************/

	template<typename KerCate, int ME, int NE, class Arg>
	struct unary_ewise_percol_vstate
	{
		typedef matrix_visit_setting<percol_vis, KerCate, ME, NE> setting_t;

		typedef typename matrix_visitor_state<
					typename matrix_vismap<Arg, setting_t>::type
				>::type type;
	};

	template<typename KerCate, int ME, int NE, class Arg1, class Arg2>
	struct binary_ewise_percol_vstate
	{
		typedef matrix_visit_setting<percol_vis, KerCate, ME, NE> setting_t;

		typedef typename matrix_visitor_state<
					typename matrix_vismap<Arg1, setting_t>::type
				>::type arg1_state_t;

		typedef typename matrix_visitor_state<
					typename matrix_vismap<Arg2, setting_t>::type
				>::type arg2_state_t;

		typedef std::pair<arg1_state_t, arg2_state_t> type;
	};


	template<typename VisSetting, class Fun, class Arg>
	struct matrix_visitor_state<unary_ewise_mvisitor<VisSetting, Fun, Arg> >
	{
		typedef typename unary_ewise_percol_vstate<
				typename VisSetting::kernel_category,
				VisSetting::ct_rows, VisSetting::ct_cols, Arg>::type type;
	};

	template<typename VisSetting, class Fun, class Arg1, class Arg2>
	struct matrix_visitor_state<binary_ewise_mvisitor<VisSetting, Fun, Arg1, Arg2> >
	{
		typedef typename binary_ewise_percol_vstate<
				typename VisSetting::kernel_category,
				VisSetting::ct_rows, VisSetting::ct_cols, Arg1, Arg2>::type type;
	};


	/********************************************
	 *
	 *  Visitors
	 *
	 ********************************************/

	// unary

	template<int ME, int NE, class Fun, class Arg>
	class unary_ewise_mvisitor<
				matrix_visit_setting<linear_vis, scalar_kernel_t, ME, NE>, Fun, Arg>
	: public ILinearMatrixScalarVisitor<unary_ewise_mvisitor<
	  	  	  	matrix_visit_setting<linear_vis, scalar_kernel_t, ME, NE>, Fun, Arg>,
	  	  	  	typename Fun::result_type>
	{
	public:
		typedef matrix_visit_setting<linear_vis, scalar_kernel_t, ME, NE> setting_t;
		typedef typename matrix_vismap<Arg, setting_t>::type arg_visitor_t;

		template<typename Arg_HP>
		LMAT_ENSURE_INLINE
		unary_ewise_mvisitor(const unary_ewise_expr<Fun, Arg_HP, Arg>& expr)
		: m_fun(expr.fun()), m_arg_visitor(expr.arg())
		{
		}

		LMAT_ENSURE_INLINE
		typename Fun::result_type get_scalar(const index_t i) const
		{
			return m_fun(m_arg_visitor.get_scalar(i));
		}

	private:
		Fun m_fun;
		arg_visitor_t m_arg_visitor;
	};




	template<int ME, int NE, class Fun, class Arg>
	class unary_ewise_mvisitor<
				matrix_visit_setting<percol_vis, scalar_kernel_t, ME, NE>, Fun, Arg>
	: public IPerColMatrixScalarVisitor<unary_ewise_mvisitor<
	  	  	  	matrix_visit_setting<percol_vis, scalar_kernel_t, ME, NE>, Fun, Arg>,
	  	  	  	typename Fun::result_type>
	{
	public:
		typedef matrix_visit_setting<percol_vis, scalar_kernel_t, ME, NE> setting_t;
		typedef typename matrix_vismap<Arg, setting_t>::type arg_visitor_t;
		typedef typename unary_ewise_percol_vstate<scalar_kernel_t, ME, NE, Arg>::type state_t;

		template<typename Arg_HP>
		LMAT_ENSURE_INLINE
		unary_ewise_mvisitor(const unary_ewise_expr<Fun, Arg_HP, Arg>& expr)
		: m_fun(expr.fun()), m_arg_visitor(expr.arg())
		{
		}

		LMAT_ENSURE_INLINE
		typename Fun::result_type get_scalar(const index_t i, const state_t& s) const
		{
			return m_fun(m_arg_visitor.get_scalar(i, s));
		}

		LMAT_ENSURE_INLINE
		state_t col_state(const index_t j) const
		{
			return m_arg_visitor.col_state(j);
		}

	private:
		Fun m_fun;
		arg_visitor_t m_arg_visitor;
	};


	// binary

	template<int ME, int NE, class Fun, class Arg1, class Arg2>
	class binary_ewise_mvisitor<
				matrix_visit_setting<linear_vis, scalar_kernel_t, ME, NE>, Fun, Arg1, Arg2>
	: public ILinearMatrixScalarVisitor<binary_ewise_mvisitor<
	  	  	  	matrix_visit_setting<linear_vis, scalar_kernel_t, ME, NE>, Fun, Arg1, Arg2>,
	  	  	  	typename Fun::result_type>
	{
	public:
		typedef matrix_visit_setting<linear_vis, scalar_kernel_t, ME, NE> setting_t;
		typedef typename matrix_vismap<Arg1, setting_t>::type arg1_visitor_t;
		typedef typename matrix_vismap<Arg2, setting_t>::type arg2_visitor_t;

		template<typename Arg1_HP, typename Arg2_HP>
		LMAT_ENSURE_INLINE
		binary_ewise_mvisitor(const binary_ewise_expr<Fun, Arg1_HP, Arg1, Arg2_HP, Arg2>& expr)
		: m_fun(expr.fun()), m_arg1_visitor(expr.first_arg()), m_arg2_visitor(expr.second_arg())
		{
		}

		LMAT_ENSURE_INLINE
		typename Fun::result_type get_scalar(const index_t i) const
		{
			return m_fun(
					m_arg1_visitor.get_scalar(i),
					m_arg2_visitor.get_scalar(i));
		}

	private:
		Fun m_fun;
		arg1_visitor_t m_arg1_visitor;
		arg2_visitor_t m_arg2_visitor;
	};


	template<int ME, int NE, class Fun, class Arg1, class Arg2>
	class binary_ewise_mvisitor<
				matrix_visit_setting<percol_vis, scalar_kernel_t, ME, NE>, Fun, Arg1, Arg2>
	: public IPerColMatrixScalarVisitor<binary_ewise_mvisitor<
	  	  	  	matrix_visit_setting<percol_vis, scalar_kernel_t, ME, NE>, Fun, Arg1, Arg2>,
	  	  	  	typename Fun::result_type>
	{
	public:
		typedef matrix_visit_setting<percol_vis, scalar_kernel_t, ME, NE> setting_t;
		typedef typename matrix_vismap<Arg1, setting_t>::type arg1_visitor_t;
		typedef typename matrix_vismap<Arg2, setting_t>::type arg2_visitor_t;
		typedef typename binary_ewise_percol_vstate<scalar_kernel_t, ME, NE, Arg1, Arg2>::type state_t;

		template<typename Arg1_HP, typename Arg2_HP>
		LMAT_ENSURE_INLINE
		binary_ewise_mvisitor(const binary_ewise_expr<Fun, Arg1_HP, Arg1, Arg2_HP, Arg2>& expr)
		: m_fun(expr.fun()), m_arg1_visitor(expr.first_arg()), m_arg2_visitor(expr.second_arg())
		{
		}

		LMAT_ENSURE_INLINE
		typename Fun::result_type get_scalar(const index_t i, const state_t& s) const
		{
			return m_fun(
					m_arg1_visitor.get_scalar(i, s.first),
					m_arg2_visitor.get_scalar(i, s.second));
		}

		LMAT_ENSURE_INLINE
		state_t col_state(const index_t j) const
		{
			return state_t(m_arg1_visitor.col_state(j), m_arg2_visitor.col_state(j));
		}

	private:
		Fun m_fun;
		arg1_visitor_t m_arg1_visitor;
		arg2_visitor_t m_arg2_visitor;
	};



	/********************************************
	 *
	 *  Evaluator dispatch
	 *
	 ********************************************/

	template<typename Fun, typename Arg_HP, class Arg, typename Setting>
	struct matrix_vismap<
		unary_ewise_expr<Fun, Arg_HP, Arg>, Setting>
	{
		typedef unary_ewise_mvisitor<Setting, Fun, Arg> type;
		static const int cost = matrix_vismap<Arg, Setting>::cost;
	};


	template<typename Fun, typename Arg1_HP, class Arg1, typename Arg2_HP, class Arg2, typename Setting>
	struct matrix_vismap<
		binary_ewise_expr<Fun, Arg1_HP, Arg1, Arg2_HP, Arg2>, Setting>
	{
		typedef binary_ewise_mvisitor<Setting, Fun, Arg1, Arg2> type;

		static const int cost =
				matrix_vismap<Arg1, Setting>::cost +
				matrix_vismap<Arg2, Setting>::cost;
	};


	template<typename Fun, typename Arg_HP, class Arg, class Dst>
	struct default_matrix_eval_policy<unary_ewise_expr<Fun, Arg_HP, Arg>, Dst>
	{
		typedef unary_ewise_expr<Fun, Arg_HP, Arg> expr_t;
		typedef typename default_matrix_visit_policy<expr_t, Dst>::type type;
	};

	template<typename Fun, typename Arg1_HP, class Arg1, typename Arg2_HP, class Arg2, class Dst>
	struct default_matrix_eval_policy<binary_ewise_expr<Fun, Arg1_HP, Arg1, Arg2_HP, Arg2>, Dst>
	{
		typedef binary_ewise_expr<Fun, Arg1_HP, Arg1, Arg2_HP, Arg2> expr_t;
		typedef typename default_matrix_visit_policy<expr_t, Dst>::type type;
	};

}

#endif /* MATRIX_EWISE_EVAL_H_ */
