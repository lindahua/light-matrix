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

#include <light_mat/matexpr/matrix_ewise_expr.h>
#include <light_mat/matexpr/dense_accessors.h>
#include <utility>

namespace lmat
{

	// forward declaration

	template<typename Acc, typename Ker, typename Tag, class Arg>
	class unary_ewise_accessor;

	template<typename Acc, typename Ker, typename Tag, class Arg1, class Arg2>
	class binary_ewise_accessor;


	/********************************************
	 *
	 *  Evaluation scheme maps
	 *
	 ********************************************/

	// default evaluate scheme

	template<typename Tag, typename Arg_HP, class Arg, class DMat, typename DT>
	LMAT_ENSURE_INLINE
	inline typename default_macc_scheme<unary_ewise_expr<Tag, Arg_HP, Arg>, DMat>::type
	get_default_eval_scheme(
			const unary_ewise_expr<Tag, Arg_HP, Arg>& sexpr,
			IDenseMatrix<DMat, DT>& dmat)
	{
		typedef unary_ewise_expr<Tag, Arg_HP, Arg> expr_t;
		return default_macc_scheme<expr_t, DMat>::get(sexpr, dmat.derived());
	}

	template<typename Tag, typename Arg1_HP, class Arg1, typename Arg2_HP, class Arg2, class DMat, typename DT>
	LMAT_ENSURE_INLINE
	inline typename default_macc_scheme<binary_ewise_expr<Tag, Arg1_HP, Arg1, Arg2_HP, Arg2>, DMat>::type
	get_default_eval_scheme(
			const binary_ewise_expr<Tag, Arg1_HP, Arg1, Arg2_HP, Arg2>& sexpr,
			IDenseMatrix<DMat, DT>& dmat)
	{
		typedef binary_ewise_expr<Tag, Arg1_HP, Arg1, Arg2_HP, Arg2> expr_t;
		return default_macc_scheme<expr_t, DMat>::get(sexpr, dmat.derived());
	}


	// accessors

	template<typename Tag, typename Arg_HP, class Arg, typename Acc, typename Ker>
	struct macc_accessor_map<unary_ewise_expr<Tag, Arg_HP, Arg>, Acc, Ker>
	{
		typedef unary_ewise_accessor<Acc, Ker, Tag, Arg> type;
	};

	template<typename Tag, typename Arg1_HP, class Arg1, typename Arg2_HP, class Arg2, typename Acc, typename Ker>
	struct macc_accessor_map<binary_ewise_expr<Tag, Arg1_HP, Arg1, Arg2_HP, Arg2>, Acc, Ker>
	{
		typedef binary_ewise_accessor<Acc, Ker, Tag, Arg1, Arg2> type;
	};


	/********************************************
	 *
	 *  cost model
	 *
	 ********************************************/

	template<typename Tag, typename Arg_HP, class Arg,
		typename AccCate, typename KerCate>
	struct macc_cost<unary_ewise_expr<Tag, Arg_HP, Arg>, AccCate, KerCate>
	{
		static const int value = macc_cost<Arg, AccCate, KerCate>::value;
	};

	template<typename Tag, typename Arg1_HP, class Arg1, typename Arg2_HP, class Arg2,
		typename AccCate, typename KerCate>
	struct macc_cost<binary_ewise_expr<Tag, Arg1_HP, Arg1, Arg2_HP, Arg2>, AccCate, KerCate>
	{
		static const int value =
				macc_cost<Arg1, AccCate, KerCate>::value +
				macc_cost<Arg2, AccCate, KerCate>::value;
	};


	/********************************************
	 *
	 *  helpers
	 *
	 ********************************************/

	template<typename Acc, typename Ker, typename Tag, class Arg>
	struct unary_ewise_accbase
	{
		typedef typename matrix_traits<Arg>::value_type arg_value_type;
		typedef typename unary_op_result<Tag, arg_value_type>::type value_type;

		typedef typename
				if_<is_same<Acc, linear_macc>,
					ILinearMatrixScalarAccessor<unary_ewise_accessor<Acc, Ker, Tag, Arg>, value_type>,
					IPerColMatrixScalarAccessor<unary_ewise_accessor<Acc, Ker, Tag, Arg>, value_type>
				>::type base_type;
	};

	template<typename Acc, typename Ker, typename Tag, class Arg1, class Arg2>
	struct binary_ewise_accbase
	{
		typedef typename matrix_traits<Arg1>::value_type arg1_value_type;
		typedef typename matrix_traits<Arg2>::value_type arg2_value_type;
		typedef typename binary_op_result<Tag, arg1_value_type, arg2_value_type>::type value_type;

		typedef typename
				if_<is_same<Acc, linear_macc>,
					ILinearMatrixScalarAccessor<binary_ewise_accessor<Acc, Ker, Tag, Arg1, Arg2>, value_type>,
					IPerColMatrixScalarAccessor<binary_ewise_accessor<Acc, Ker, Tag, Arg1, Arg2>, value_type>
				>::type base_type;
	};

	/********************************************
	 *
	 *  Linear scalar accessors
	 *
	 ********************************************/

	template<typename Tag, class Arg>
	class unary_ewise_accessor<linear_macc, scalar_ker, Tag, Arg>
	: public unary_ewise_accbase<linear_macc, scalar_ker, Tag, Arg>::base_type
	{
		typedef typename matrix_traits<Arg>::value_type arg_value_type;
		typedef typename unary_op_result<Tag, arg_value_type>::type value_type;

		typedef typename unary_op_fun<Tag, scalar_ker, arg_value_type>::type fun_t;
		typedef typename macc_accessor_map<Arg, linear_macc, scalar_ker>::type arg_accessor_t;

	public:

		template<typename Arg_HP>
		unary_ewise_accessor(const unary_ewise_expr<Tag, Arg_HP, Arg>& expr)
		: m_fun(expr.tag()), m_arg_acc(expr.arg()) { }

		LMAT_ENSURE_INLINE
		value_type get_scalar(const index_t i) const
		{
			return m_fun(m_arg_acc.get_scalar(i));
		}

	private:
		fun_t m_fun;
		arg_accessor_t m_arg_acc;
	};


	template<typename Tag, class Arg1, class Arg2>
	class binary_ewise_accessor<linear_macc, scalar_ker, Tag, Arg1, Arg2>
	: public binary_ewise_accbase<linear_macc, scalar_ker, Tag, Arg1, Arg2>::base_type
	{
		typedef typename matrix_traits<Arg1>::value_type arg1_value_type;
		typedef typename matrix_traits<Arg2>::value_type arg2_value_type;
		typedef typename binary_op_result<Tag, arg1_value_type, arg2_value_type>::type value_type;

		typedef typename binary_op_fun<Tag, scalar_ker, arg1_value_type, arg2_value_type>::type fun_t;
		typedef typename macc_accessor_map<Arg1, linear_macc, scalar_ker>::type arg1_accessor_t;
		typedef typename macc_accessor_map<Arg2, linear_macc, scalar_ker>::type arg2_accessor_t;

	public:

		template<typename Arg1_HP, typename Arg2_HP>
		binary_ewise_accessor(const binary_ewise_expr<Tag, Arg1_HP, Arg1, Arg2_HP, Arg2>& expr)
		: m_fun(expr.tag()), m_arg1_acc(expr.first_arg()), m_arg2_acc(expr.second_arg()) { }

		LMAT_ENSURE_INLINE
		value_type get_scalar(const index_t i) const
		{
			return m_fun(m_arg1_acc.get_scalar(i), m_arg2_acc.get_scalar(i));
		}

	private:
		fun_t m_fun;
		arg1_accessor_t m_arg1_acc;
		arg2_accessor_t m_arg2_acc;
	};


	/********************************************
	 *
	 *  Percol scalar accessors
	 *
	 ********************************************/

	template<typename Tag, class Arg>
	struct percol_macc_state_map<
		unary_ewise_accessor<percol_macc, scalar_ker, Tag, Arg> >
	{
		typedef typename macc_accessor_map<Arg, percol_macc, scalar_ker>::type arg_accessor_t;
		typedef typename percol_macc_state_map<arg_accessor_t>::type type;
	};

	template<typename Tag, class Arg>
	class unary_ewise_accessor<percol_macc, scalar_ker, Tag, Arg>
	: public unary_ewise_accbase<percol_macc, scalar_ker, Tag, Arg>::base_type
	{
		typedef typename matrix_traits<Arg>::value_type arg_value_type;
		typedef typename unary_op_result<Tag, arg_value_type>::type value_type;

		typedef typename unary_op_fun<Tag, scalar_ker, arg_value_type>::type fun_t;
		typedef typename macc_accessor_map<Arg, percol_macc, scalar_ker>::type arg_accessor_t;

		typedef typename percol_macc_state_map<arg_accessor_t>::type col_state_t;

	public:

		template<typename Arg_HP>
		unary_ewise_accessor(const unary_ewise_expr<Tag, Arg_HP, Arg>& expr)
		: m_fun(expr.tag()), m_arg_acc(expr.arg()) { }

		LMAT_ENSURE_INLINE
		value_type get_scalar(const index_t i, const col_state_t& s) const
		{
			return m_fun(m_arg_acc.get_scalar(i, s));
		}

		LMAT_ENSURE_INLINE
		col_state_t col_state(const index_t j) const
		{
			return m_arg_acc.col_state(j);
		}

	private:
		fun_t m_fun;
		arg_accessor_t m_arg_acc;
	};


	template<typename Tag, class Arg1, class Arg2>
	struct percol_macc_state_map<
		binary_ewise_accessor<percol_macc, scalar_ker, Tag, Arg1, Arg2> >
	{
		typedef typename macc_accessor_map<Arg1, percol_macc, scalar_ker>::type arg1_accessor_t;
		typedef typename macc_accessor_map<Arg2, percol_macc, scalar_ker>::type arg2_accessor_t;
		typedef typename std::pair<
				typename percol_macc_state_map<arg1_accessor_t>::type,
				typename percol_macc_state_map<arg2_accessor_t>::type> type;
	};


	template<typename Tag, class Arg1, class Arg2>
	class binary_ewise_accessor<percol_macc, scalar_ker, Tag, Arg1, Arg2>
	: public binary_ewise_accbase<percol_macc, scalar_ker, Tag, Arg1, Arg2>::base_type
	{
		typedef typename matrix_traits<Arg1>::value_type arg1_value_type;
		typedef typename matrix_traits<Arg2>::value_type arg2_value_type;
		typedef typename binary_op_result<Tag, arg1_value_type, arg2_value_type>::type value_type;

		typedef typename binary_op_fun<Tag, scalar_ker, arg1_value_type, arg2_value_type>::type fun_t;
		typedef typename macc_accessor_map<Arg1, percol_macc, scalar_ker>::type arg1_accessor_t;
		typedef typename macc_accessor_map<Arg2, percol_macc, scalar_ker>::type arg2_accessor_t;

		typedef typename percol_macc_state_map<arg1_accessor_t>::type arg1_col_state_t;
		typedef typename percol_macc_state_map<arg2_accessor_t>::type arg2_col_state_t;
		typedef typename std::pair<arg1_col_state_t, arg2_col_state_t> col_state_t;

	public:

		template<typename Arg1_HP, typename Arg2_HP>
		binary_ewise_accessor(const binary_ewise_expr<Tag, Arg1_HP, Arg1, Arg2_HP, Arg2>& expr)
		: m_fun(expr.tag()), m_arg1_acc(expr.first_arg()), m_arg2_acc(expr.second_arg()) { }

		LMAT_ENSURE_INLINE
		value_type get_scalar(const index_t i, const col_state_t& s) const
		{
			return m_fun(m_arg1_acc.get_scalar(i, s.first), m_arg2_acc.get_scalar(i, s.second));
		}

		LMAT_ENSURE_INLINE
		col_state_t col_state(const index_t j) const
		{
			return col_state_t(m_arg1_acc.col_state(j), m_arg2_acc.col_state(j));
		}

	private:
		fun_t m_fun;
		arg1_accessor_t m_arg1_acc;
		arg2_accessor_t m_arg2_acc;
	};

}

#endif



