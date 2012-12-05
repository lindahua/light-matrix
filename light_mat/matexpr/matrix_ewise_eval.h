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

	template<typename Acc, typename Ker, typename Tag, class ArgList>
	class ewise_accessor;


	/********************************************
	 *
	 *  Evaluation scheme maps
	 *
	 ********************************************/

	// default evaluate scheme

	template<typename Tag, class QLst, class Dst>
	LMAT_ENSURE_INLINE
	inline void evaluate(const ewise_expr<Tag, QLst>& expr,
			IDenseMatrix<Dst, typename matrix_traits<ewise_expr<Tag, QLst> >::value_type >& dst)
	{
		typedef typename default_macc_policy<ewise_expr<Tag, QLst>, Dst>::type policy_t;
		evaluate(expr, dst.derived(), policy_t());
	}


	// accessors

	template<typename Tag, class QArgList, typename Acc, typename Ker>
	struct macc_accessor_map<ewise_expr<Tag, QArgList>, macc_policy<Acc, Ker> >
	{
		typedef typename meta::transform_<meta::argument_of, QArgList>::type arg_list;
		typedef ewise_accessor<Acc, Ker, Tag, arg_list> type;
	};


	/********************************************
	 *
	 *  cost model
	 *
	 ********************************************/

	template<typename Tag, class QArg1, typename Acc, typename Ker>
	struct macc_cost< ewise_expr<Tag, LMAT_TYPELIST_1(QArg1) >, Acc, Ker>
	{
		static const int value = macc_cost<typename QArg1::argument, Acc, Ker>::value;
	};

	template<typename Tag, class QArg1, class QArg2, typename Acc, typename Ker>
	struct macc_cost< ewise_expr<Tag, LMAT_TYPELIST_2(QArg1, QArg2) >, Acc, Ker>
	{
		static const int value =
				macc_cost<typename QArg1::argument, Acc, Ker>::value +
				macc_cost<typename QArg2::argument, Acc, Ker>::value;
	};


	/********************************************
	 *
	 *  helpers
	 *
	 ********************************************/

	template<typename Acc, typename Ker, typename Tag, class ArgList>
	struct ewise_accbase
	{
		typedef typename meta::transform_<meta::value_type_of, ArgList>::type value_type_list;
		typedef typename op_result<Tag, value_type_list>::type value_type;

		typedef typename
				meta::if_<is_same<Acc, linear_macc>,
					ILinearMatrixScalarAccessor<ewise_accessor<Acc, Ker, Tag, ArgList>, value_type>,
					IPerColMatrixScalarAccessor<ewise_accessor<Acc, Ker, Tag, ArgList>, value_type>
				>::type base_type;
	};


	/********************************************
	 *
	 *  Linear scalar accessors
	 *
	 ********************************************/

	template<typename Tag, class Arg>
	class ewise_accessor<linear_macc, scalar_ker, Tag, LMAT_TYPELIST_1(Arg) >
	: public ewise_accbase<linear_macc, scalar_ker, Tag, LMAT_TYPELIST_1(Arg) >::base_type
	{
		typedef typename matrix_traits<Arg>::value_type arg_value_type;
		typedef LMAT_TYPELIST_1(arg_value_type) arg_value_types;

		typedef typename op_result<Tag, arg_value_types>::type value_type;
		typedef typename op_fun<Tag, scalar_ker, arg_value_types>::type fun_t;

		typedef macc_policy<linear_macc, scalar_ker> policy_t;
		typedef typename macc_accessor_map<Arg, policy_t>::type arg_accessor_t;

	public:
		template<class QArg>
		LMAT_ENSURE_INLINE
		ewise_accessor(const ewise_expr<Tag, LMAT_TYPELIST_1(QArg) >& expr)
		: m_fun(expr.tag()), m_arg_acc(expr.args().arg1()) { }

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
	class ewise_accessor<linear_macc, scalar_ker, Tag, LMAT_TYPELIST_2(Arg1, Arg2) >
	: public ewise_accbase<linear_macc, scalar_ker, Tag, LMAT_TYPELIST_2(Arg1, Arg2) >::base_type
	{
		typedef typename matrix_traits<Arg1>::value_type arg1_value_type;
		typedef typename matrix_traits<Arg2>::value_type arg2_value_type;
		typedef LMAT_TYPELIST_2(arg1_value_type, arg2_value_type) arg_value_types;

		typedef typename op_result<Tag, arg_value_types>::type value_type;
		typedef typename op_fun<Tag, scalar_ker, arg_value_types>::type fun_t;

		typedef macc_policy<linear_macc, scalar_ker> policy_t;
		typedef typename macc_accessor_map<Arg1, policy_t>::type arg1_accessor_t;
		typedef typename macc_accessor_map<Arg2, policy_t>::type arg2_accessor_t;

	public:
		template<class QArg1, class QArg2>
		LMAT_ENSURE_INLINE
		ewise_accessor(const ewise_expr<Tag, LMAT_TYPELIST_2(QArg1, QArg2) >& expr)
		: m_fun(expr.tag())
		, m_arg1_acc(expr.args().arg1())
		, m_arg2_acc(expr.args().arg2()) { }

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
		ewise_accessor<percol_macc, scalar_ker, Tag, LMAT_TYPELIST_1(Arg) > >
	{
		typedef typename macc_accessor_map<Arg, macc_policy<percol_macc, scalar_ker> >::type arg_accessor_t;
		typedef typename percol_macc_state_map<arg_accessor_t>::type type;
	};

	template<typename Tag, class Arg>
	class ewise_accessor<percol_macc, scalar_ker, Tag, LMAT_TYPELIST_1(Arg) >
	: public ewise_accbase<percol_macc, scalar_ker, Tag, LMAT_TYPELIST_1(Arg) >::base_type
	{
		typedef typename matrix_traits<Arg>::value_type arg_value_type;
		typedef LMAT_TYPELIST_1(arg_value_type) arg_value_types;

		typedef typename op_result<Tag, arg_value_types>::type value_type;
		typedef typename op_fun<Tag, scalar_ker, arg_value_types>::type fun_t;

		typedef macc_policy<percol_macc, scalar_ker> policy_t;
		typedef typename macc_accessor_map<Arg, policy_t>::type arg_accessor_t;

		typedef typename percol_macc_state_map<arg_accessor_t>::type col_state_t;

	public:

		template<class QArg>
		LMAT_ENSURE_INLINE
		ewise_accessor(const ewise_expr<Tag, LMAT_TYPELIST_1(QArg)>& expr)
		: m_fun(expr.tag()), m_arg_acc(expr.args().arg1()) { }

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
		ewise_accessor<percol_macc, scalar_ker, Tag, LMAT_TYPELIST_2(Arg1, Arg2) > >
	{
		typedef macc_policy<percol_macc, scalar_ker> policy_t;
		typedef typename macc_accessor_map<Arg1, policy_t>::type arg1_accessor_t;
		typedef typename macc_accessor_map<Arg2, policy_t>::type arg2_accessor_t;
		typedef typename std::pair<
				typename percol_macc_state_map<arg1_accessor_t>::type,
				typename percol_macc_state_map<arg2_accessor_t>::type> type;
	};


	template<typename Tag, class Arg1, class Arg2>
	class ewise_accessor<percol_macc, scalar_ker, Tag, LMAT_TYPELIST_2(Arg1, Arg2) >
	: public ewise_accbase<percol_macc, scalar_ker, Tag, LMAT_TYPELIST_2(Arg1, Arg2) >::base_type
	{
		typedef typename matrix_traits<Arg1>::value_type arg1_value_type;
		typedef typename matrix_traits<Arg2>::value_type arg2_value_type;
		typedef typename LMAT_TYPELIST_2( arg1_value_type, arg2_value_type ) arg_value_types;

		typedef typename op_result<Tag, arg_value_types >::type value_type;
		typedef typename op_fun<Tag, scalar_ker, arg_value_types >::type fun_t;

		typedef macc_policy<percol_macc, scalar_ker> policy_t;
		typedef typename macc_accessor_map<Arg1, policy_t>::type arg1_accessor_t;
		typedef typename macc_accessor_map<Arg2, policy_t>::type arg2_accessor_t;

		typedef typename percol_macc_state_map<arg1_accessor_t>::type arg1_col_state_t;
		typedef typename percol_macc_state_map<arg2_accessor_t>::type arg2_col_state_t;
		typedef typename std::pair<arg1_col_state_t, arg2_col_state_t> col_state_t;

	public:

		template<class QArg1, class QArg2>
		ewise_accessor(const ewise_expr<Tag, LMAT_TYPELIST_2(QArg1, QArg2) >& expr)
		: m_fun(expr.tag()), m_arg1_acc(expr.args().arg1()), m_arg2_acc(expr.args().arg2()) { }

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



