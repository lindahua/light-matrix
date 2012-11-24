/*
 * @file matrix_ewise_expr.h
 *
 * Generic element-wise matrix expression
 *
 * @author Dahua Lin
 */

#ifndef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_EWISE_EXPR_H_
#define LIGHTMAT_MATRIX_EWISE_EXPR_H_

#include <light_mat/matexpr/matexpr_fwd.h>

namespace lmat
{

	/********************************************
	 *
	 *  Expression classes
	 *
	 ********************************************/

	template<class Fun, typename Arg_HP, class Arg>
	struct matrix_traits<unary_ewise_expr<Fun, Arg_HP, Arg> >
	{
		static const int num_dimensions = 2;
		static const int ct_num_rows = ct_rows<Arg>::value;
		static const int ct_num_cols = ct_cols<Arg>::value;

		static const bool is_readonly = true;

		typedef typename Fun::result_type value_type;
		typedef typename matrix_traits<Arg>::domain domain;
	};

	template<class Fun, typename Arg1_HP, class Arg1, class Arg2_HP, class Arg2>
	struct matrix_traits<binary_ewise_expr<Fun, Arg1_HP, Arg1, Arg2_HP, Arg2> >
	{
#ifdef LMAT_USE_STATIC_ASSERT
		static_assert(has_same_domain<Arg1, Arg2>::value,
				"Arg1 and Arg2 should be in the same domain");
#endif

		static const int num_dimensions = 2;
		static const int ct_num_rows = binary_ct_rows<Arg1, Arg2>::value;
		static const int ct_num_cols = binary_ct_cols<Arg1, Arg2>::value;

		static const bool is_readonly = true;

		typedef typename Fun::result_type value_type;
		typedef typename matrix_traits<Arg1>::domain domain;
	};


	template<class Fun, typename Arg_HP, class Arg>
	class unary_ewise_expr
	: public unary_expr_base<Arg_HP, Arg>
	, public IMatrixXpr<unary_ewise_expr<Fun, Arg_HP, Arg>, typename Fun::result_type>
	{
#ifdef LMAT_USE_STATIC_ASSERT
		static_assert(is_unary_ewise_functor<Fun>::value, "Fun must be a unary_ewise_functor");
		static_assert(is_mat_xpr<Arg>::value, "Arg must be a matrix expression class.");
#endif

	public:
		typedef unary_expr_base<Arg_HP, Arg> base_t;
		typedef typename Fun::result_type value_type;

		LMAT_ENSURE_INLINE
		unary_ewise_expr(const Fun& fun,
				const arg_forwarder<Arg_HP, Arg>& arg_fwd)
		: base_t(arg_fwd), m_fun(fun) { }

		LMAT_ENSURE_INLINE const Fun& fun() const
		{
			return m_fun;
		}

		LMAT_ENSURE_INLINE index_t nelems() const
		{
			return this->arg().nelems();
		}

		LMAT_ENSURE_INLINE index_t nrows() const
		{
			return this->arg().nrows();
		}

		LMAT_ENSURE_INLINE index_t ncolumns() const
		{
			return this->arg().ncolumns();
		}

	private:
		Fun m_fun;
	};


	template<class Fun, typename Arg1_HP, class Arg1, typename Arg2_HP, class Arg2>
	class binary_ewise_expr
	: public binary_expr_base<Arg1_HP, Arg1, Arg2_HP, Arg2>
	, public IMatrixXpr<binary_ewise_expr<Fun, Arg1_HP, Arg1, Arg2_HP, Arg2>, typename Fun::result_type>
	{
#ifdef LMAT_USE_STATIC_ASSERT
		static_assert(is_binary_ewise_functor<Fun>::value, "Fun must be a binary_ewise_functor");
		static_assert(is_mat_xpr<Arg1>::value, "Arg1 must be a matrix expression class.");
		static_assert(is_mat_xpr<Arg2>::value, "Arg2 must be a matrix expression class.");
#endif

	public:
		typedef binary_expr_base<Arg1_HP, Arg1, Arg2_HP, Arg2> base_t;
		typedef typename Fun::result_type value_type;

		LMAT_ENSURE_INLINE
		binary_ewise_expr(const Fun& fun,
				const arg_forwarder<Arg1_HP, Arg1>& arg1,
				const arg_forwarder<Arg2_HP, Arg2>& arg2)
		: base_t(arg1, arg2), m_fun(fun)
		{
			LMAT_CHECK_SAME_SHAPE(this->first_arg(), this->second_arg())
		}

		LMAT_ENSURE_INLINE const Fun& fun() const
		{
			return m_fun;
		}

		LMAT_ENSURE_INLINE index_t nelems() const
		{
			return this->first_arg().nelems();
		}

		LMAT_ENSURE_INLINE index_t nrows() const
		{
			return this->first_arg().nrows();
		}

		LMAT_ENSURE_INLINE index_t ncolumns() const
		{
			return this->first_arg().ncolumns();
		}

	private:
		Fun m_fun;
	};


	/********************************************
	 *
	 *  Expression mapping and construction
	 *
	 ********************************************/

	template<class Fun>
	struct ewise_t
	{
		const Fun& fun;

		LMAT_ENSURE_INLINE
		ewise_t(const Fun& f) : fun(f) { }
	};

	template<class Fun>
	LMAT_ENSURE_INLINE
	inline ewise_t<Fun> ewise(const Fun& f)
	{
		return ewise_t<Fun>(f);
	}

	template<class Fun, class Arg>
	struct unary_expr_verifier<ewise_t<Fun>, Arg>
	{
		static const bool value = is_mat_xpr<Arg>::value;
	};

	template<class Fun, class Arg1, class Arg2>
	struct binary_expr_verifier<ewise_t<Fun>, Arg1, Arg2>
	{
		static const bool value = is_mat_xpr<Arg1>::value && is_mat_xpr<Arg2>::value;
	};

	template<class Fun, typename Arg_HP, class Arg>
	struct unary_expr_map<ewise_t<Fun>, Arg_HP, Arg>
	{
		typedef unary_ewise_expr<Fun, Arg_HP, Arg> type;

		LMAT_ENSURE_INLINE
		static type get(const ewise_t<Fun>& spec,
				const arg_forwarder<Arg_HP, Arg>& arg_fwd)
		{
			return type(spec.fun, arg_fwd);
		}
	};

	template<class Fun, typename Arg1_HP, class Arg1, typename Arg2_HP, class Arg2>
	struct binary_expr_map<ewise_t<Fun>, Arg1_HP, Arg1, Arg2_HP, Arg2>
	{
		typedef binary_ewise_expr<Fun, Arg1_HP, Arg1, Arg2_HP, Arg2> type;

		LMAT_ENSURE_INLINE
		static type get(const ewise_t<Fun>& spec,
				const arg_forwarder<Arg1_HP, Arg1>& arg1_fwd,
				const arg_forwarder<Arg2_HP, Arg2>& arg2_fwd)
		{
			return type(spec.fun, arg1_fwd, arg2_fwd);
		}
	};




	template<class Fun, class Arg>
	LMAT_ENSURE_INLINE
	inline typename unary_expr_map<ewise_t<Fun>,
		ref_arg_t, Arg
	>::type
	ewise(const Fun& fun, const IMatrixXpr<Arg, typename Fun::arg_type>& arg)
	{
		return make_expr(ewise(fun), ref_arg(arg.derived()) );
	}

	template<class Fun, class Arg1, class Arg2>
	LMAT_ENSURE_INLINE
	inline typename binary_expr_map<ewise_t<Fun>,
		ref_arg_t, Arg1,
		ref_arg_t, Arg2
	>::type
	ewise(  const Fun& fun,
			const IMatrixXpr<Arg1, typename Fun::first_arg_type>& arg1,
			const IMatrixXpr<Arg2, typename Fun::second_arg_type>& arg2 )
	{
		return make_expr(ewise(fun),
				ref_arg(arg1.derived()),
				ref_arg(arg2.derived()));
	}


}

#endif 




