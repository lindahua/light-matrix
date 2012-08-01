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

#include <light_mat/matrix/matrix_expr_base.h>

namespace lmat
{

	/********************************************
	 *
	 *  Auxiliary devices
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


	/********************************************
	 *
	 *  Expression classes
	 *
	 ********************************************/

	template<class Fun, class Arg_Holder>
	struct matrix_traits<unary_ewise_expr<Fun, Arg_Holder> >
	{
		typedef typename Arg_Holder::arg_type arg_type;

		static const int num_dimensions = 2;
		static const int compile_time_num_rows = ct_rows<arg_type>::value;
		static const int compile_time_num_cols = ct_cols<arg_type>::value;

		static const bool is_readonly = true;

		typedef typename Fun::result_type value_type;
		typedef typename matrix_traits<arg_type>::domain domain;
	};

	template<class Fun, class Arg1_Holder, class Arg2_Holder>
	struct matrix_traits<binary_ewise_expr<Fun, Arg1_Holder, Arg2_Holder> >
	{
		typedef typename Arg1_Holder::arg_type arg1_type;
		typedef typename Arg2_Holder::arg_type arg2_type;

#ifdef LMAT_USE_STATIC_ASSERT
		static_assert(has_same_domain<arg1_type, arg2_type>::value,
				"Arg1 and Arg2 should be in the same domain");
#endif

		static const int num_dimensions = 2;
		static const int compile_time_num_rows = binary_ct_rows<arg1_type, arg2_type>::value;
		static const int compile_time_num_cols = binary_ct_cols<arg1_type, arg2_type>::value;

		static const bool is_readonly = true;

		typedef typename Fun::result_type value_type;
		typedef typename matrix_traits<arg1_type>::domain domain;
	};


	template<class Fun, class Arg_Holder>
	class unary_ewise_expr
	: public unary_expr<Arg_Holder>
	, public IMatrixXpr<unary_ewise_expr<Fun, Arg_Holder>, typename Fun::result_type>
	{
		typedef unary_expr<Arg_Holder> base_t;

#ifdef LMAT_USE_STATIC_ASSERT
		static_assert(is_unary_ewise_functor<Fun>::value, "Fun must be a unary_ewise_functor");
		static_assert(is_mat_xpr<typename base_t::arg_type>::value, "Arg must be a matrix expression class.");
#endif

	public:
		typedef typename Fun::result_type value_type;

		LMAT_ENSURE_INLINE
		unary_ewise_expr(const Fun& fun,
				const typename base_t::arg_forwarder& arg_fwd)
		: base_t(arg_fwd), m_fun(fun) { }

		LMAT_ENSURE_INLINE const Fun& fun() const
		{
			return m_fun;
		}

		LMAT_ENSURE_INLINE index_t nelems() const
		{
			return this->arg().nelems();
		}

		LMAT_ENSURE_INLINE size_t size() const
		{
			return this->arg().size();
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


	template<class Fun, class Arg1_Holder, class Arg2_Holder>
	class binary_ewise_expr
	: public binary_expr<Arg1_Holder, Arg2_Holder>
	, public IMatrixXpr<binary_ewise_expr<Fun, Arg1_Holder, Arg2_Holder>, typename Fun::result_type>
	{
		typedef binary_expr<Arg1_Holder, Arg2_Holder> base_t;

#ifdef LMAT_USE_STATIC_ASSERT
		static_assert(is_binary_ewise_functor<Fun>::value, "Fun must be a binary_ewise_functor");
		static_assert(is_mat_xpr<typename base_t::arg1_type>::value, "Arg1 must be a matrix expression class.");
		static_assert(is_mat_xpr<typename base_t::arg2_type>::value, "Arg2 must be a matrix expression class.");
#endif

	public:
		typedef typename Fun::result_type value_type;

		LMAT_ENSURE_INLINE
		binary_ewise_expr(const Fun& fun,
				const typename base_t::arg1_forwarder& arg1,
				const typename base_t::arg2_forwarder& arg2)
		: base_t(arg1, arg2), m_fun(fun)
		{
			check_same_size(
					this->first_arg(),
					this->second_arg(), "arg1 and arg2 must be of the same size.");
		}

		LMAT_ENSURE_INLINE const Fun& fun() const
		{
			return m_fun;
		}

		LMAT_ENSURE_INLINE index_t nelems() const
		{
			return this->first_arg().nelems();
		}

		LMAT_ENSURE_INLINE size_t size() const
		{
			return this->first_arg().size();
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

	template<class Fun, class Arg_Holder>
	struct unary_expr_map<ewise_t<Fun>, Arg_Holder>
	{
		typedef unary_ewise_expr<Fun, Arg_Holder> type;

		LMAT_ENSURE_INLINE
		static type get(const ewise_t<Fun>& spec,
				const typename arg_forwarder<Arg_Holder>::type& arg_fwd)
		{
			return type(spec.fun, arg_fwd);
		}
	};

	template<class Fun, class Arg1_Holder, class Arg2_Holder>
	struct binary_expr_map<ewise_t<Fun>, Arg1_Holder, Arg2_Holder>
	{
		typedef binary_ewise_expr<Fun, Arg1_Holder, Arg2_Holder> type;

		LMAT_ENSURE_INLINE
		static type get(const ewise_t<Fun>& spec,
				const typename arg_forwarder<Arg1_Holder>::type& arg1_fwd,
				const typename arg_forwarder<Arg2_Holder>::type& arg2_fwd)
		{
			return type(spec.fun, arg1_fwd, arg2_fwd);
		}
	};

	template<class Fun, class Arg1_Holder>
	struct binary_fix2_ewise_expr_map
	{
		typedef typename Fun::second_arg_type T2;

		typedef typename Arg1_Holder::arg_type Arg1;
		typedef const_matrix<T2, ct_rows<Arg1>::value, ct_cols<Arg1>::value> Arg2;

		typedef typename binary_expr_map<ewise_t<Fun>,
				Arg1_Holder, copy_arg_holder<Arg2> >::type type;
	};

	template<class Fun, class Arg2_Holder>
	struct binary_fix1_ewise_expr_map
	{
		typedef typename Fun::first_arg_type T1;

		typedef typename Arg2_Holder::arg_type Arg2;
		typedef const_matrix<T1, ct_rows<Arg2>::value, ct_cols<Arg2>::value> Arg1;

		typedef typename binary_expr_map<ewise_t<Fun>,
				copy_arg_holder<Arg1>, Arg2_Holder>::type type;
	};


	template<class Fun, class Arg>
	LMAT_ENSURE_INLINE
	inline typename unary_expr_map<ewise_t<Fun>,
		ref_arg_holder<Arg>
	>::type
	ewise(const Fun& fun, const IMatrixXpr<Arg, typename Fun::arg_type>& arg)
	{
		return make_expr(ewise(fun), ref_arg(arg.derived()) );
	}

	template<class Fun, class Arg1, class Arg2>
	LMAT_ENSURE_INLINE
	inline typename binary_expr_map<ewise_t<Fun>,
		ref_arg_holder<Arg1>,
		ref_arg_holder<Arg2>
	>::type
	ewise(  const Fun& fun,
			const IMatrixXpr<Arg1, typename Fun::first_arg_type>& arg1,
			const IMatrixXpr<Arg2, typename Fun::second_arg_type>& arg2 )
	{
		return make_expr(ewise(fun),
				ref_arg(arg1.derived()),
				ref_arg(arg2.derived()));
	}


	template<class Fun, class Arg1>
	LMAT_ENSURE_INLINE
	inline typename binary_fix2_ewise_expr_map<Fun, ref_arg_holder<Arg1> >::type
	ewise(  const Fun& fun,
			const IMatrixXpr<Arg1, typename Fun::first_arg_type>& arg1,
			const typename Fun::second_arg_type& arg2v )
	{
		const_matrix<typename Fun::second_arg_type,
			ct_rows<Arg1>::value,
			ct_cols<Arg1>::value>
		arg2(arg1.nrows(), arg1.ncolumns(), arg2v);

		return make_expr(ewise(fun),
				ref_arg(arg1.derived()),
				copy_arg(arg2) );
	}


	template<class Fun, class Arg2>
	LMAT_ENSURE_INLINE
	inline typename binary_fix1_ewise_expr_map<Fun, ref_arg_holder<Arg2> >::type
	ewise(  const Fun& fun,
			const typename Fun::first_arg_type& arg1v,
			const IMatrixXpr<Arg2, typename Fun::second_arg_type>& arg2 )
	{
		const_matrix<typename Fun::first_arg_type,
			ct_rows<Arg2>::value,
			ct_cols<Arg2>::value>
		arg1(arg2.nrows(), arg2.ncolumns(), arg1v);

		return make_expr(ewise(fun),
				copy_arg(arg1),
				ref_arg(arg2.derived()));
	}

}

#endif 




