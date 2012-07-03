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

#include <light_mat/matrix/matrix_concepts.h>
#include <light_mat/math/functor_base.h>

namespace lmat
{

	/********************************************
	 *
	 *  Expression traits
	 *
	 ********************************************/
	template<class Fun, typename Arg>
	struct matrix_traits<unary_ewise_expr<Fun, Arg> >
	{
		static const int num_dimensions = 2;
		static const int compile_time_num_rows = ct_rows<Arg>::value;
		static const int compile_time_num_cols = ct_cols<Arg>::value;

		static const bool is_readonly = true;

		typedef typename Fun::result_type value_type;
	};

	template<class Fun, typename Arg1, typename Arg2>
	struct matrix_traits<binary_ewise_expr<Fun, Arg1, Arg2> >
	{
		static const int num_dimensions = 2;
		static const int compile_time_num_rows = binary_ct_rows<Arg1, Arg2>::value;
		static const int compile_time_num_cols = binary_ct_cols<Arg1, Arg2>::value;

		static const bool is_readonly = true;

		typedef typename Fun::result_type value_type;
	};

	template<class Fun, typename Arg>
	struct matrix_traits<binary_fix1_ewise_expr<Fun, Arg> >
	{
		static const int num_dimensions = 2;
		static const int compile_time_num_rows = ct_rows<Arg>::value;
		static const int compile_time_num_cols = ct_cols<Arg>::value;

		static const bool is_readonly = true;

		typedef typename Fun::result_type value_type;
	};

	template<class Fun, typename Arg>
	struct matrix_traits<binary_fix2_ewise_expr<Fun, Arg> >
	{
		static const int num_dimensions = 2;
		static const int compile_time_num_rows = ct_rows<Arg>::value;
		static const int compile_time_num_cols = ct_cols<Arg>::value;

		static const bool is_readonly = true;

		typedef typename Fun::result_type value_type;
	};


	/********************************************
	 *
	 *  Expression classes
	 *
	 ********************************************/

	template<class Fun, typename Arg>
	class unary_ewise_expr
	: public IMatrixXpr<unary_ewise_expr<Fun, Arg>, typename Fun::result_type>
	{
#ifdef LMAT_USE_STATIC_ASSERT
		static_assert(is_unary_ewise_functor<Fun>::value, "Fun must be a unary_ewise_functor");
		static_assert(is_mat_xpr<Arg>::value, "Arg must be a matrix expression class.");
#endif

	public:
		typedef typename Fun::result_type value_type;

		LMAT_ENSURE_INLINE
		unary_ewise_expr(const Fun& fun, const Arg& a)
		: m_fun(fun), m_arg(a) { }

		LMAT_ENSURE_INLINE const Fun& fun() const
		{
			return m_fun;
		}

		LMAT_ENSURE_INLINE const Arg& arg() const
		{
			return m_arg;
		}

		LMAT_ENSURE_INLINE index_t nelems() const
		{
			return m_arg.nelems();
		}

		LMAT_ENSURE_INLINE size_t size() const
		{
			return m_arg.size();
		}

		LMAT_ENSURE_INLINE index_t nrows() const
		{
			return m_arg.nrows();
		}

		LMAT_ENSURE_INLINE index_t ncolumns() const
		{
			return m_arg.ncolumns();
		}

	private:
		Fun m_fun;
		const Arg& m_arg;
	};

	template<class Fun, typename Arg1, typename Arg2>
	class binary_ewise_expr
	: public IMatrixXpr<binary_ewise_expr<Fun, Arg1, Arg2>, typename Fun::result_type>
	{
#ifdef LMAT_USE_STATIC_ASSERT
		static_assert(is_binary_ewise_functor<Fun>::value, "Fun must be a binary_ewise_functor");
		static_assert(is_mat_xpr<Arg1>::value, "Arg1 must be a matrix expression class.");
		static_assert(is_mat_xpr<Arg2>::value, "Arg2 must be a matrix expression class.");
#endif

	public:
		typedef typename Fun::result_type value_type;

		LMAT_ENSURE_INLINE
		binary_ewise_expr(const Fun& fun, const Arg1& arg1, const Arg2& arg2)
		: m_fun(fun), m_arg1(arg1), m_arg2(arg2)
		{
			check_same_size(arg1, arg2, "arg1 and arg2 must be of the same size.");
		}

		LMAT_ENSURE_INLINE const Fun& fun() const
		{
			return m_fun;
		}

		LMAT_ENSURE_INLINE const Arg1& first_arg() const
		{
			return m_arg1;
		}

		LMAT_ENSURE_INLINE const Arg2& second_arg() const
		{
			return m_arg2;
		}

		LMAT_ENSURE_INLINE index_t nelems() const
		{
			return m_arg1.nelems();
		}

		LMAT_ENSURE_INLINE size_t size() const
		{
			return m_arg1.size();
		}

		LMAT_ENSURE_INLINE index_t nrows() const
		{
			return m_arg1.nrows();
		}

		LMAT_ENSURE_INLINE index_t ncolumns() const
		{
			return m_arg1.ncolumns();
		}

	private:
		Fun m_fun;
		const Arg1& m_arg1;
		const Arg2& m_arg2;
	};


	template<class Fun, typename Arg2>
	class binary_fix1_ewise_expr
	: public IMatrixXpr<binary_fix1_ewise_expr<Fun, Arg2>, typename Fun::result_type>
	{
#ifdef LMAT_USE_STATIC_ASSERT
		static_assert(is_binary_ewise_functor<Fun>::value, "Fun must be a binary_ewise_functor");
		static_assert(is_mat_xpr<Arg2>::value, "Arg2 must be a matrix expression class.");
#endif

	public:
		typedef typename Fun::result_type value_type;
		typedef typename Fun::first_arg_type arg1_vtype;

		LMAT_ENSURE_INLINE
		binary_fix1_ewise_expr(const Fun& fun, const arg1_vtype& arg1v, const Arg2& arg2)
		: m_fun(fun), m_arg1v(arg1v), m_arg2(arg2)
		{
		}

		LMAT_ENSURE_INLINE const Fun& fun() const
		{
			return m_fun;
		}

		LMAT_ENSURE_INLINE const arg1_vtype& first_arg_value() const
		{
			return m_arg1v;
		}

		LMAT_ENSURE_INLINE const Arg2& second_arg() const
		{
			return m_arg2;
		}

		LMAT_ENSURE_INLINE index_t nelems() const
		{
			return m_arg2.nelems();
		}

		LMAT_ENSURE_INLINE size_t size() const
		{
			return m_arg2.size();
		}

		LMAT_ENSURE_INLINE index_t nrows() const
		{
			return m_arg2.nrows();
		}

		LMAT_ENSURE_INLINE index_t ncolumns() const
		{
			return m_arg2.ncolumns();
		}

	private:
		Fun m_fun;
		const arg1_vtype m_arg1v;
		const Arg2& m_arg2;
	};


	template<class Fun, typename Arg1>
	class binary_fix2_ewise_expr
	: public IMatrixXpr<binary_fix2_ewise_expr<Fun, Arg1>, typename Fun::result_type>
	{
#ifdef LMAT_USE_STATIC_ASSERT
		static_assert(is_binary_ewise_functor<Fun>::value, "Fun must be a binary_ewise_functor");
		static_assert(is_mat_xpr<Arg1>::value, "Arg1 must be a matrix expression class.");
#endif

	public:
		typedef typename Fun::result_type value_type;
		typedef typename Fun::second_arg_type arg2_vtype;

		LMAT_ENSURE_INLINE
		binary_fix2_ewise_expr(const Fun& fun, const Arg1& arg1, const arg2_vtype& arg2v)
		: m_fun(fun), m_arg1(arg1), m_arg2v(arg2v)
		{
		}

		LMAT_ENSURE_INLINE const Fun& fun() const
		{
			return m_fun;
		}

		LMAT_ENSURE_INLINE const Arg1& first_arg() const
		{
			return m_arg1;
		}

		LMAT_ENSURE_INLINE const arg2_vtype& second_arg_value() const
		{
			return m_arg2v;
		}

		LMAT_ENSURE_INLINE index_t nelems() const
		{
			return m_arg1.nelems();
		}

		LMAT_ENSURE_INLINE size_t size() const
		{
			return m_arg1.size();
		}

		LMAT_ENSURE_INLINE index_t nrows() const
		{
			return m_arg1.nrows();
		}

		LMAT_ENSURE_INLINE index_t ncolumns() const
		{
			return m_arg1.ncolumns();
		}

	private:
		Fun m_fun;
		const Arg1& m_arg1;
		const arg2_vtype m_arg2v;
	};


	/********************************************
	 *
	 *  Expression type mapping
	 *
	 ********************************************/

	template<class Fun, class Arg>
	struct unary_ewise_expr_map
	{
		typedef unary_ewise_expr<Fun, Arg> type;
	};

	template<class Fun, class Arg1, class Arg2>
	struct binary_ewise_expr_map
	{
		typedef binary_ewise_expr<Fun, Arg1, Arg2> type;
	};

	template<class Fun, class Arg1>
	struct binary_fix2_ewise_expr_map
	{
		typedef binary_fix2_ewise_expr<Fun, Arg1> type;
	};

	template<class Fun, class Arg2>
	struct binary_fix1_ewise_expr_map
	{
		typedef binary_fix1_ewise_expr<Fun, Arg2> type;
	};


	/********************************************
	 *
	 *  Expression construction
	 *
	 ********************************************/

	template<class Fun, class Arg>
	LMAT_ENSURE_INLINE
	inline typename unary_ewise_expr_map<Fun, Arg>::type
	ewise(  const Fun& fun,
			const IMatrixXpr<Arg, typename Fun::arg_type>& arg )
	{
		return unary_ewise_expr<Fun, Arg>(fun, arg.derived());
	}

	template<class Fun, class Arg1, class Arg2>
	LMAT_ENSURE_INLINE
	inline typename binary_ewise_expr_map<Fun, Arg1, Arg2>::type
	ewise(  const Fun& fun,
			const IMatrixXpr<Arg1, typename Fun::first_arg_type>& arg1,
			const IMatrixXpr<Arg2, typename Fun::second_arg_type>& arg2 )
	{
		return binary_ewise_expr<Fun, Arg1, Arg2>(fun,
				arg1.derived(), arg2.derived());
	}

	template<class Fun, class Arg1>
	LMAT_ENSURE_INLINE
	inline typename binary_fix2_ewise_expr_map<Fun, Arg1>::type
	ewise(  const Fun& fun,
			const IMatrixXpr<Arg1, typename Fun::first_arg_type>& arg1,
			const typename Fun::second_arg_type& arg2v )
	{
		return binary_fix2_ewise_expr<Fun, Arg1>(fun,
				arg1.derived(), arg2v);
	}


	template<class Fun, class Arg2>
	LMAT_ENSURE_INLINE
	inline typename binary_fix1_ewise_expr_map<Fun, Arg2>::type
	ewise(  const Fun& fun,
			const typename Fun::first_arg_type& arg1v,
			const IMatrixXpr<Arg2, typename Fun::second_arg_type>& arg2 )
	{
		return binary_fix1_ewise_expr<Fun, Arg2>(fun,
				arg1v, arg2.derived());
	}

	/********************************************
	 *
	 *  Conversion expression
	 *
	 ********************************************/

	template<class SMat, typename T>
	struct convert_expr_map
	{
		typedef typename matrix_traits<SMat>::value_type S;
		typedef type_converter<S, T> converter_t;
		typedef typename unary_ewise_expr_map<converter_t, SMat>::type type;
	};

	template<class SMat, typename S, typename T>
	LMAT_ENSURE_INLINE
	inline typename convert_expr_map<SMat, T>::type
	cast(const IMatrixXpr<SMat, S>& smat, type<T> )
	{
		return ewise(type_converter<S, T>(), smat);
	}

}

#endif 




