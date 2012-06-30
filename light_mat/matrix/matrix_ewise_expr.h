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
	 *  Expression classes
	 *
	 ********************************************/

	template<class Fun, typename Arg>
	class simple_unary_ewise_expr
	: public IMatrixXpr<simple_unary_ewise_expr<Fun, Arg>, typename Fun::result_type>
	{
	public:
		typedef typename Fun::result_type value_type;

		LMAT_ENSURE_INLINE
		explicit simple_unary_ewise_expr(const Arg& a) : m_arg(a) { }

		LMAT_ENSURE_INLINE const Arg& arg() const
		{
			return m_arg;
		}

	private:
		const Arg& m_arg;
	};

	template<class Fun, typename Arg>
	class unary_ewise_expr
	: public IMatrixXpr<unary_ewise_expr<Fun, Arg>, typename Fun::result_type>
	{
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

	private:
		const Fun& m_fun;
		const Arg& m_arg;
	};

	template<class Fun, typename Arg1, typename Arg2>
	class simple_binary_ewise_expr
	: public IMatrixXpr<simple_binary_ewise_expr<Fun, Arg1, Arg2>, typename Fun::result_type>
	{
	public:
		typedef typename Fun::result_type value_type;

		LMAT_ENSURE_INLINE
		simple_binary_ewise_expr(const Arg1& arg1, const Arg2& arg2)
		: m_arg1(arg1), m_arg2(arg2) { }

		LMAT_ENSURE_INLINE const Arg1& first_arg() const
		{
			return m_arg1;
		}

		LMAT_ENSURE_INLINE const Arg2& second_arg() const
		{
			return m_arg2;
		}

	private:
		const Arg1& m_arg1;
		const Arg2& m_arg2;
	};

	template<class Fun, typename Arg1, typename Arg2>
	class binary_ewise_expr
	: public IMatrixXpr<binary_ewise_expr<Fun, Arg1, Arg2>, typename Fun::result_type>
	{
	public:
		typedef typename Fun::result_type value_type;

		LMAT_ENSURE_INLINE
		binary_ewise_expr(const Fun& fun, const Arg1& arg1, const Arg2& arg2)
		: m_fun(fun), m_arg1(arg1), m_arg2(arg2) { }

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

	private:
		const Fun& m_fun;
		const Arg1& m_arg1;
		const Arg2& m_arg2;
	};


	/********************************************
	 *
	 *  Expression construction
	 *
	 ********************************************/

	template<class Fun, class Arg>
	LMAT_ENSURE_INLINE
	inline simple_unary_ewise_expr<Fun, Arg> simple_ewise_expr(
			Fun,
			const IMatrixXpr<Arg, typename Fun::arg_type>& arg)
	{
		return simple_unary_ewise_expr<Fun, Arg>(arg.derived());
	}

	template<class Fun, class Arg1, class Arg2>
	LMAT_ENSURE_INLINE
	inline simple_binary_ewise_expr<Fun, Arg1, Arg2> simple_ewise_expr(
			Fun,
			const IMatrixXpr<Arg1, typename Fun::first_arg_type>& arg1,
			const IMatrixXpr<Arg2, typename Fun::second_arg_type>& arg2)
	{
		return simple_binary_ewise_expr<Fun, Arg1, Arg2>(
				arg1.derived(), arg2.derived());
	}


	template<class Fun, class Arg>
	LMAT_ENSURE_INLINE
	inline unary_ewise_expr<Fun, Arg> ewise_expr(
			const Fun& fun,
			const IMatrixXpr<Arg, typename Fun::arg_type>& arg)
	{
		return unary_ewise_expr<Fun, Arg>(fun, arg.derived());
	}

	template<class Fun, class Arg1, class Arg2>
	LMAT_ENSURE_INLINE
	inline binary_ewise_expr<Fun, Arg1, Arg2> simple_ewise_expr(
			const Fun& fun,
			const IMatrixXpr<Arg1, typename Fun::first_arg_type>& arg1,
			const IMatrixXpr<Arg2, typename Fun::second_arg_type>& arg2)
	{
		return binary_ewise_expr<Fun, Arg1, Arg2>(fun,
				arg1.derived(), arg2.derived());
	}

}

#endif 




