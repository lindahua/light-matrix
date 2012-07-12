/**
 * @file partial_reduce.h
 *
 * Partial reduction expression and evaluation
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_PARTIAL_REDUCE_H_
#define LIGHTMAT_PARTIAL_REDUCE_H_

#include <light_mat/matrix/matrix_arith.h>
#include <light_mat/math/reduction_functors.h>
#include "bits/partial_reduce_internal.h"

namespace lmat
{
	/********************************************
	 *
	 *  Expression classes
	 *
	 ********************************************/

	template<class Fun, class Arg>
	struct matrix_traits<colwise_reduce_expr<Fun, Arg> >
	{
		static const int num_dimensions = 2;
		static const int compile_time_num_rows = 1;
		static const int compile_time_num_cols = ct_cols<Arg>::value;

		static const bool is_readonly = true;

		typedef typename Fun::result_type value_type;
	};


	template<class Fun, class Arg>
	struct matrix_traits<rowwise_reduce_expr<Fun, Arg> >
	{
		static const int num_dimensions = 2;
		static const int compile_time_num_rows = ct_rows<Arg>::value;
		static const int compile_time_num_cols = 1;

		static const bool is_readonly = true;

		typedef typename Fun::result_type value_type;
	};


	template<class Fun, class Arg>
	class colwise_reduce_expr
	: public IMatrixXpr<colwise_reduce_expr<Fun, Arg>, typename Fun::result_type>
	{
#ifdef LMAT_USE_STATIC_ASSERT
		static_assert(is_reduction_functor<Fun>::value, "Fun must be a reduction_functor");
		static_assert(is_mat_xpr<Arg>::value, "Arg must be a matrix expression class.");
#endif

	public:
		typedef typename Fun::result_type value_type;

		LMAT_ENSURE_INLINE
		colwise_reduce_expr(const Fun& fun, const Arg& a)
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
			return m_arg.ncolumns();
		}

		LMAT_ENSURE_INLINE size_t size() const
		{
			return static_cast<size_t>(nelems());
		}

		LMAT_ENSURE_INLINE index_t nrows() const
		{
			return 1;
		}

		LMAT_ENSURE_INLINE index_t ncolumns() const
		{
			return m_arg.ncolumns();
		}

	private:
		Fun m_fun;
		const Arg& m_arg;
	};


	template<class Fun, class Arg>
	class rowwise_reduce_expr
	: public IMatrixXpr<rowwise_reduce_expr<Fun, Arg>, typename Fun::result_type>
	{
#ifdef LMAT_USE_STATIC_ASSERT
		static_assert(is_reduction_functor<Fun>::value, "Fun must be a reduction_functor");
		static_assert(is_mat_xpr<Arg>::value, "Arg must be a matrix expression class.");
#endif

	public:
		typedef typename Fun::result_type value_type;

		LMAT_ENSURE_INLINE
		rowwise_reduce_expr(const Fun& fun, const Arg& a)
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
			return m_arg.nrows();
		}

		LMAT_ENSURE_INLINE size_t size() const
		{
			return static_cast<size_t>(nelems());
		}

		LMAT_ENSURE_INLINE index_t nrows() const
		{
			return m_arg.nrows();
		}

		LMAT_ENSURE_INLINE index_t ncolumns() const
		{
			return 1;
		}

	private:
		Fun m_fun;
		const Arg& m_arg;
	};


	/********************************************
	 *
	 *  Generic expressions
	 *
	 ********************************************/

	template<class Fun, class Arg>
	struct colwise_reduce_expr_map
	{
		typedef colwise_reduce_expr<Fun, Arg> type;
	};

	template<class Fun, class Arg>
	struct rowwise_reduce_expr_map
	{
		typedef rowwise_reduce_expr<Fun, Arg> type;
	};

	template<class Fun, class Arg>
	LMAT_ENSURE_INLINE
	inline typename colwise_reduce_expr_map<Fun, Arg>::type
	reduce(  const Fun& fun,
			const IMatrixXpr<Arg, typename Fun::arg_type>& arg,
			colwise )
	{
		return colwise_reduce_expr<Fun, Arg>(fun, arg.derived());
	}

	template<class Fun, class Arg>
	LMAT_ENSURE_INLINE
	inline typename rowwise_reduce_expr_map<Fun, Arg>::type
	reduce(  const Fun& fun,
			const IMatrixXpr<Arg, typename Fun::arg_type>& arg,
			rowwise )
	{
		return rowwise_reduce_expr<Fun, Arg>(fun, arg.derived());
	}


	template<class Fun, class Arg, class DMat>
	inline void evaluate_to(const colwise_reduce_expr<Fun, Arg>& expr,
			IDenseMatrix<DMat, typename Fun::result_type>& dst)
	{
		detail::colwise_reduce_internal::eval(expr.fun(), expr.arg(), dst.derived());
	}


	template<class Fun, class Arg, class DMat>
	inline void evaluate_to(const rowwise_reduce_expr<Fun, Arg>& expr,
			IDenseMatrix<DMat, typename Fun::result_type>& dst)
	{
		detail::rowwise_reduce_internal::eval(expr.fun(), expr.arg(), dst.derived());
	}


	/********************************************
	 *
	 *  Specific expressions
	 *
	 ********************************************/

	// sum

	template<typename T, class Arg>
	struct colwise_sum_t
	{
		typedef typename colwise_reduce_expr_map<
				sum_fun<T>, Arg>::type type;
	};

	template<typename T, class Arg>
	struct rowwise_sum_t
	{
		typedef typename rowwise_reduce_expr_map<
				sum_fun<T>, Arg>::type type;
	};

	template<typename T, class Arg>
	LMAT_ENSURE_INLINE
	inline typename colwise_sum_t<T, Arg>::type
	sum(const IMatrixXpr<Arg, T>& arg, colwise)
	{
		return reduce(sum_fun<T>(), arg, colwise());
	}

	template<typename T, class Arg>
	LMAT_ENSURE_INLINE
	inline typename rowwise_sum_t<T, Arg>::type
	sum(const IMatrixXpr<Arg, T>& arg, rowwise)
	{
		return reduce(sum_fun<T>(), arg, rowwise());
	}

	// mean

	template<typename T, class Arg>
	struct colwise_mean_t
	{
		typedef typename
				binary_fix2_ewise_expr_map<
					mul_op<T>,
					typename colwise_sum_t<T, Arg>::type
				>::type type;
	};

	template<typename T, class Arg>
	struct rowwise_mean_t
	{
		typedef typename
				binary_fix2_ewise_expr_map<
					mul_op<T>,
					typename rowwise_sum_t<T, Arg>::type
				>::type type;
	};

}

#endif
