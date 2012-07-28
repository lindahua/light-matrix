/**
 * @file matrix_reduce_internal.h
 *
 * Internal implementation of matrix reduction
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_REDUCE_INTERNAL_H_
#define LIGHTMAT_MATRIX_REDUCE_INTERNAL_H_

#include <light_mat/matrix/matrix_properties.h>
#include <light_mat/matrix/matrix_veval.h>

namespace lmat { namespace detail {

	template<int CTLen>
	struct single_vec_reduce
	{
		template<class Fun, class Vec>
		LMAT_ENSURE_INLINE
		static typename Fun::result_type eval(const Fun& fun, const index_t, const Vec& vec)
		{
			typename Fun::result_type r = fun(vec.get_value(0));
			for (index_t i = 1; i < CTLen; ++i)
			{
				r = fun(r, vec.get_value(i));
			}
			return r;
		}
	};

	template<>
	struct single_vec_reduce<0>
	{
		template<class Fun, class Vec>
		LMAT_ENSURE_INLINE
		static typename Fun::result_type eval(const Fun& fun, const index_t len, const Vec& vec)
		{
			if (len > 0)
			{
				typename Fun::result_type r = fun(vec.get_value(0));
				for (index_t i = 1; i < len; ++i)
				{
					r = fun(r, vec.get_value(i));
				}
				return r;
			}
			else
			{
				return fun();
			}
		}
	};

	template<>
	struct single_vec_reduce<1>
	{
		template<class Fun, class Vec>
		LMAT_ENSURE_INLINE
		static typename Fun::result_type eval(const Fun& fun, const index_t, const Vec& vec)
		{
			return fun(vec.get_value(0));
		}
	};


	struct full_reduce_linear_internal
	{
		template<class Fun, class Expr>
		LMAT_ENSURE_INLINE
		static typename Fun::result_type evaluate(const Fun& fun, const Expr& expr)
		{
			typename linear_eval<Expr>::evaluator_type evaluator(expr);
			return single_vec_reduce<ct_size<Expr>::value>::eval(
					fun, expr.nelems(), evaluator);
		}
	};


	struct full_reduce_percol_internal
	{
		template<class Fun, class Expr>
		inline
		static typename Fun::result_type evaluate(const Fun& fun, const Expr& expr)
		{
			typename percol_eval<Expr>::evaluator_type evaluator(expr);
			const index_t m = expr.nrows();
			const index_t n = expr.ncolumns();

			typedef typename Fun::result_type RT;

			RT r = single_vec_reduce<ct_rows<Expr>::value>::eval(
					fun, m, evaluator);

			for (index_t j = 1; j < n; ++j)
			{
				evaluator.next_column();

				RT rj = single_vec_reduce<ct_rows<Expr>::value>::eval(
						fun, m, evaluator);

				r = fun(r, rj);
			}

			return r;
		}
	};


	template<class Expr>
	struct reduce_internal_map
	{
		static const int linear_cost = linear_eval<Expr>::cost;

		static const bool is_shortv = ct_rows<Expr>::value <= SHORTVEC_LENGTH_THRESHOLD;
		static const int percol_cost = is_shortv ?
				percol_eval<Expr>::shortv_cost : percol_eval<Expr>::normal_cost;

		static const bool use_percol = percol_cost < linear_cost;

		typedef typename
				if_c<use_percol,
					full_reduce_percol_internal,
					full_reduce_linear_internal>::type type;
	};


} }

#endif
