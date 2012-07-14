/**
 * @file partial_reduce_internal.h
 *
 * Internal implementation of partial reduction evaluation
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_PARTIAL_REDUCE_INTERNAL_H_
#define LIGHTMAT_PARTIAL_REDUCE_INTERNAL_H_

#include "matrix_reduce_internal.h"

namespace lmat { namespace detail {

	struct colwise_reduce_internal
	{
		template<class Fun, class Arg, class Dst>
		inline
		static void eval(const Fun& fun, const Arg& arg, Dst& dst)
		{
			const index_t m = arg.nrows();
			const index_t n = arg.ncolumns();

			typename percol_eval<Arg>::evaluator_type evaluator(arg);

			for (index_t j = 0; j < n; ++j, evaluator.next_column())
			{
				dst.elem(0, j) = single_vec_reduce<ct_rows<Arg>::value>::eval(
						fun, m, evaluator);
			}
		}
	};


	template<int CTLen>
	struct accum_to_row
	{
		template<class Fun, class Vec>
		LMAT_ENSURE_INLINE
		static void eval(const Fun& fun, const index_t, const Vec& vec, typename Fun::result_type *pd)
		{
			for (index_t i = 0; i < CTLen; ++i)
			{
				pd[i] = fun(pd[i], vec.get_value(i));
			}
		}

		template<class Fun, class Vec>
		LMAT_ENSURE_INLINE
		static void init_eval(const Fun& fun, const index_t, const Vec& vec, typename Fun::result_type *pd)
		{
			for (index_t i = 0; i < CTLen; ++i)
			{
				pd[i] = fun(vec.get_value(i));
			}
		}

		template<class Fun>
		LMAT_ENSURE_INLINE
		static void empty_eval(const Fun& fun, const index_t, typename Fun::result_type *pd)
		{
			for (index_t i = 0; i < CTLen; ++i)
			{
				pd[i] = fun();
			}
		}
	};

	template<>
	struct accum_to_row<0>
	{
		template<class Fun, class Vec>
		LMAT_ENSURE_INLINE
		static void eval(const Fun& fun, const index_t m, const Vec& vec, typename Fun::result_type *pd)
		{
			for (index_t i = 0; i < m; ++i)
			{
				pd[i] = fun(pd[i], vec.get_value(i));
			}
		}

		template<class Fun, class Vec>
		LMAT_ENSURE_INLINE
		static void init_eval(const Fun& fun, const index_t m, const Vec& vec, typename Fun::result_type *pd)
		{
			for (index_t i = 0; i < m; ++i)
			{
				pd[i] = fun(vec.get_value(i));
			}
		}

		template<class Fun>
		LMAT_ENSURE_INLINE
		static void empty_eval(const Fun& fun, const index_t m, typename Fun::result_type *pd)
		{
			for (index_t i = 0; i < m; ++i)
			{
				pd[i] = fun();
			}
		}
	};

	template<>
	struct accum_to_row<1>
	{
		template<class Fun, class Vec>
		LMAT_ENSURE_INLINE
		static void eval(const Fun& fun, const index_t, const Vec& vec, typename Fun::result_type *pd)
		{
			pd[0] = fun(pd[0], vec.get_value(0));
		}

		template<class Fun, class Vec>
		LMAT_ENSURE_INLINE
		static void init_eval(const Fun& fun, const index_t, const Vec& vec, typename Fun::result_type *pd)
		{
			pd[0] = fun(vec.get_value(0));
		}

		template<class Fun>
		LMAT_ENSURE_INLINE
		static void empty_eval(const Fun& fun, const index_t, typename Fun::result_type *pd)
		{
			pd[0] = fun();
		}
	};



	struct rowwise_reduce_internal
	{
		template<class Fun, class Arg, class Dst>
		inline
		static void eval(const Fun& fun, const Arg& arg, Dst& dst)
		{
			typedef accum_to_row<ct_rows<Arg>::value> accum_t;
			typename matrix_traits<Dst>::value_type* pd = dst.ptr_data();

			const index_t m = arg.nrows();
			const index_t n = arg.ncolumns();

			if (n > 0)
			{
				typename percol_eval<Arg>::evaluator_type evaluator(arg);

				accum_t::init_eval(fun, m, evaluator, pd);

				if (n > 1)
				{
					evaluator.next_column();
					for (index_t j = 1; j < n; ++j)
					{
						accum_t::eval(fun, m, evaluator, pd);
						evaluator.next_column();
					}
				}
			}
			else
			{
				accum_t::empty_eval(fun, m, pd);
			}

		}
	};


} }

#endif
