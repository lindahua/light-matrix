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

			typename percol_vector_evaluator<Arg>::type evaluator(arg);

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
	};



	struct rowwise_reduce_internal
	{
		template<class Fun, class Arg, class Dst>
		inline
		static void eval(const Fun& fun, const Arg& arg, Dst& dst)
		{
			const index_t m = arg.nrows();
			const index_t n = arg.ncolumns();

			typename percol_vector_evaluator<Arg>::type evaluator(arg);
			typename matrix_traits<Dst>::value_type* pd = dst.ptr_data();

			for (index_t j = 0; j < n; ++j, evaluator.next_column())
			{
				accum_to_row<ct_rows<Arg>::value>::eval(fun, m, evaluator, pd);
			}
		}
	};


} }

#endif
