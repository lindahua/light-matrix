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

	template<typename KerCate> struct colwise_reduce_internal;
	template<typename KerCate> struct rowwise_reduce_internal;


	template<>
	struct colwise_reduce_internal<scalar_kernel_t>
	{
		template<class Fun, class Arg, class Dst>
		inline
		static void eval(const Fun& fun, const Arg& arg, Dst& dst)
		{
			typedef matrix_visit_setting<
					percol_vis,
					scalar_kernel_t,
					ct_rows<Arg>::value,
					binary_ct_cols<Arg, Dst>::value> setting_t;

			const index_t m = arg.nrows();
			const index_t n = arg.ncolumns();

			typedef typename matrix_vismap<Arg, setting_t>::type visitor_t;
			visitor_t visitor(arg);

			for (index_t j = 0; j < n; ++j)
			{
				typename matrix_visitor_state<visitor_t>::type s = visitor.col_state(j);

				dst.elem(0, j) = single_vec_reduce<
							ct_rows<Arg>::value,
							scalar_kernel_t>::eval(fun, m, visitor, s);
			}
		}
	};



	template<int CTLen, typename KerCate> struct accum_col;

	template<int CTLen>
	struct accum_col<CTLen, scalar_kernel_t>
	{
		template<class Fun, class Visitor, class State>
		LMAT_ENSURE_INLINE
		static void eval(const Fun& fun, const index_t, const Visitor& vis, const State& s,
				typename Fun::result_type *pd)
		{
			for (index_t i = 0; i < CTLen; ++i)
			{
				pd[i] = fun(pd[i], vis.get_scalar(i, s));
			}
		}

		template<class Fun, class Visitor, class State>
		LMAT_ENSURE_INLINE
		static void init_eval(const Fun& fun, const index_t, const Visitor& vis, const State& s,
				typename Fun::result_type *pd)
		{
			for (index_t i = 0; i < CTLen; ++i)
			{
				pd[i] = fun(vis.get_scalar(i, s));
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
	struct accum_col<0, scalar_kernel_t>
	{
		template<class Fun, class Visitor, typename State>
		LMAT_ENSURE_INLINE
		static void eval(const Fun& fun, const index_t m, const Visitor& vis, const State& s,
				typename Fun::result_type *pd)
		{
			for (index_t i = 0; i < m; ++i)
			{
				pd[i] = fun(pd[i], vis.get_scalar(i, s));
			}
		}

		template<class Fun, class Visitor, typename State>
		LMAT_ENSURE_INLINE
		static void init_eval(const Fun& fun, const index_t m, const Visitor& vis, const State& s,
				typename Fun::result_type *pd)
		{
			for (index_t i = 0; i < m; ++i)
			{
				pd[i] = fun(vis.get_scalar(i, s));
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
	struct accum_col<1, scalar_kernel_t>
	{
		template<class Fun, class Visitor, typename State>
		LMAT_ENSURE_INLINE
		static void eval(const Fun& fun, const index_t, const Visitor& vis, const State& s,
				typename Fun::result_type *pd)
		{
			pd[0] = fun(pd[0], vis.get_scalar(0, s));
		}

		template<class Fun, class Visitor, typename State>
		LMAT_ENSURE_INLINE
		static void init_eval(const Fun& fun, const index_t, const Visitor& vis, const State& s,
				typename Fun::result_type *pd)
		{
			pd[0] = fun(vis.get_scalar(0, s));
		}

		template<class Fun>
		LMAT_ENSURE_INLINE
		static void empty_eval(const Fun& fun, const index_t, typename Fun::result_type *pd)
		{
			pd[0] = fun();
		}
	};


	template<>
	struct rowwise_reduce_internal<scalar_kernel_t>
	{
		template<class Fun, class Arg, class Dst>
		inline
		static void eval(const Fun& fun, const Arg& arg, Dst& dst)
		{
			typedef accum_col<ct_rows<Arg>::value, scalar_kernel_t> accum_t;

			typedef matrix_visit_setting<
					percol_vis,
					scalar_kernel_t,
					binary_ct_rows<Arg, Dst>::value,
					ct_cols<Arg>::value> setting_t;

			typedef typename matrix_vismap<Arg, setting_t>::type visitor_t;
			typedef typename matrix_visitor_state<visitor_t>::type state_t;

			typename matrix_traits<Dst>::value_type* pd = dst.ptr_data();

			const index_t m = arg.nrows();
			const index_t n = arg.ncolumns();

			if (n > 0)
			{
				visitor_t visitor(arg);

				state_t s0 = visitor.col_state(0);
				accum_t::init_eval(fun, m, visitor, s0, pd);

				if (n > 1)
				{
					for (index_t j = 1; j < n; ++j)
					{
						state_t s = visitor.col_state(j);
						accum_t::eval(fun, m, visitor, s, pd);
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
