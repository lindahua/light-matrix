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

#include <light_mat/matexpr/matrix_access_base.h>
#include <light_mat/matrix/matrix_fill.h>

#include "macc_reduce_core.h"

namespace lmat { namespace internal {

	template<typename KerCate> struct colwise_reduce_internal;
	template<typename InterT, typename KerCate> struct rowwise_reduce_internal;


	template<>
	struct colwise_reduce_internal<scalar_kernel_t>
	{
		template<class Fun, class Arg, class Dst>
		static void eval(const Fun& fun, const Arg& a, Dst& dst)
		{
			typedef typename matrix_traits<Arg>::value_type T;
			typedef typename matrix_traits<Dst>::value_type RT;
			typedef typename macc_accessor_map<Arg, percol_macc, scalar_kernel_t>::type accessor_t;
			typedef typename percol_macc_state_map<accessor_t>::type col_state_t;

			typedef percol_to_linear_accessor<accessor_t, scalar_kernel_t, T> acc1_wrapper;

			const index_t m = a.nrows();
			const index_t n = a.ncolumns();

			accessor_t acc(a);

			if (m > 0)
			{
				if (m == 1)
				{
					for (index_t j = 0; j < n; ++j)
					{
						dst(0, j) = fun.get(
								fun.transform(acc.get_scalar(0, acc.col_state(j))), 1);
					}
				}
				else
				{
					typedef vec_reduce<RT, ct_rows<Arg>::value, scalar_kernel_t> vimpl_t;
					for (index_t j = 0; j < n; ++j)
					{
						dst(0, j) = fun.get(vimpl_t::eval_s(fun, m, acc, acc.col_state(j)), m);
					}
				}
			}
			else
			{
				if (n > 0)
					fill(dst, fun.empty_value());
			}
		}


		template<class Fun, class Arg1, class Arg2, class Dst>
		static void eval(const Fun& fun, const Arg1& a1, const Arg2& a2, Dst& dst)
		{
			typedef typename matrix_traits<Dst>::value_type RT;
			typedef typename macc_accessor_map<Arg1, percol_macc, scalar_kernel_t>::type accessor1_t;
			typedef typename macc_accessor_map<Arg2, percol_macc, scalar_kernel_t>::type accessor2_t;
			typedef typename percol_macc_state_map<accessor1_t>::type col_state1_t;
			typedef typename percol_macc_state_map<accessor2_t>::type col_state2_t;

			const index_t m = a1.nrows();
			const index_t n = a2.ncolumns();

			accessor1_t acc1(a1);
			accessor2_t acc2(a2);

			if (m > 0)
			{
				if (m == 1)
				{
					for (index_t j = 0; j < n; ++j)
					{
						dst(0, j) = fun.get(
								fun.transform(
										acc1.get_scalar(0, acc1.col_state(j)),
										acc2.get_scalar(0, acc2.col_state(j))), 1);
					}
				}
				else
				{
					typedef vec_reduce<RT, common_ctrows<Arg1, Arg2>::value, scalar_kernel_t> vimpl_t;
					for (index_t j = 0; j < n; ++j)
					{
						dst(0, j) = fun.get(vimpl_t::eval_s(fun, m,
								acc1, acc1.col_state(j), acc2, acc2.col_state(j)), m);
					}
				}
			}
			else
			{
				if (n > 0)
					fill(dst, fun.empty_value());
			}
		}
	};


	template<int M, typename ColStateT, typename InterT, typename T>
	struct rowwise_reduce_mcol_impl
	{
		template<class Fun, class Accessor, class Dst>
		static void eval(const Fun& fun, const Accessor& acc,
				index_t m, index_t n,
				const T *pd, index_t rs, index_t cs)
		{
			dense_col<InterT, M> mvs(m);

			ColStateT s0 = acc.col_state(0);
			for (index_t i = 0; i < m; ++i)
			{
				mvs[i] = fun.transform(acc.get_scalar(i, s0));
			}

			for (index_t j = 1; j < n; ++j, pd += cs)
			{
				ColStateT s = acc.col_state(j);
				for (index_t i = 0; i < m; ++i)
				{
					mvs[i] = fun.combine(mvs[i], fun.transform(acc.get_scalar(i, s)));
				}
			}

			if (rs == 1)
			{
				for (index_t i = 0; i < m; ++i)
				{
					pd[i] = fun.get(mvs[i], n);
				}
			}
			else
			{
				for (index_t i = 0; i < m; ++i)
				{
					pd[i * rs] = fun.get(mvs[i], n);
				}
			}
		}
	};


	template<int M, typename ColStateT, typename T>
	struct rowwise_reduce_mcol_impl<M, ColStateT, T, T>
	{
		template<class Fun, class Accessor, class Dst>
		static void eval(const Fun& fun, const Accessor& acc,
				index_t m, index_t n,
				const T *pd, index_t rs, index_t cs)
		{

			if (rs == 1)
			{
				ColStateT s0 = acc.col_state(0);

				for (index_t i = 0; i < m; ++i)
				{
					pd[i] = fun.transform(acc.get_scalar(i, s0));
				}

				for (index_t j = 1; j < n; ++j)
				{
					ColStateT s = acc.col_state(j);
					for (index_t i = 0; i < m; ++i)
					{
						pd[i] = fun.combine(pd[i], fun.transform(acc.get_scalar(i, s)));
					}
				}

				for (index_t i = 0; i < m; ++i)
				{
					pd[i] = fun.get(pd[i], n);
				}
			}
			else
			{
				ColStateT s0 = acc.col_state(0);

				for (index_t i = 0; i < m; ++i)
				{
					pd[i * rs] = fun.transform(acc.get_scalar(i, s0));
				}

				for (index_t j = 1; j < n; ++j)
				{
					ColStateT s = acc.col_state(j);
					for (index_t i = 0; i < m; ++i)
					{
						pd[i * rs] = fun.combine(pd[i * rs], fun.transform(acc.get_scalar(i, s)));
					}
				}

				for (index_t i = 0; i < m; ++i)
				{
					pd[i * rs] = fun.get(pd[i * rs], n);
				}
			}

		}
	};




	template<typename InterT>
	struct rowwise_reduce_internal<InterT, scalar_kernel_t>
	{
		template<class Fun, class Arg, class Dst>
		static void eval(const Fun& fun, const Arg& a, Dst& dst)
		{
			typedef typename matrix_traits<Dst>::value_type RT;
			typedef typename macc_accessor_map<Arg, percol_macc, scalar_kernel_t>::type accessor_t;
			typedef typename percol_macc_state_map<accessor_t>::type col_state_t;
			typedef rowwise_reduce_mcol_impl<common_ctrows<Arg, Dst>::value, col_state_t, InterT, RT> mcol_impl_t;

			const index_t m = a.nrows();
			const index_t n = a.ncolumns();

			accessor_t acc(a);
			RT *pd = dst.ptr_data();


			if (n > 0)
			{
				if (n == 1)
				{
					col_state_t s0 = acc.col_state(0);

					if (m == 1)
					{
						*pd = fun.get(fun.transform(acc.get_scalar(0, s0)));
					}
					else
					{
						const index_t rs = dst.row_stride();
						if (rs == 1)
						{
							for (index_t i = 0; i < m; ++i)
							{
								pd[i] = fun.get(fun.transform(acc.get_scalar(i, s0)));
							}
						}
						else
						{
							for (index_t i = 0; i < m; ++i)
							{
								pd[i * rs] = fun.get(fun.transform(acc.get_scalar(i, s0)));
							}
						}
					}
				}
				else
				{
					const index_t cs = dst.col_stride();

					if (m == 1)
					{
						InterT mv = fun.transform(acc.get_scalar(0, acc.col_state(0)));
						for (index_t j = 1; j < n; ++j)
						{
							mv = fun.combine(mv, acc.get_scalar(0, acc.col_state(j)));
						}
					}
					else
					{
						mcol_impl_t::eval(fun, acc, m, n, pd, dst.row_stride(), dst.col_stride());
					}
				}
			}
			else
			{
				if (m > 0)
					fill(dst, fun.empty_value());
			}
		}
	};


} }

#endif
