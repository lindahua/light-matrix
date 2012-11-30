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

#include "macc_eval_core.h"
#include "macc_reduce_core.h"

namespace lmat { namespace internal {

	template<typename KerCate> struct colwise_reduce_internal;
	template<typename InterT, typename KerCate> struct rowwise_reduce_internal;

	template<typename Tag, class Arg, class Dst>
	static void _colwise_reduce(const Tag& tag, const Arg& a, Dst& dst, scalar_kernel_t)
	{
		typedef typename matrix_traits<Arg>::value_type T;

		typedef typename reduction_result<Tag, T>::intermediate_type media_t;
		typedef typename reduction_fun<Tag, scalar_kernel_t, T>::type fun_t;

		typedef typename macc_accessor_map<Arg, percol_macc, scalar_kernel_t>::type accessor_t;
		typedef typename percol_macc_state_map<accessor_t>::type col_state_t;

		const index_t m = a.nrows();
		const index_t n = a.ncolumns();

		fun_t fun(tag);
		accessor_t acc(a);

		if (m > 0)
		{
			if (m == 1)
			{
				for (index_t j = 0; j < n; ++j)
				{
					media_t u = internal::reduc_single_term<Tag, T>::get(fun, acc, 0, 0);
					dst(0, j) = fun.get(u, 1);
				}
			}
			else
			{
				typedef vec_reduce<Tag, ct_rows<Arg>::value, scalar_kernel_t, T> impl_t;
				for (index_t j = 0; j < n; ++j)
				{
					media_t u = impl_t::eval_col(fun, m, acc, j);
					dst(0, j) = fun.get(u, m);
				}
			}
		}
		else
		{
			if (n > 0)
				fill(dst, fun.empty_value());
		}
	}


	template<typename Tag, class Arg1, class Arg2, class Dst>
	static void _colwise_reduce(const Tag& tag, const Arg1& a, class Arg2& b, Dst& dst, scalar_kernel_t)
	{
		typedef typename matrix_traits<Arg1>::value_type T1;
		typedef typename matrix_traits<Arg2>::value_type T2;

		typedef typename reduction_result<Tag, T1, T2>::intermediate_type media_t;
		typedef typename reduction_fun<Tag, scalar_kernel_t, T1, T2>::type fun_t;

		typedef typename macc_accessor_map<Arg1, percol_macc, scalar_kernel_t>::type accessor1_t;
		typedef typename macc_accessor_map<Arg2, percol_macc, scalar_kernel_t>::type accessor2_t;
		typedef typename percol_macc_state_map<accessor1_t>::type col_state1_t;
		typedef typename percol_macc_state_map<accessor2_t>::type col_state2_t;

		const index_t m = a.nrows();
		const index_t n = a.ncolumns();

		fun_t fun(tag);
		accessor1_t acc1(a);
		accessor2_t acc2(b);

		if (m > 0)
		{
			if (m == 1)
			{
				for (index_t j = 0; j < n; ++j)
				{
					media_t u = internal::reduc_single_term<Tag, T1, T2>::get(fun, acc1, acc2, 0, 0);
					dst(0, j) = fun.get(u, 1);
				}
			}
			else
			{
				typedef vec_reduce<Tag, common_ctrows<Arg1, Arg2>::value, scalar_kernel_t, T1, T2> impl_t;
				for (index_t j = 0; j < n; ++j)
				{
					media_t u = impl_t::eval_col(fun, m, acc1, acc2, j);
					dst(0, j) = fun.get(u, m);
				}
			}
		}
		else
		{
			if (n > 0)
				fill(dst, fun.empty_value());
		}
	}


	template<int M, typename Ker, typename InterT, typename ResultT>
	struct _rowwise_reduc_helper;

	template<int M, typename InterT, typename ResultT>
	struct _rowwise_reduc_helper<M, scalar_kernel_t, InterT, ResultT>
	{
		template<class Fun, class Acc>
		void eval(const Fun& fun, const Acc& acc, index_t m, index_t n, ResultT* dst, index_t dst_rs)
		{
			dense_col<InterT, M> temp(m);

			ContVecRW<scalar_kernel_t, InterT> out(temp.ptr_data());

			reduc_terms_accessor_pcol<InterT, Fun, Acc> tacc0(acc, 0);
			macc_vec_copy<scalar_kernel_t, M>::eval(tacc0, out);

			for (index_t j = 1; j < n; ++j)
			{
				reduc_terms_accessor_pcol<InterT, Fun, Acc> tacc(acc, j);
				vec_accum_core<InterT, M, scalar_kernel_t>::eval(out, tacc);
			}

			if (dst_rs == 1)
			{
				for (index_t i = 0; i < m; ++i)	dst[i] = fun.get(temp[i], n);
			}
			else
			{
				for (index_t i = 0; i < m; ++i)	dst[i * rs] = fun.get(temp[i], n);
			}
		}
	};





	template<typename Tag, class Arg, class Dst>
	static void _rowwise_reduce(const Tag& tag, const Arg& a, Dst& dst, scalar_kernel_t)
	{
		typedef typename matrix_traits<Arg>::value_type T;

		typedef typename reduction_result<Tag, T>::intermediate_type media_t;
		typedef typename reduction_result<Tag, T>::type result_t;
		typedef typename reduction_fun<Tag, scalar_kernel_t, T>::type fun_t;

		typedef typename macc_accessor_map<Arg, percol_macc, scalar_kernel_t>::type accessor_t;
		typedef typename percol_macc_state_map<accessor_t>::type col_state_t;

		typedef percol_to_linear_accessor<accessor_t, scalar_kernel_t, T> acc1_wrapper;

		const index_t M = common_ctrows<Arg, Dst>::value;
		const index_t m = a.nrows();
		const index_t n = a.ncolumns();

		fun_t fun(tag);
		accessor_t acc(a);
		result_t *pd = dst.ptr_data();
		const index_t rs = dst.row_stride();

		if (n == 1)
		{
			if (rs == 1)
			{
				ContVecRW<scalar_kernel_t, result_t> out(pd);
				internal::single_elem_reduc_core<Tag, M, scalar_kernel_t, T>::eval(fun, m, acc, out);
			}
			else
			{
				StepVecRW<scalar_kernel_t, result_t> out(pd, rs);
				internal::single_elem_reduc_core<Tag, M, scalar_kernel_t, T>::eval(fun, m, acc, out);
			}
		}
		else
		{
			_rowwise_reduc_helper<M, scalar_kernel_t, media_t, result_t>::eval(fun, acc, m, n, pd, rs);
		}
	}














} }

#endif
