/**
 * @file full_reduce.h
 *
 * Full matrix reduction expression
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_REDUCE_H_
#define LIGHTMAT_MATRIX_REDUCE_H_

#include <light_mat/matexpr/dense_accessors.h>
#include "bits/matrix_reduce_internal.h"

namespace lmat
{
	/********************************************
	 *
	 *  generic reduction functions
	 *
	 ********************************************/

	template<typename Op, typename T, class Xpr>
	inline typename unary_reduc_result<Op, T>::type
	reduce(Op op, const IMatrixXpr<Xpr, T>& xpr, scalar_kernel_t, linear_macc)
	{
		typedef typename macc_accessor_map<Xpr, scalar_kernel_t, linear_macc>::type accessor_t;
		typedef typename unary_reduc_result<Op, T>::type result_t;
		typedef typename unary_reduc_fun<Op, scalar_kernel_t, T>::type fun_t;

		fun_t fun(op);
		accessor_t acc(xpr.derived());
		const index_t n = xpr.nelems();

		if (n > 0)
		{
			return detail::vec_reduce<
					result_t, ct_size<Xpr>::value, scalar_kernel_t>::eval(fun, n, acc);
		}
		else
		{
			return fun.empty_value();
		}
	}


	template<typename Op, typename T1, class Xpr1, typename T2, class Xpr2>
	inline typename binary_reduc_result<Op, T1, T2>::type
	reduce(Op op, const IMatrixXpr<Xpr1, T1>& xpr1, const IMatrixXpr<Xpr2, T2>& xpr2, scalar_kernel_t, linear_macc)
	{
		typedef typename macc_accessor_map<Xpr1, scalar_kernel_t, linear_macc>::type accessor1_t;
		typedef typename macc_accessor_map<Xpr2, scalar_kernel_t, linear_macc>::type accessor2_t;

		typedef typename binary_reduc_result<Op, T1, T2>::type result_t;
		typedef typename binary_reduc_fun<Op, scalar_kernel_t, T1, T2>::type fun_t;

		fun_t fun(op);
		accessor1_t acc1(xpr1.derived());
		accessor1_t acc2(xpr2.derived());

		const index_t n = xpr1.nelems();

		if (n > 0)
		{
			return detail::vec_reduce<
					result_t, binary_ct_size<Xpr1, Xpr2>::value, scalar_kernel_t>::eval(fun, n, acc1, acc2);
		}
		else
		{
			return fun.empty_value();
		}
	}


	template<typename Op, typename T, class Xpr>
	inline typename unary_reduc_result<Op, T>::type
	reduce(Op op, const IMatrixXpr<Xpr, T>& xpr, scalar_kernel_t, percol_macc)
	{
		typedef typename macc_accessor_map<Xpr, scalar_kernel_t, percol_macc>::type accessor_t;
		typedef typename percol_macc_state_map<accessor_t>::type col_state_t;

		typedef typename unary_reduc_result<Op, T>::type result_t;
		typedef typename unary_reduc_fun<Op, scalar_kernel_t, T>::type fun_t;

		typedef detail::vec_reduce<result_t, ct_rows<Xpr>::value, scalar_kernel_t> impl_t;

		fun_t fun(op);
		accessor_t acc(xpr.derived());

		const index_t m = xpr.nrows();
		const index_t n = xpr.ncolumns();

		if (n > 0)
		{
			col_state_t s0 = acc.col_state(0);
			result_t r = impl_t::eval_s(fun, m, acc, s0);

			if (n > 1)
			{
				for (index_t j = 1; j < n; ++j)
				{
					col_state_t s = acc.col_state(j);
					r = fun.combine(r, impl_t::eval_s(fun, m, acc, s));
				}
			}

			return r;
		}
		else
		{
			return fun.empty_value();
		}
	}


	template<typename Op, typename T1, class Xpr1, typename T2, class Xpr2>
	inline typename binary_reduc_result<Op, T1, T2>::type
	reduce(Op op, const IMatrixXpr<Xpr1, T1>& xpr1, class IMatrixXpr<Xpr2, T2>& xpr2, scalar_kernel_t, percol_macc)
	{
		typedef typename macc_accessor_map<Xpr1, scalar_kernel_t, percol_macc>::type accessor1_t;
		typedef typename percol_macc_state_map<accessor1_t>::type col_state1_t;
		typedef typename macc_accessor_map<Xpr2, scalar_kernel_t, percol_macc>::type accessor2_t;
		typedef typename percol_macc_state_map<accessor2_t>::type col_state2_t;

		typedef typename binary_reduc_result<Op, T1, T2>::type result_t;
		typedef typename binary_reduc_fun<Op, scalar_kernel_t, T1, T2>::type fun_t;

		typedef detail::vec_reduce<result_t, binary_ct_rows<Xpr1, Xpr2>::value, scalar_kernel_t> impl_t;

		fun_t fun(op);
		accessor1_t acc1(xpr1.derived());
		accessor1_t acc2(xpr2.derived());

		const index_t m = xpr1.nrows();
		const index_t n = xpr2.ncolumns();

		if (n > 0)
		{
			col_state1_t s0_1 = acc1.col_state(0);
			col_state1_t s0_2 = acc2.col_state(0);
			result_t r = impl_t::eval_s(fun, m, acc1, s0_1, acc2, s0_2);

			if (n > 1)
			{
				for (index_t j = 1; j < n; ++j)
				{
					col_state1_t s1 = acc1.col_state(j);
					col_state1_t s2 = acc2.col_state(j);

					r = fun.combine(r, impl_t::eval_s(fun, m, acc1, s1, acc2, s2));
				}
			}

			return r;
		}
		else
		{
			return fun.empty_value();
		}
	}


	template<typename Op, typename T, class Xpr>
	LMAT_ENSURE_INLINE
	inline typename unary_reduc_result<Op, T>::type
	reduce(Op op, const IMatrixXpr<Xpr, T>& xpr)
	{
		typedef dense_matrix<T, ct_rows<Xpr>::value, ct_cols<Xpr>::value> dmat_t;
		typedef default_macc_scheme<Xpr, dmat_t> scheme_t;

		typedef typename scheme_t::kernel_category ker_t;
		typedef typename scheme_t::access_category acc_t;

		return reduce(op, xpr, ker_t(), acc_t());
	}


}

#endif /* MATRIX_REDUC_EXPR_H_ */


