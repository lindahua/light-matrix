/**
 * @file full_reduce_internal.h
 *
 * @brief Internal implementation of full matrix reduction
 *
 * @author Dahua Lin
 */

#ifndef LIGHTMAT_FULL_REDUCE_INTERNAL_H_
#define LIGHTMAT_FULL_REDUCE_INTERNAL_H_

#include <light_mat/matexpr/dense_accessors.h>
#include "macc_reduce_core.h"

namespace lmat { namespace internal {

	template<typename Tag, typename T, class Xpr>
	inline typename reduction_result<Tag, T>::type
	_reduce(Tag tag, const IMatrixXpr<Xpr, T>& a, linear_macc, scalar_kernel_t)
	{
		typedef typename macc_accessor_map<Xpr, linear_macc, scalar_kernel_t>::type accessor_t;
		typedef typename reduction_result<Tag, T>::intermediate_type media_t;
		typedef typename reduction_fun<Tag, scalar_kernel_t, T>::type fun_t;

		typedef internal::vec_reduce<Tag, ct_size<Xpr>::value, scalar_kernel_t, T> impl_t;

		fun_t fun(tag);
		accessor_t acc(a.derived());
		const index_t n = a.nelems();

		if (n > 0)
		{
			media_t mv = impl_t::eval(fun, n, acc);
			return fun.get(mv, n);
		}
		else
		{
			return fun.empty_value();
		}
	}


	template<typename Tag, typename T1, class Xpr1, typename T2, class Xpr2>
	inline typename reduction_result<Tag, T1, T2>::type
	_reduce(Tag tag, const IMatrixXpr<Xpr1, T1>& a, const IMatrixXpr<Xpr2, T2>& b, linear_macc, scalar_kernel_t)
	{
		typedef typename macc_accessor_map<Xpr1, linear_macc, scalar_kernel_t>::type accessor1_t;
		typedef typename macc_accessor_map<Xpr2, linear_macc, scalar_kernel_t>::type accessor2_t;

		typedef typename reduction_result<Tag, T1, T2>::intermediate_type media_t;
		typedef typename reduction_fun<Tag, scalar_kernel_t, T1, T2>::type fun_t;

		typedef internal::vec_reduce<Tag, common_ctsize<Xpr1, Xpr2>::value, scalar_kernel_t, T1, T2> impl_t;

		fun_t fun(tag);
		accessor1_t acc1(a.derived());
		accessor1_t acc2(b.derived());

		const index_t n = a.nelems();

		if (n > 0)
		{
			media_t mv = impl_t::eval(fun, n, acc1, acc2);
			return fun.get(mv, n);
		}
		else
		{
			return fun.empty_value();
		}
	}


	template<typename Tag, typename T, class Xpr>
	inline typename reduction_result<Tag, T>::type
	_reduce(Tag tag, const IMatrixXpr<Xpr, T>& a, percol_macc, scalar_kernel_t)
	{
		typedef typename macc_accessor_map<Xpr, percol_macc, scalar_kernel_t>::type accessor_t;
		typedef typename percol_macc_state_map<accessor_t>::type col_state_t;

		typedef typename reduction_result<Tag, T>::intermediate_type media_t;
		typedef typename reduction_fun<Tag, scalar_kernel_t, T>::type fun_t;

		typedef percol_to_linear_accessor<accessor_t, scalar_kernel_t, T> acc_wrapper;
		typedef internal::vec_reduce<Tag, ct_rows<Xpr>::value, scalar_kernel_t, T> impl_t;

		fun_t fun(tag);
		accessor_t acc(a.derived());

		const index_t m = a.nrows();
		const index_t n = a.ncolumns();

		if (n > 0)
		{
			media_t mv = impl_t::eval_col(fun, m, acc, 0);

			if (n > 1)
			{
				for (index_t j = 1; j < n; ++j)
					mv = fun.combine(mv, impl_t::eval_col(fun, m, acc, j));
			}

			return fun.get(mv, a.nelems());
		}
		else
		{
			return fun.empty_value();
		}
	}


	template<typename Tag, typename T1, class Xpr1, typename T2, class Xpr2>
	inline typename reduction_result<Tag, T1, T2>::type
	_reduce(Tag tag, const IMatrixXpr<Xpr1, T1>& a, const IMatrixXpr<Xpr2, T2>& b, percol_macc, scalar_kernel_t)
	{
		typedef typename macc_accessor_map<Xpr1, percol_macc, scalar_kernel_t>::type accessor1_t;
		typedef typename percol_macc_state_map<accessor1_t>::type col_state1_t;
		typedef typename macc_accessor_map<Xpr2, percol_macc, scalar_kernel_t>::type accessor2_t;
		typedef typename percol_macc_state_map<accessor2_t>::type col_state2_t;

		typedef typename reduction_result<Tag, T1, T2>::intermediate_type media_t;
		typedef typename reduction_fun<Tag, scalar_kernel_t, T1, T2>::type fun_t;

		typedef percol_to_linear_accessor<accessor1_t, scalar_kernel_t, T1> acc1_wrapper;
		typedef percol_to_linear_accessor<accessor2_t, scalar_kernel_t, T2> acc2_wrapper;
		typedef internal::vec_reduce<Tag, common_ctrows<Xpr1, Xpr2>::value, scalar_kernel_t, T1, T2> impl_t;

		fun_t fun(tag);
		accessor1_t acc1(a.derived());
		accessor1_t acc2(b.derived());

		const index_t m = a.nrows();
		const index_t n = a.ncolumns();

		if (n > 0)
		{
			media_t mv = impl_t::eval_col(fun, m, acc1, acc2, 0);

			if (n > 1)
			{
				for (index_t j = 1; j < n; ++j)
					mv = fun.combine(mv, impl_t::eval_col(fun, m, acc1, acc2, j));
			}

			return fun.get(mv, a.nelems());
		}
		else
		{
			return fun.empty_value();
		}
	}

} }

#endif /* FULL_REDUCE_INTERNAL_H_ */
