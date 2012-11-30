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

		fun_t fun(tag);
		accessor_t acc(a.derived());
		const index_t n = a.nelems();

		if (n > 0)
		{
			media_t mv = internal::vec_reduce<
					media_t, ct_size<Xpr>::value, scalar_kernel_t>::eval(fun, n, acc);
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

		fun_t fun(tag);
		accessor1_t acc1(a.derived());
		accessor1_t acc2(b.derived());

		const index_t n = a.nelems();

		if (n > 0)
		{
			media_t mv = internal::vec_reduce<
					media_t, common_ctsize<Xpr1, Xpr2>::value, scalar_kernel_t>::eval(fun, n, acc1, acc2);
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
		typedef internal::vec_reduce<media_t, ct_rows<Xpr>::value, scalar_kernel_t> impl_t;

		fun_t fun(tag);
		accessor_t acc(a.derived());

		const index_t m = a.nrows();
		const index_t n = a.ncolumns();

		if (n > 0)
		{
			col_state_t s0 = acc.col_state(0);
			media_t mv = impl_t::eval(fun, m, acc_wrapper(acc, s0));

			if (n > 1)
			{
				for (index_t j = 1; j < n; ++j)
				{
					col_state_t s = acc.col_state(j);
					mv = fun.combine(mv, impl_t::eval(fun, m, acc_wrapper(acc, s)));
				}
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
		typedef internal::vec_reduce<media_t, common_ctrows<Xpr1, Xpr2>::value, scalar_kernel_t> impl_t;

		fun_t fun(tag);
		accessor1_t acc1(a.derived());
		accessor1_t acc2(b.derived());

		const index_t m = a.nrows();
		const index_t n = a.ncolumns();

		if (n > 0)
		{
			col_state1_t s0_1 = acc1.col_state(0);
			col_state2_t s0_2 = acc2.col_state(0);
			media_t mv = impl_t::eval(fun, m, acc1_wrapper(acc1, s0_1), acc2_wrapper(acc2, s0_2));

			if (n > 1)
			{
				for (index_t j = 1; j < n; ++j)
				{
					col_state1_t s1 = acc1.col_state(j);
					col_state1_t s2 = acc2.col_state(j);

					mv = fun.combine(mv, impl_t::eval(fun, m,
							acc1_wrapper(acc1, s1), acc2_wrapper(acc2, s2) ));
				}
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
