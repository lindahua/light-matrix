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
	_reduce(Tag tag, const IMatrixXpr<Xpr, T>& a, linear_macc, scalar_ker)
	{
		typedef vec_reduce_core_setting<Tag, scalar_ker, T> setting_t;
		typedef typename setting_t::media_t media_t;

		typedef typename macc_accessor_map<Xpr, linear_macc, scalar_ker>::type accessor_t;
		typedef typename reduction_fun<Tag, scalar_ker, T>::type fun_t;
		typedef internal::vec_reduce_linear_scalar_term_accessor<
				fun_t, media_t, accessor_t> term_accessor_t;

		typedef vec_reduce_core<setting_t, ct_size<Xpr>::value> impl_t;

		fun_t fun(tag);
		accessor_t acc(a.derived());
		term_accessor_t tacc(fun, acc);
		const index_t n = a.nelems();

		if (n > 0)
		{
			return impl_t::reduce(fun, n, tacc);
		}
		else
		{
			return fun.empty_value();
		}
	}


	template<typename Tag, typename T1, class Xpr1, typename T2, class Xpr2>
	inline typename reduction_result<Tag, T1, T2>::type
	_reduce(Tag tag, const IMatrixXpr<Xpr1, T1>& a, const IMatrixXpr<Xpr2, T2>& b, linear_macc, scalar_ker)
	{
		typedef vec_reduce_core_setting<Tag, scalar_ker, T1, T2> setting_t;
		typedef typename setting_t::media_t media_t;

		typedef typename macc_accessor_map<Xpr1, linear_macc, scalar_ker>::type accessor1_t;
		typedef typename macc_accessor_map<Xpr2, linear_macc, scalar_ker>::type accessor2_t;
		typedef typename reduction_fun<Tag, scalar_ker, T1, T2>::type fun_t;
		typedef internal::vec_reduce_linear_scalar_term_accessor<
				fun_t, media_t, accessor1_t, accessor2_t> term_accessor_t;

		typedef vec_reduce_core<setting_t, common_ctsize<Xpr1, Xpr2>::value> impl_t;

		fun_t fun(tag);
		accessor1_t acc1(a.derived());
		accessor2_t acc2(b.derived());
		term_accessor_t tacc(fun, acc1, acc2);
		const index_t n = a.nelems();

		if (n > 0)
		{
			return impl_t::reduce(fun, n, tacc);
		}
		else
		{
			return fun.empty_value();
		}
	}


	template<typename Tag, typename T, class Xpr>
	inline typename reduction_result<Tag, T>::type
	_reduce(Tag tag, const IMatrixXpr<Xpr, T>& a, percol_macc, scalar_ker)
	{
		typedef vec_reduce_core_setting<Tag, scalar_ker, T> setting_t;
		typedef typename setting_t::media_t media_t;

		typedef typename macc_accessor_map<Xpr, percol_macc, scalar_ker>::type accessor_t;
		typedef typename reduction_fun<Tag, scalar_ker, T>::type fun_t;
		typedef internal::vec_reduce_percol_scalar_term_accessor<
				fun_t, media_t, accessor_t> term_accessor_t;

		typedef vec_reduce_core<setting_t, ct_rows<Xpr>::value> impl_t;

		fun_t fun(tag);
		accessor_t acc(a.derived());

		const index_t m = a.nrows();
		const index_t n = a.ncolumns();

		if (n > 0)
		{
			term_accessor_t tacc0(fun, 0, acc);
			media_t mv = impl_t::accumulate(fun, m, tacc0);

			if (n > 1)
			{
				for (index_t j = 1; j < n; ++j)
				{
					term_accessor_t tacc(fun, j, acc);
					media_t u = impl_t::accumulate(fun, m, tacc);
					mv = fun.combine(mv, u);
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
	_reduce(Tag tag, const IMatrixXpr<Xpr1, T1>& a, const IMatrixXpr<Xpr2, T2>& b, percol_macc, scalar_ker)
	{
		typedef vec_reduce_core_setting<Tag, scalar_ker, T1, T2> setting_t;
		typedef typename setting_t::media_t media_t;

		typedef typename macc_accessor_map<Xpr1, percol_macc, scalar_ker>::type accessor1_t;
		typedef typename macc_accessor_map<Xpr2, percol_macc, scalar_ker>::type accessor2_t;
		typedef typename reduction_fun<Tag, scalar_ker, T1, T2>::type fun_t;
		typedef internal::vec_reduce_percol_scalar_term_accessor<
				fun_t, media_t, accessor1_t, accessor2_t> term_accessor_t;

		typedef vec_reduce_core<setting_t, common_ctrows<Xpr1, Xpr2>::value> impl_t;

		fun_t fun(tag);
		accessor1_t acc1(a.derived());
		accessor2_t acc2(b.derived());

		const index_t m = get_common_nrows(a, b);
		const index_t n = get_common_ncolumns(a, b);

		if (n > 0)
		{
			term_accessor_t tacc0(fun, 0, acc1, acc2);
			media_t mv = impl_t::accumulate(fun, m, tacc0);

			if (n > 1)
			{
				for (index_t j = 1; j < n; ++j)
				{
					term_accessor_t tacc(fun, j, acc1, acc2);
					media_t u = impl_t::accumulate(fun, m, tacc);
					mv = fun.combine(mv, u);
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
