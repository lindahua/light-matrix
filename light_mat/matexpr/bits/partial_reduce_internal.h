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


	/********************************************
	 *
	 *  col-wise reduction
	 *
	 ********************************************/

	template<typename Tag, class Xpr, typename T, class Dst>
	static void _partial_reduce(const Tag& tag, colwise,
			const IMatrixXpr<Xpr, T>& a, Dst& dst, scalar_ker)
	{
		typedef vec_reduce_core_setting<Tag, scalar_ker, T> setting_t;
		typedef typename setting_t::media_t media_t;

		typedef typename macc_accessor_map<Xpr, linear_macc, scalar_ker>::type accessor_t;
		typedef typename reduction_fun<Tag, scalar_ker, T>::type fun_t;
		typedef internal::vec_reduce_percol_scalar_term_accessor<
				fun_t, media_t, accessor_t> term_accessor_t;

		typedef vec_reduce_core<setting_t, ct_rows<Xpr>::value> impl_t;

		fun_t fun(tag);
		accessor_t acc(a.derived());

		const index_t m = a.nrows();
		const index_t n = a.ncolumns();

		if (m > 0)
		{
			if (m == 1)
			{
				for (index_t j = 0; j < n; ++j)
				{
					media_t u = fun.transform(acc.get_scalar(0, acc.col_state(0)));
					dst(0, j) = fun.get(u, 1);
				}
			}
			else
			{
				for (index_t j = 0; j < n; ++j)
				{
					term_accessor_t tacc(fun, j, acc);
					dst(0, j) = impl_t::reduce(fun, m, tacc);
				}
			}
		}
		else
		{
			if (n > 0)
				fill(dst, fun.empty_value());
		}
	}


	template<typename Tag, class Xpr1, typename T1, class Xpr2, typename T2, class Dst>
	static void _partial_reduce(const Tag& tag, colwise,
			const IMatrixXpr<Xpr1, T1>& a, class IMatrixXpr<Xpr2, T2>& b, Dst& dst, scalar_ker)
	{
		typedef vec_reduce_core_setting<Tag, scalar_ker, T1, T2> setting_t;
		typedef typename setting_t::media_t media_t;

		typedef typename macc_accessor_map<Xpr1, linear_macc, scalar_ker>::type accessor1_t;
		typedef typename macc_accessor_map<Xpr2, linear_macc, scalar_ker>::type accessor2_t;
		typedef typename reduction_fun<Tag, scalar_ker, T1, T2>::type fun_t;
		typedef internal::vec_reduce_percol_scalar_term_accessor<
				fun_t, media_t, accessor1_t, accessor2_t> term_accessor_t;

		typedef vec_reduce_core<setting_t, common_ctrows<Xpr1, Xpr2>::value> impl_t;

		fun_t fun(tag);
		accessor1_t acc1(a.derived());
		accessor2_t acc2(b.derived());

		const index_t m = a.nrows();
		const index_t n = a.ncolumns();

		if (m > 0)
		{
			if (m == 1)
			{
				for (index_t j = 0; j < n; ++j)
				{
					media_t u = fun.transform(
							acc1.get_scalar(0, acc1.col_state(0)),
							acc2.get_scalar(0, acc2.col_state(0)));
					dst(0, j) = fun.get(u, 1);
				}
			}
			else
			{
				for (index_t j = 0; j < n; ++j)
				{
					term_accessor_t tacc(fun, j, acc1, acc2);
					dst(0, j) = impl_t::reduce(fun, m, tacc);
				}
			}
		}
		else
		{
			if (n > 0)
				fill(dst, fun.empty_value());
		}
	}



	/********************************************
	 *
	 *  row-wise reduction
	 *
	 ********************************************/

	template<class Impl, int M, typename Ker, typename InterT, typename ResultT>
	struct _rowwise_reduc_helper;

	template<class Fun, class TAcc, typename ResultT>
	LMAT_ENSURE_INLINE
	inline void _rowwise_single_reduce(const Fun& fun, index_t m, const TAcc& acc, ResultT *dst, index_t dst_rs)
	{
		if (dst_rs == 1)
		{
			for (index_t i = 0; i < m; ++i)
			{
				dst[i] = fun.get(acc.get_scalar(i), 1);
			}
		}
		else
		{
			for (index_t i = 0; i < m; ++i)
			{
				dst[i * dst_rs] = fun.get(acc.get_scalar(i), 1);
			}
		}
	}



	template<class Fun, class In, typename ResultT>
	LMAT_ENSURE_INLINE
	inline void _rowwise_post_reduce(const Fun& fun, index_t len, const In& in, ResultT* dst, index_t dst_rs)
	{
		if (dst_rs == 1)
		{
			for (index_t i = 0; i < len; ++i) dst[i] = fun.get(in[i]);
		}
		else
		{
			for (index_t i = 0; i < len; ++i) dst[i * dst_rs] = fun.get(in[i]);
		}
	}


	template<class Impl, int M, typename InterT, typename ResultT>
	struct _rowwise_reduc_helper<Impl, M, scalar_ker, InterT, ResultT>
	{
		template<class Fun, class Acc>
		void eval(const Fun& fun, const Acc& acc, index_t m, index_t n, ResultT* dst, index_t dst_rs)
		{
			typedef internal::vec_reduce_percol_scalar_term_accessor<Fun, InterT, Acc> term_accessor_t;
			term_accessor_t tacc0(fun, 0, acc);

			if (n == 1)
			{
				_rowwise_single_reduce(fun, m, tacc0, dst, dst_rs);
			}
			else
			{
				dense_col<InterT, M> temp(m);
				ContVecRW<scalar_ker, InterT> out(temp.ptr_data());

				Impl::init_terms(fun, m, tacc0, out);

				for (index_t j = 1; j < n; ++j)
				{
					term_accessor_t tacc(fun, j, acc);
					Impl::parallel_accum(fun, m, tacc, out);
				}

				_rowwise_post_reduce(fun, m, temp, dst, dst_rs);
			}
		}

		template<class Fun, class Acc1, class Acc2>
		void eval(const Fun& fun, const Acc1& acc1, const Acc2& acc2, index_t m, index_t n, ResultT* dst, index_t dst_rs)
		{
			typedef internal::vec_reduce_percol_scalar_term_accessor<Fun, InterT, Acc1, Acc2> term_accessor_t;
			term_accessor_t tacc0(fun, 0, acc1, acc2);

			if (n == 1)
			{
				_rowwise_single_reduce(fun, m, tacc0, dst, dst_rs);
			}
			else
			{
				dense_col<InterT, M> temp(m);
				ContVecRW<scalar_ker, InterT> out(temp.ptr_data());

				Impl::init_terms(fun, m, tacc0, out);

				for (index_t j = 1; j < n; ++j)
				{
					term_accessor_t tacc(fun, j, acc1, acc2);
					Impl::parallel_accum(fun, m, tacc, out);
				}

				_rowwise_post_reduce(fun, m, temp, dst, dst_rs);
			}
		}
	};

	template<class Impl, int M, typename RT>
	struct _rowwise_reduc_helper<Impl, M, scalar_ker, RT, RT>
	{
		template<class Fun, class Acc>
		void eval(const Fun& fun, const Acc& acc, index_t m, index_t n, RT* dst, index_t dst_rs)
		{
			typedef internal::vec_reduce_percol_scalar_term_accessor<Fun, RT, Acc> term_accessor_t;
			term_accessor_t tacc0(fun, 0, acc);

			if (dst_rs == 1)
			{
				if (n == 1)
				{
					_rowwise_single_reduce(fun, m, tacc0, dst, dst_rs);
				}
				else
				{
					ContVecRW<scalar_ker, RT> out(dst);
					Impl::init_terms(fun, m, tacc0, out);

					for (index_t j = 1; j < n; ++j)
					{
						term_accessor_t tacc(fun, j, acc);
						Impl::parallel_accum(fun, m, tacc, out);
					}

					if (Impl::setting_t::with_post)
					{
						_rowwise_post_reduce(fun, m, out, dst, 1);
					}
				}
			}
			else
			{
				if (n == 1)
				{
					_rowwise_single_reduce(fun, m, tacc0, dst, dst_rs);
				}
				else
				{
					StepVecRW<scalar_ker, RT> out(dst, dst_rs);

					term_accessor_t tacc(fun, 0, acc);
					Impl::init_terms(fun, m, tacc, out);

					for (index_t j = 1; j < n; ++j)
					{
						term_accessor_t tacc(fun, j, acc);
						Impl::parallel_accum(fun, m, tacc, out);
					}

					if (Impl::setting_t::with_post)
					{
						_rowwise_post_reduce(fun, m, out, dst, dst_rs);
					}
				}
			}
		}


		template<class Fun, class Acc1, class Acc2>
		void eval(const Fun& fun, const Acc1& acc1, const Acc2& acc2, index_t m, index_t n, RT* dst, index_t dst_rs)
		{
			typedef internal::vec_reduce_percol_scalar_term_accessor<Fun, RT, Acc1, Acc2> term_accessor_t;
			term_accessor_t tacc0(fun, 0, acc1, acc2);

			if (dst_rs == 1)
			{
				if (n == 1)
				{
					_rowwise_single_reduce(fun, m, tacc0, dst, dst_rs);
				}
				else
				{
					ContVecRW<scalar_ker, RT> out(dst);
					Impl::init_terms(fun, m, tacc0, out);

					for (index_t j = 1; j < n; ++j)
					{
						term_accessor_t tacc(fun, j, acc1, acc2);
						Impl::parallel_accum(fun, m, tacc, out);
					}

					if (Impl::setting_t::with_post)
					{
						_rowwise_post_reduce(fun, m, out, dst, 1);
					}
				}
			}
			else
			{
				if (n == 1)
				{
					_rowwise_single_reduce(fun, m, tacc0, dst, dst_rs);
				}
				else
				{
					StepVecRW<scalar_ker, RT> out(dst, dst_rs);

					term_accessor_t tacc(fun, 0, acc1, acc2);
					Impl::init_terms(fun, m, tacc, out);

					for (index_t j = 1; j < n; ++j)
					{
						term_accessor_t tacc(fun, j, acc1, acc2);
						Impl::parallel_accum(fun, m, tacc, out);
					}

					if (Impl::setting_t::with_post)
					{
						_rowwise_post_reduce(fun, m, out, dst, dst_rs);
					}
				}
			}
		}

	};


	template<typename Tag, class Xpr, typename T, class Dst>
	inline void _partial_reduce(const Tag& tag, rowwise,
			const IMatrixXpr<Xpr, T>& a, Dst& dst, scalar_ker)
	{
		typedef vec_reduce_core_setting<Tag, scalar_ker, T> setting_t;
		typedef typename setting_t::media_t media_t;
		typedef typename setting_t::result_t result_t;

		typedef typename macc_accessor_map<Xpr, linear_macc, scalar_ker>::type accessor_t;
		typedef typename reduction_fun<Tag, scalar_ker, T>::type fun_t;

		const index_t m = a.nrows();
		const index_t n = a.ncolumns();
		accessor_t acc(a);

		typedef vec_reduce_core<setting_t, ct_rows<Xpr>::value> impl_t;
		typedef _rowwise_reduc_helper<impl_t, ct_rows<Xpr>::value, scalar_ker, media_t, result_t> helper_t;

		helper_t::eval(fun_t(tag), acc, m, n, dst.ptr_data(), dst.row_stride());
	}

	template<typename Tag, class Xpr1, typename T1, class Xpr2, typename T2, class Dst>
	inline void _partial_reduce(const Tag& tag, rowwise,
			const IMatrixXpr<Xpr1, T1>& a, const IMatrixXpr<Xpr2, T2>& b, Dst& dst, scalar_ker)
	{
		typedef vec_reduce_core_setting<Tag, scalar_ker, T1, T2> setting_t;
		typedef typename setting_t::media_t media_t;
		typedef typename setting_t::result_t result_t;

		typedef typename macc_accessor_map<Xpr1, linear_macc, scalar_ker>::type accessor1_t;
		typedef typename macc_accessor_map<Xpr2, linear_macc, scalar_ker>::type accessor2_t;
		typedef typename reduction_fun<Tag, scalar_ker, T1, T2>::type fun_t;

		const index_t m = get_common_nrows(a, b);
		const index_t n = get_common_ncolumns(a, b);
		accessor1_t acc1(a);
		accessor2_t acc2(b);

		typedef vec_reduce_core<setting_t, common_ctrows<Xpr1, Xpr2>::value> impl_t;
		typedef _rowwise_reduc_helper<impl_t, common_ctrows<Xpr1, Xpr2>::value, scalar_ker, media_t, result_t> helper_t;

		helper_t::eval(fun_t(tag), acc1, acc2, m, n, dst.ptr_data(), dst.row_stride());
	}

} }

#endif
