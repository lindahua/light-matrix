/**
 * @file macc_reduce_core.h
 *
 * Internal implementation of access-based reduction
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MACC_REDUCE_CORE_H_
#define LIGHTMAT_MACC_REDUCE_CORE_H_

#include <light_mat/matrix/matrix_properties.h>
#include <light_mat/math/reduction_functors.h>

namespace lmat { namespace internal {


	/********************************************
	 *
	 *  vector access
	 *
	 ********************************************/

	template<typename TermT, class Fun, class Acc1, class Acc2=nil_t>
	struct reduc_terms_accessor
	{
		const Fun& fun;
		const Acc1& acc1;
		const Acc2& acc2;

		LMAT_ENSURE_INLINE
		reduc_terms_accessor(const Fun& f, const Acc1& a1, const Acc2& a2)
		: fun(f), acc1(a1), acc2(a2)
		{ }

		LMAT_ENSURE_INLINE
		TermT get_term(index_t i) const
		{
			return fun.transform(acc1.get_scalar(i), acc2.get_scalar(i));
		}
	};


	template<typename TermT, class Fun, class Acc1>
	struct reduc_terms_accessor<TermT, Fun, Acc1, nil_t>
	{
		const Fun& fun;
		const Acc1& acc1;

		LMAT_ENSURE_INLINE
		reduc_terms_accessor(const Fun& f, const Acc1& a1)
		: fun(f), acc1(a1)
		{ }

		LMAT_ENSURE_INLINE
		TermT get_term(index_t i) const
		{
			return fun.transform(acc1.get_scalar(i));
		}
	};




	/********************************************
	 *
	 *  vector reduction
	 *
	 ********************************************/

	template<typename RT, int CTLen, typename KerCate> struct vec_reduce_core;

	template<typename RT, int CTLen>
	struct vec_reduce_core<RT, CTLen, scalar_kernel_t>
	{
		template<class Fun, class Vec>
		LMAT_ENSURE_INLINE
		static RT eval(const Fun& fun, const index_t, const Vec& vec)
		{
			RT r = vec.get_term(0);
			for (index_t i = 1; i < CTLen; ++i)
			{
				r = fun.combine(r, vec.get_term(i));
			}
			return r;
		}
	};


	template<typename RT>
	struct vec_reduce_core<RT, 0, scalar_kernel_t>
	{
		template<class Fun, class Vec>
		LMAT_ENSURE_INLINE
		static RT eval(const Fun& fun, const index_t len, const Vec& vec)
		{
			RT r = vec.get_term(0);
			for (index_t i = 1; i < len; ++i)
			{
				r = fun.combine(r, vec.get_term(i));
			}
			return r;
		}
	};


	template<typename RT>
	struct vec_reduce_core<RT, 1, scalar_kernel_t>
	{
		template<class Fun, class Vec>
		LMAT_ENSURE_INLINE
		static RT eval(const Fun& fun, const index_t len, const Vec& vec)
		{
			return vec.get_term(0);
		}
	};


	template<typename Tag, int CTLen, typename KerCate, typename T1, typename T2=nil_t> struct vec_reduce;


	template<typename Tag, int CTLen, typename T>
	struct vec_reduce<Tag, CTLen, scalar_kernel_t, T, nil_t>
	{
		typedef typename reduction_result<Tag, T>::intermediate_type RT;
		typedef typename reduction_fun<Tag, scalar_kernel_t, T>::type Fun;

		template<class Acc>
		LMAT_ENSURE_INLINE
		static RT eval(const Fun& fun, const index_t len, const Acc& a)
		{
			typedef reduc_terms_accessor<RT, Fun, Acc> term_accessor_t;
			term_accessor_t term_vec(fun, a);
			return vec_reduce_core<RT, CTLen, scalar_kernel_t>::eval(fun, len, term_vec);
		}

		template<class PerColAcc>
		LMAT_ENSURE_INLINE
		static RT eval_col(const Fun& fun, const index_t len, const PerColAcc& a, const index_t j)
		{
			typedef percol_to_linear_accessor<PerColAcc, scalar_kernel_t, T> Acc;
			return eval(fun, len, Acc(a, a.col_state(j)));
		}
	};


	template<typename Tag, int CTLen, typename T1, typename T2>
	struct vec_reduce<Tag, CTLen, scalar_kernel_t, T1, T2>
	{
		typedef typename reduction_result<Tag, T1, T2>::intermediate_type RT;
		typedef typename reduction_fun<Tag, scalar_kernel_t, T1, T2>::type Fun;

		template<class Acc1, class Acc2>
		LMAT_ENSURE_INLINE
		static RT eval(const Fun& fun, const index_t len, const Acc1& a1, const Acc2& a2)
		{
			typedef reduc_terms_accessor<RT, Fun, Acc1, Acc2> term_accessor_t;
			term_accessor_t term_vec(fun, a1, a2);
			return vec_reduce_core<RT, CTLen, scalar_kernel_t>::eval(fun, len, term_vec);
		}

		template<class PerColAcc1, class PerColAcc2>
		LMAT_ENSURE_INLINE
		static RT eval_col(const Fun& fun, const index_t len, const PerColAcc1& a1, const PerColAcc2& a2, const index_t j)
		{
			typedef percol_to_linear_accessor<PerColAcc1, scalar_kernel_t, T1> Acc1;
			typedef percol_to_linear_accessor<PerColAcc2, scalar_kernel_t, T2> Acc2;
			return eval(fun, len,
					Acc1(a1, a1.col_state(j)),
					Acc2(a2, a2.col_state(j)) );
		}
	};



} }

#endif
