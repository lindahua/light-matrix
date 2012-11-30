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
		TermT get_scalar(index_t i) const
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


	template<typename RT, int CTLen, typename KerCate>
	struct vec_reduce
	{
		template<class Fun, class Acc>
		LMAT_ENSURE_INLINE
		static RT eval(const Fun& fun, const index_t len, const Acc& a)
		{
			typedef reduc_terms_accessor<RT, Fun, Acc> term_accessor_t;
			term_accessor_t term_vec(fun, a);
			vec_reduce_core<RT, CTLen, KerCate>::eval(fun, len, term_vec);
		}

		template<class Fun, class Acc1, class Acc2>
		LMAT_ENSURE_INLINE
		static RT eval(const Fun& fun, const index_t len, const Acc1& a1, const Acc2& a2)
		{
			typedef reduc_terms_accessor<RT, Fun, Acc1, Acc2> term_accessor_t;
			term_accessor_t term_vec(fun, a1, a2);
			vec_reduce_core<RT, CTLen, KerCate>::eval(fun, len, term_vec);
		}
	};


} }

#endif
