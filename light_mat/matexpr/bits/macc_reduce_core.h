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
	 *  term accessors
	 *
	 ********************************************/

	template<class Fun, typename RT, class Acc1, class Acc2=nil_t>
	struct vec_reduce_linear_scalar_term_accessor;

	template<class Fun, typename RT, class Acc1, class Acc2=nil_t>
	struct vec_reduce_percol_scalar_term_accessor;

	template<class Fun, typename RT, class Acc>
	struct vec_reduce_linear_scalar_term_accessor<Fun, RT, Acc, nil_t>
	{
		const Fun& m_fun;
		const Acc& m_acc;

		LMAT_ENSURE_INLINE
		vec_reduce_linear_scalar_term_accessor(const Fun& fun, const Acc& acc)
		: m_fun(fun), m_acc(acc)
		{ }

		LMAT_ENSURE_INLINE
		RT get_scalar(index_t i) const
		{
			return m_fun.transform(m_acc.get_scalar(i));
		}
	};


	template<class Fun, typename RT, class Acc1, class Acc2>
	struct vec_reduce_linear_scalar_term_accessor
	{
		const Fun& m_fun;
		const Acc1& m_acc1;
		const Acc2& m_acc2;

		LMAT_ENSURE_INLINE
		vec_reduce_linear_scalar_term_accessor(const Fun& fun, const Acc1& acc1, const Acc2& acc2)
		: m_fun(fun), m_acc1(acc1), m_acc2(acc2)
		{ }

		LMAT_ENSURE_INLINE
		RT get_scalar(index_t i) const
		{
			return m_fun.transform(m_acc1.get_scalar(i), m_acc2.get_scalar(i));
		}
	};


	template<class Fun, typename RT, class Acc>
	struct vec_reduce_percol_scalar_term_accessor<Fun, RT, Acc, nil_t>
	{
		const Fun& m_fun;
		const Acc& m_acc;
		typename percol_macc_state_map<Acc>::type m_colstate;

		LMAT_ENSURE_INLINE
		vec_reduce_percol_scalar_term_accessor(const Fun& fun, index_t j, const Acc& acc)
		: m_fun(fun), m_acc(acc), m_colstate(acc.col_state(j))
		{ }

		LMAT_ENSURE_INLINE
		RT get_scalar(index_t i) const
		{
			return m_fun.transform(m_acc.get_scalar(i, m_colstate));
		}
	};


	template<class Fun, typename RT, class Acc1, class Acc2>
	struct vec_reduce_percol_scalar_term_accessor
	{
		const Fun& m_fun;
		const Acc1& m_acc1;
		const Acc2& m_acc2;
		typename percol_macc_state_map<Acc1>::type m_colstate1;
		typename percol_macc_state_map<Acc2>::type m_colstate2;

		LMAT_ENSURE_INLINE
		vec_reduce_percol_scalar_term_accessor(const Fun& fun, index_t j, const Acc1& acc1, const Acc2& acc2)
		: m_fun(fun), m_acc1(acc1), m_acc2(acc2), m_colstate1(acc1.col_state(j)), m_colstate2(acc2.col_state(j))
		{ }

		LMAT_ENSURE_INLINE
		RT get_scalar(index_t i) const
		{
			return m_fun.transform(m_acc1.get_scalar(i, m_colstate1), m_acc2.get_scalar(i, m_colstate2));
		}
	};


	/********************************************
	 *
	 *  core classes
	 *
	 ********************************************/

	template<typename Tag, typename Ker, typename T1, typename T2=nil_t>
	struct vec_reduce_core_setting;

	template<typename Tag, typename T1, typename T2>
	struct vec_reduce_core_setting<Tag, scalar_ker, T1, T2>
	{
		typedef Tag tag_t;
		typedef scalar_ker ker_t;

		typedef typename reduction_result<Tag, T1, T2>::type result_t;
		typedef typename reduction_result<Tag, T1, T2>::intermediate_type media_t;
		typedef typename reduction_fun<Tag, scalar_ker, T1, T2>::type fun_t;

		static const bool with_post = reduction_with_post<Tag>::value;
	};



	template<class Setting, int Len>
	struct vec_reduce_core;

	template<typename Tag, int Len, typename T1, typename T2>
	struct vec_reduce_core<vec_reduce_core_setting<Tag, scalar_ker, T1, T2>, Len>
	{
		typedef vec_reduce_core_setting<Tag, scalar_ker, T1, T2> setting_t;

		typedef typename setting_t::result_t result_t;
		typedef typename setting_t::media_t media_t;
		typedef typename setting_t::fun_t fun_t;

		template<class TAcc>
		LMAT_ENSURE_INLINE
		static media_t accumulate(const fun_t& fun, index_t, const TAcc& tacc)
		{
			media_t u = tacc.get_scalar(0);

			for (index_t i = 1; i < Len; ++i)
			{
				u = fun.combine(u, tacc.get_scalar(i));
			}

			return u;
		}

		template<class TAcc>
		LMAT_ENSURE_INLINE
		static result_t reduce(const fun_t& fun, index_t, const TAcc& tacc)
		{
			return fun.get(accumulate(fun, Len, tacc), Len);
		}


		template<class TAcc, class Out>
		LMAT_ENSURE_INLINE
		static void parallel_accum(const fun_t& fun, index_t, const TAcc& tacc, Out& out)
		{
			for (index_t i = 0; i < Len; ++i)
			{
				out[i] = fun.combine(out[i], tacc.get_scalar(i));
			}
		}

		template<class TAcc, class Out>
		LMAT_ENSURE_INLINE
		static void init_terms(const fun_t& fun, index_t, const TAcc& tacc, Out& out)
		{
			for (index_t i = 0; i < Len; ++i)
			{
				out[i] = tacc.get_scalar(i);
			}
		}

		template<class TAcc, class Out>
		LMAT_ENSURE_INLINE
		static void single_reduce(const fun_t& fun, index_t, const TAcc& tacc, Out& out)
		{
			for (index_t i = 0; i < Len; ++i)
			{
				out[i] = fun.get(tacc.get_scalar(i), 1);
			}
		}
	};


	template<typename Tag, typename T1, typename T2>
	struct vec_reduce_core<vec_reduce_core_setting<Tag, scalar_ker, T1, T2>, 0>
	{
		typedef vec_reduce_core_setting<Tag, scalar_ker, T1, T2> setting_t;

		typedef typename setting_t::result_t result_t;
		typedef typename setting_t::media_t media_t;
		typedef typename setting_t::fun_t fun_t;

		template<class TAcc>
		LMAT_ENSURE_INLINE
		static media_t accumulate(const fun_t& fun, index_t len, const TAcc& tacc)
		{
			media_t u = tacc.get_scalar(0);

			for (index_t i = 1; i < len; ++i)
			{
				u = fun.combine(u, tacc.get_scalar(i));
			}

			return u;
		}

		template<class TAcc>
		LMAT_ENSURE_INLINE
		static result_t reduce(const fun_t& fun, index_t len, const TAcc& tacc)
		{
			return fun.get(accumulate(fun, len, tacc), len);
		}


		template<class TAcc, class Out>
		LMAT_ENSURE_INLINE
		static void parallel_accum(const fun_t& fun, index_t len, const TAcc& tacc, Out& out)
		{
			for (index_t i = 0; i < len; ++i)
			{
				out[i] = fun.combine(out[i], tacc.get_scalar(i));
			}
		}

		template<class TAcc, class Out>
		LMAT_ENSURE_INLINE
		static void init_terms(const fun_t& fun, index_t len, const TAcc& tacc, Out& out)
		{
			for (index_t i = 0; i < len; ++i)
			{
				out[i] = tacc.get_scalar(i);
			}
		}
	};


} }

#endif
