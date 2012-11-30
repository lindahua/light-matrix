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

	template<typename Tag, typename Ker, int Len, typename T1, typename T2=nil_t>
	struct vec_reduce_core;

	template<typename Tag, int Len, typename T>
	struct vec_reduce_core<Tag, scalar_kernel_t, Len, T, nil_t>
	{
		typedef typename reduction_result<Tag, T>::type result_t;
		typedef typename reduction_result<Tag, T>::intermediate_type media_t;
		typedef typename reduction_fun<Tag, scalar_kernel_t, T>::type fun_t;

		template<class Acc>
		LMAT_ENSURE_INLINE
		media_t accumulate(const fun_t& fun, index_t, const Acc& acc)
		{
			media_t u = fun.transform(acc.get_scalar(0));

			for (index_t i = 1; i < Len; ++i)
			{
				u = fun.combine(u, fun.transform(acc.get_scalar(i)));
			}

			return u;
		}

		template<class Acc>
		LMAT_ENSURE_INLINE
		result_t reduce(const fun_t& fun, index_t, const Acc& acc)
		{
			return fun.get(accumulate(fun, Len, acc), Len);
		}


		template<class Acc, class Out>
		LMAT_ENSURE_INLINE
		void parallel_accum(const fun_t& fun, index_t, const Acc& acc, Out& out)
		{
			for (index_t i = 0; i < Len; ++i)
			{
				out[i] = fun.combine(out[i], fun.transform(acc.get_scalar(i)));
			}
		}

		template<class Acc, class Out>
		LMAT_ENSURE_INLINE
		void init_terms(const fun_t& fun, index_t, const Acc& acc, Out& out)
		{
			for (index_t i = 0; i < Len; ++i)
			{
				out[i] = fun.transform(acc.get_scalar(i));
			}
		}
	};


	template<typename Tag, typename T>
	struct vec_reduce_core<Tag, scalar_kernel_t, 0, T, nil_t>
	{
		typedef typename reduction_result<Tag, T>::type result_t;
		typedef typename reduction_result<Tag, T>::intermediate_type media_t;
		typedef typename reduction_fun<Tag, scalar_kernel_t, T>::type fun_t;

		template<class Acc>
		LMAT_ENSURE_INLINE
		media_t accumulate(const fun_t& fun, index_t len, const Acc& acc)
		{
			media_t u = fun.transform(acc.get_scalar(0));

			for (index_t i = 1; i < len; ++i)
			{
				u = fun.combine(u, fun.transform(acc.get_scalar(i)));
			}

			return u;
		}

		template<class Acc>
		LMAT_ENSURE_INLINE
		result_t reduce(const fun_t& fun, index_t len, const Acc& acc)
		{
			return fun.get(accumulate(fun, len, acc), len);
		}

		template<class Acc, class Out>
		LMAT_ENSURE_INLINE
		void parallel_accum(const fun_t& fun, index_t len, const Acc& acc, Out& out)
		{
			for (index_t i = 0; i < len; ++i)
			{
				out[i] = fun.combine(out[i], fun.transform(acc.get_scalar(i)));
			}
		}

		template<class Acc, class Out>
		LMAT_ENSURE_INLINE
		void init_terms(const fun_t& fun, index_t len, const Acc& acc, Out& out)
		{
			for (index_t i = 0; i < len; ++i)
			{
				out[i] = fun.transform(acc.get_scalar(i));
			}
		}
	};


	template<typename Tag, int Len, typename T1, typename T2>
	struct vec_reduce_core<Tag, scalar_kernel_t, Len, T1, T2>
	{
		typedef typename reduction_result<Tag, T1, T2>::type result_t;
		typedef typename reduction_result<Tag, T1, T2>::intermediate_type media_t;
		typedef typename reduction_fun<Tag, scalar_kernel_t, T1, T2>::type fun_t;

		template<class Acc1, class Acc2>
		LMAT_ENSURE_INLINE
		media_t accumulate(const fun_t& fun, index_t, const Acc1& acc1, const Acc2& acc2)
		{
			media_t u = fun.transform(acc1.get_scalar(0), acc2.get_scalar(0));

			for (index_t i = 1; i < Len; ++i)
			{
				u = fun.combine(u, fun.transform(acc1.get_scalar(i), acc2.get_scalar(i)));
			}

			return u;
		}

		template<class Acc1, class Acc2>
		LMAT_ENSURE_INLINE
		result_t reduce(const fun_t& fun, index_t, const Acc1& acc1, const Acc2& acc2)
		{
			return fun.get(accumulate(fun, Len, acc1, acc2), Len);
		}


		template<class Acc1, class Acc2, class Out>
		LMAT_ENSURE_INLINE
		void parallel_accum(const fun_t& fun, index_t, const Acc1& acc1, const Acc2& acc2, Out& out)
		{
			for (index_t i = 0; i < Len; ++i)
			{
				out[i] = fun.combine(out[i], fun.transform(acc1.get_scalar(i), acc2.get_scalar(i)));
			}
		}

		template<class Acc1, class Acc2, class Out>
		LMAT_ENSURE_INLINE
		void init_terms(const fun_t& fun, index_t, const Acc1& acc1, const Acc2& acc2, Out& out)
		{
			for (index_t i = 0; i < Len; ++i)
			{
				out[i] = fun.transform(acc1.get_scalar(i), acc2.get_scalar(i));
			}
		}
	};


	template<typename Tag, typename T1, typename T2>
	struct vec_reduce_core<Tag, scalar_kernel_t, 0, T1, T2>
	{
		typedef typename reduction_result<Tag, T1, T2>::type result_t;
		typedef typename reduction_result<Tag, T1, T2>::intermediate_type media_t;
		typedef typename reduction_fun<Tag, scalar_kernel_t, T1, T2>::type fun_t;

		template<class Acc1, class Acc2>
		LMAT_ENSURE_INLINE
		media_t accumulate(const fun_t& fun, index_t len, const Acc1& acc1, const Acc2& acc2)
		{
			media_t u = fun.transform(acc1.get_scalar(0), acc2.get_scalar(0));

			for (index_t i = 1; i < len; ++i)
			{
				u = fun.combine(u, fun.transform(acc1.get_scalar(i), acc2.get_scalar(i)));
			}

			return u;
		}

		template<class Acc1, class Acc2>
		LMAT_ENSURE_INLINE
		result_t reduce(const fun_t& fun, index_t len, const Acc1& acc1, const Acc2& acc2)
		{
			return fun.get(accumulate(fun, len, acc1, acc2), len);
		}


		template<class Acc1, class Acc2, class Out>
		LMAT_ENSURE_INLINE
		void parallel_accum(const fun_t& fun, index_t len, const Acc1& acc1, const Acc2& acc2, Out& out)
		{
			for (index_t i = 0; i < len; ++i)
			{
				out[i] = fun.combine(out[i], fun.transform(acc1.get_scalar(i), acc2.get_scalar(i)));
			}
		}

		template<class Acc1, class Acc2, class Out>
		LMAT_ENSURE_INLINE
		void init_terms(const fun_t& fun, index_t len, const Acc1& acc1, const Acc2& acc2, Out& out)
		{
			for (index_t i = 0; i < len; ++i)
			{
				out[i] = fun.transform(acc1.get_scalar(i), acc2.get_scalar(i));
			}
		}
	};


} }

#endif
