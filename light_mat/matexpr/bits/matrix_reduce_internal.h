/**
 * @file matrix_reduce_internal.h
 *
 * Internal implementation of matrix reduction
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_REDUCE_INTERNAL_H_
#define LIGHTMAT_MATRIX_REDUCE_INTERNAL_H_

#include <light_mat/matrix/matrix_properties.h>
#include <light_mat/matrix/matrix_visitors.h>

namespace lmat { namespace detail {

	template<int CTLen, typename KerCate> struct single_vec_reduce;


	template<int CTLen>
	struct single_vec_reduce<CTLen, scalar_kernel_t>
	{
		template<class Fun, class Vec>
		LMAT_ENSURE_INLINE
		static typename Fun::result_type eval(const Fun& fun, const index_t, const Vec& vec)
		{
			typename Fun::result_type r = fun(vec.get_scalar(0));
			for (index_t i = 1; i < CTLen; ++i)
			{
				r = fun(r, vec.get_scalar(i));
			}
			return r;
		}

		template<class Fun, class Vec, typename State>
		LMAT_ENSURE_INLINE
		static typename Fun::result_type eval(const Fun& fun, const index_t, const Vec& vec, const State& s)
		{
			typename Fun::result_type r = fun(vec.get_scalar(0, s));
			for (index_t i = 1; i < CTLen; ++i)
			{
				r = fun(r, vec.get_scalar(i, s));
			}
			return r;
		}
	};

	template<>
	struct single_vec_reduce<0, scalar_kernel_t>
	{
		template<class Fun, class Vec>
		LMAT_ENSURE_INLINE
		static typename Fun::result_type eval(const Fun& fun, const index_t len, const Vec& vec)
		{
			if (len > 0)
			{
				typename Fun::result_type r = fun(vec.get_scalar(0));
				for (index_t i = 1; i < len; ++i)
				{
					r = fun(r, vec.get_scalar(i));
				}
				return r;
			}
			else
			{
				return fun();
			}
		}

		template<class Fun, class Vec, typename State>
		LMAT_ENSURE_INLINE
		static typename Fun::result_type eval(const Fun& fun, const index_t len, const Vec& vec, const State& s)
		{
			if (len > 0)
			{
				typename Fun::result_type r = fun(vec.get_scalar(0, s));
				for (index_t i = 1; i < len; ++i)
				{
					r = fun(r, vec.get_scalar(i, s));
				}
				return r;
			}
			else
			{
				return fun();
			}
		}
	};

	template<>
	struct single_vec_reduce<1, scalar_kernel_t>
	{
		template<class Fun, class Vec>
		LMAT_ENSURE_INLINE
		static typename Fun::result_type eval(const Fun& fun, const index_t, const Vec& vec)
		{
			return fun(vec.get_scalar(0));
		}

		template<class Fun, class Vec, typename State>
		LMAT_ENSURE_INLINE
		static typename Fun::result_type eval(const Fun& fun, const index_t, const Vec& vec, const State& s)
		{
			return fun(vec.get_scalar(0, s));
		}
	};


} }

#endif
