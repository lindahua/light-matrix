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
#include <light_mat/matrix/matrix_vector_eval.h>

namespace lmat { namespace detail {

	template<int CTLen, typename Means> struct single_vec_reduce;


	template<int CTLen>
	struct single_vec_reduce<CTLen, by_scalars>
	{
		template<class Fun, class Vec>
		LMAT_ENSURE_INLINE
		static typename Fun::result_type eval(const Fun& fun, const index_t, const Vec& vec)
		{
			typename Fun::result_type r = fun(vec.get_value(0));
			for (index_t i = 1; i < CTLen; ++i)
			{
				r = fun(r, vec.get_value(i));
			}
			return r;
		}

		template<class Fun, class Vec, typename State>
		LMAT_ENSURE_INLINE
		static typename Fun::result_type eval(const Fun& fun, const index_t, const Vec& vec, const State& s)
		{
			typename Fun::result_type r = fun(vec.get_value(s, 0));
			for (index_t i = 1; i < CTLen; ++i)
			{
				r = fun(r, vec.get_value(s, i));
			}
			return r;
		}
	};

	template<>
	struct single_vec_reduce<0, by_scalars>
	{
		template<class Fun, class Vec>
		LMAT_ENSURE_INLINE
		static typename Fun::result_type eval(const Fun& fun, const index_t len, const Vec& vec)
		{
			if (len > 0)
			{
				typename Fun::result_type r = fun(vec.get_value(0));
				for (index_t i = 1; i < len; ++i)
				{
					r = fun(r, vec.get_value(i));
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
				typename Fun::result_type r = fun(vec.get_value(s, 0));
				for (index_t i = 1; i < len; ++i)
				{
					r = fun(r, vec.get_value(s, i));
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
	struct single_vec_reduce<1, by_scalars>
	{
		template<class Fun, class Vec>
		LMAT_ENSURE_INLINE
		static typename Fun::result_type eval(const Fun& fun, const index_t, const Vec& vec)
		{
			return fun(vec.get_value(0));
		}

		template<class Fun, class Vec, typename State>
		LMAT_ENSURE_INLINE
		static typename Fun::result_type eval(const Fun& fun, const index_t, const Vec& vec, const State& s)
		{
			return fun(vec.get_value(s, 0));
		}
	};


} }

#endif
