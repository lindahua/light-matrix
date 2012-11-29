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
#include <light_mat/math/reduction_functors.h>

namespace lmat { namespace internal {

	template<typename RT, int CTLen, typename KerCate> struct vec_reduce;

	template<typename RT, int CTLen>
	struct vec_reduce<RT, CTLen, scalar_kernel_t>
	{
		template<class Fun, class Vec>
		LMAT_ENSURE_INLINE
		static RT eval(const Fun& fun, const index_t, const Vec& vec)
		{
			RT r = fun.transform(vec.get_scalar(0));
			for (index_t i = 1; i < CTLen; ++i)
			{
				r = fun.combine(r, fun.transform(vec.get_scalar(i)));
			}
			return r;
		}

		template<class Fun, class Vec1, class Vec2>
		LMAT_ENSURE_INLINE
		static RT eval(const Fun& fun, const index_t, const Vec1& vec1, const Vec2& vec2)
		{
			RT r = fun.transform(vec1.get_scalar(0), vec2.get_scalar(0));
			for (index_t i = 1; i < CTLen; ++i)
			{
				r = fun.combine(r, fun.transform(vec1.get_scalar(i), vec2.get_scalar(i)));
			}
			return r;
		}

		template<class Fun, class Vec, typename State>
		LMAT_ENSURE_INLINE
		static RT eval_s(const Fun& fun, const index_t, const Vec& vec, const State& s)
		{
			RT r = fun.transform(vec.get_scalar(0, s));
			for (index_t i = 1; i < CTLen; ++i)
			{
				r = fun.combine(r, fun.transform(vec.get_scalar(i, s)));
			}
			return r;
		}

		template<class Fun, class Vec1, typename State1, class Vec2, typename State2>
		LMAT_ENSURE_INLINE
		static RT eval_s(const Fun& fun, const index_t,
				const Vec1& vec1, const State1& s1, const Vec2& vec2, const State2& s2)
		{
			RT r = fun.transform(vec1.get_scalar(0, s1), vec2.get_scalar(0, s2));
			for (index_t i = 1; i < CTLen; ++i)
			{
				r = fun.combine(r, fun.transform(vec1.get_scalar(i, s1), vec2.get_scalar(i, s2)));
			}
			return r;
		}
	};


	template<typename RT>
	struct vec_reduce<RT, 0, scalar_kernel_t>
	{
		template<class Fun, class Vec>
		LMAT_ENSURE_INLINE
		static RT eval(const Fun& fun, const index_t len, const Vec& vec)
		{
			RT r = fun.transform(vec.get_scalar(0));
			for (index_t i = 1; i < len; ++i)
			{
				r = fun.combine(r, fun.transform(vec.get_scalar(i)));
			}
			return r;
		}

		template<class Fun, class Vec1, class Vec2>
		LMAT_ENSURE_INLINE
		static RT eval(const Fun& fun, const index_t len, const Vec1& vec1, const Vec2& vec2)
		{
			RT r = fun.transform(vec1.get_scalar(0), vec2.get_scalar(0));
			for (index_t i = 1; i < len; ++i)
			{
				r = fun.combine(r, fun.transform(vec1.get_scalar(i), vec2.get_scalar(i)));
			}
			return r;
		}

		template<class Fun, class Vec, typename State>
		LMAT_ENSURE_INLINE
		static RT eval_s(const Fun& fun, const index_t len, const Vec& vec, const State& s)
		{
			RT r = fun.transform(vec.get_scalar(0, s));
			for (index_t i = 1; i < len; ++i)
			{
				r = fun.combine(r, fun.transform(vec.get_scalar(i, s)));
			}
			return r;
		}

		template<class Fun, class Vec1, typename State1, class Vec2, typename State2>
		LMAT_ENSURE_INLINE
		static RT eval_s(const Fun& fun, const index_t len,
				const Vec1& vec1, const State1& s1, const Vec2& vec2, const State2& s2)
		{
			RT r = fun.transform(vec1.get_scalar(0, s1), vec2.get_scalar(0, s2));
			for (index_t i = 1; i < len; ++i)
			{
				r = fun.combine(r, fun.transform(vec1.get_scalar(i, s1), vec2.get_scalar(i, s2)));
			}
			return r;
		}
	};


	template<typename RT>
	struct vec_reduce<RT, 1, scalar_kernel_t>
	{
		template<class Fun, class Vec>
		LMAT_ENSURE_INLINE
		static RT eval(const Fun& fun, const index_t len, const Vec& vec)
		{
			return fun.transform(vec.get_scalar(0));
		}

		template<class Fun, class Vec1, class Vec2>
		LMAT_ENSURE_INLINE
		static RT eval(const Fun& fun, const index_t len, const Vec1& vec1, const Vec2& vec2)
		{
			return fun.transform(vec1.get_scalar(0), vec2.get_scalar(0));
		}

		template<class Fun, class Vec, typename State>
		LMAT_ENSURE_INLINE
		static RT eval_s(const Fun& fun, const index_t len, const Vec& vec, const State& s)
		{
			return fun.transform(vec.get_scalar(0, s));
		}

		template<class Fun, class Vec1, typename State1, class Vec2, typename State2>
		LMAT_ENSURE_INLINE
		static RT eval_s(const Fun& fun, const index_t len,
				const Vec1& vec1, const State1& s1, const Vec2& vec2, const State2& s2)
		{
			return fun.transform(vec1.get_scalar(0, s1), vec2.get_scalar(0, s2));
		}
	};



} }

#endif
