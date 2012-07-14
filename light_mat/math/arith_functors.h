/*
 * @file arith_functors.h
 *
 * Arithmetic functors
 *
 * @author Dahua Lin
 */

#ifndef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_ARITH_FUNCTORS_H_
#define LIGHTMAT_ARITH_FUNCTORS_H_

#include <light_mat/math/functor_base.h>
#include <light_mat/math/math_base.h>

namespace lmat
{
	template<typename T>
	struct add_op : public binary_numeric_ewise_functor<T>
	{
		LMAT_ENSURE_INLINE T operator() (const T& x, const T& y) const
		{
			return x + y;
		}
	};

	template<typename T>
	struct sub_op : public binary_numeric_ewise_functor<T>
	{
		LMAT_ENSURE_INLINE T operator() (const T& x, const T& y) const
		{
			return x - y;
		}
	};

	template<typename T>
	struct mul_op : public binary_numeric_ewise_functor<T>
	{
		LMAT_ENSURE_INLINE T operator() (const T& x, const T& y) const
		{
			return x * y;
		}
	};

	template<typename T>
	struct div_op : public binary_numeric_ewise_functor<T>
	{
		LMAT_ENSURE_INLINE T operator() (const T& x, const T& y) const
		{
			return x / y;
		}
	};

	template<typename T>
	struct neg_op : public unary_numeric_ewise_functor<T>
	{
		LMAT_ENSURE_INLINE T operator() (const T& x) const
		{
			return - x;
		}
	};

	template<typename T>
	struct abs_fun : public unary_numeric_ewise_functor<T>
	{
		LMAT_ENSURE_INLINE T operator() (const T& x) const
		{
			return math::abs(x);
		}
	};

	template<typename T>
	struct sqr_fun : public unary_numeric_ewise_functor<T>
	{
		LMAT_ENSURE_INLINE T operator() (const T& x) const
		{
			return math::sqr(x);
		}
	};

	template<typename T>
	struct rcp_fun : public unary_numeric_ewise_functor<T>
	{
		LMAT_ENSURE_INLINE T operator() (const T& x) const
		{
			return math::rcp(x);
		}
	};

	template<typename T>
	struct sqrt_fun : public unary_numeric_ewise_functor<T>
	{
		LMAT_ENSURE_INLINE T operator() (const T& x) const
		{
			return math::sqrt(x);
		}
	};

	template<typename T>
	struct rsqrt_fun : public unary_numeric_ewise_functor<T>
	{
		LMAT_ENSURE_INLINE T operator() (const T& x) const
		{
			return math::rsqrt(x);
		}
	};

	// declaration as ewise functors

	LMAT_DECLARE_AS_BINARY_EWISE_TFUNCTOR( add_op, true )
	LMAT_DECLARE_AS_BINARY_EWISE_TFUNCTOR( sub_op, true )
	LMAT_DECLARE_AS_BINARY_EWISE_TFUNCTOR( mul_op, true )
	LMAT_DECLARE_AS_BINARY_EWISE_TFUNCTOR( div_op, true )

	LMAT_DECLARE_AS_UNARY_EWISE_TFUNCTOR( neg_op, true )
	LMAT_DECLARE_AS_UNARY_EWISE_TFUNCTOR( abs_fun, true )
	LMAT_DECLARE_AS_UNARY_EWISE_TFUNCTOR( sqr_fun, true )
	LMAT_DECLARE_AS_UNARY_EWISE_TFUNCTOR( sqrt_fun, true )
	LMAT_DECLARE_AS_UNARY_EWISE_TFUNCTOR( rcp_fun, true )
	LMAT_DECLARE_AS_UNARY_EWISE_TFUNCTOR( rsqrt_fun, true )

	LMAT_DEFINE_BINARY_NUMERIC_EWISE_TFUNCTOR( min_fun, (math::min), true )
	LMAT_DEFINE_BINARY_NUMERIC_EWISE_TFUNCTOR( max_fun, (math::max), true )
}


#endif 





