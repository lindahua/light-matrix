/**
 * @file reduction_functors.h
 *
 * Basic reduction functors
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_REDUCTION_FUNCTORS_H_
#define LIGHTMAT_REDUCTION_FUNCTORS_H_

#include <light_mat/math/functor_base.h>
#include <light_mat/math/math_base.h>

namespace lmat
{
	template<typename T>
	struct sum_fun : public reduction_functor<T,T>
	{
		LMAT_ENSURE_INLINE
		T operator()() const
		{
			return T(0);
		}

		LMAT_ENSURE_INLINE
		T operator()(const T& x) const
		{
			return x;
		}

		LMAT_ENSURE_INLINE
		T operator()(const T& x, const T& y) const
		{
			return x + y;
		}
	};

	template<typename T>
	struct prod_fun : public reduction_functor<T,T>
	{
		LMAT_ENSURE_INLINE
		T operator()() const
		{
			return T(1);
		}

		LMAT_ENSURE_INLINE
		T operator()(const T& x) const
		{
			return x;
		}

		LMAT_ENSURE_INLINE
		T operator()(const T& x, const T& y) const
		{
			return x * y;
		}
	};


	template<typename T>
	struct maximum_fun : public reduction_functor<T,T>
	{
		LMAT_ENSURE_INLINE
		T operator()() const
		{
			throw invalid_operation("Taking maximum of an empty array is not allowed.");
		}

		LMAT_ENSURE_INLINE
		T operator()(const T& x) const
		{
			return x;
		}

		LMAT_ENSURE_INLINE
		T operator()(const T& x, const T& y) const
		{
			return (math::max)(x, y);
		}
	};

	template<typename T>
	struct minimum_fun : public reduction_functor<T,T>
	{
		LMAT_ENSURE_INLINE
		T operator()() const
		{
			throw invalid_operation("Taking minimum of an empty array is not allowed.");
		}

		LMAT_ENSURE_INLINE
		T operator()(const T& x) const
		{
			return x;
		}

		LMAT_ENSURE_INLINE
		T operator()(const T& x, const T& y) const
		{
			return (math::min)(x, y);
		}
	};


	// Declaration as reduction functors

	LMAT_DECLARE_AS_REDUCTION_TFUNCTOR( sum_fun, false )
	LMAT_DECLARE_AS_REDUCTION_TFUNCTOR( prod_fun, false )
	LMAT_DECLARE_AS_REDUCTION_TFUNCTOR( maximum_fun, false )
	LMAT_DECLARE_AS_REDUCTION_TFUNCTOR( minimum_fun, false )

}

#endif





