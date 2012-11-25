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
	// arithmetic operations

	LMAT_DEFINE_NUMERIC_BINARY_OP( add_t )
	LMAT_DEFINE_NUMERIC_BINARY_OP( subtract_t )
	LMAT_DEFINE_NUMERIC_BINARY_OP( multiply_t )
	LMAT_DEFINE_NUMERIC_BINARY_OP( divide_t )

	LMAT_DEFINE_NUMERIC_UNARY_OP( negate_t )
	LMAT_DEFINE_NUMERIC_UNARY_OP( abs_t )

	LMAT_DEFINE_NUMERIC_BINARY_OP( max_t )
	LMAT_DEFINE_NUMERIC_BINARY_OP( min_t )

	// define functors

	template<typename T>
	struct add_fun
	{
		LMAT_ENSURE_INLINE add_fun() { }

		LMAT_ENSURE_INLINE add_fun( add_t ) { }

		LMAT_ENSURE_INLINE T operator() (const T& x, const T& y) const
		{
			return x + y;
		}
	};

	template<typename T>
	struct subtract_fun
	{
		LMAT_ENSURE_INLINE subtract_fun() { }

		LMAT_ENSURE_INLINE subtract_fun( subtract_t ) { }

		LMAT_ENSURE_INLINE T operator() (const T& x, const T& y) const
		{
			return x - y;
		}
	};

	template<typename T>
	struct multiply_fun
	{
		LMAT_ENSURE_INLINE multiply_fun() { }

		LMAT_ENSURE_INLINE multiply_fun( multiply_t ) { }

		LMAT_ENSURE_INLINE T operator() (const T& x, const T& y) const
		{
			return x * y;
		}
	};

	template<typename T>
	struct divide_fun
	{
		LMAT_ENSURE_INLINE divide_fun() { }

		LMAT_ENSURE_INLINE divide_fun( divide_t ) { }

		LMAT_ENSURE_INLINE T operator() (const T& x, const T& y) const
		{
			return x / y;
		}
	};

	template<typename T>
	struct negate_fun
	{
		LMAT_ENSURE_INLINE negate_fun() { }

		LMAT_ENSURE_INLINE negate_fun( negate_t ) { }

		LMAT_ENSURE_INLINE T operator() (const T& x) const
		{
			return -x;
		}
	};

	template<typename T>
	struct abs_fun
	{
		LMAT_ENSURE_INLINE abs_fun() { }

		LMAT_ENSURE_INLINE abs_fun( abs_t ) { }

		LMAT_ENSURE_INLINE T operator() (const T& x) const
		{
			return math::abs(x);
		}
	};

	template<typename T>
	struct max_fun
	{
		LMAT_ENSURE_INLINE max_fun() { }

		LMAT_ENSURE_INLINE max_fun( max_t ) { }

		LMAT_ENSURE_INLINE T operator() (const T& x, const T& y) const
		{
			return (math::max)(x, y);
		}
	};

	template<typename T>
	struct min_fun
	{
		LMAT_ENSURE_INLINE min_fun() { }

		LMAT_ENSURE_INLINE min_fun( min_t ) { }

		LMAT_ENSURE_INLINE T operator() (const T& x, const T& y) const
		{
			return (math::min)(x, y);
		}
	};


	// functor maps

	LMAT_DEFINE_NUMERIC_BINARY_FUNMAP( add_t,      scalar_kernel_t, add_fun )
	LMAT_DEFINE_NUMERIC_BINARY_FUNMAP( subtract_t, scalar_kernel_t, subtract_fun )
	LMAT_DEFINE_NUMERIC_BINARY_FUNMAP( multiply_t, scalar_kernel_t, multiply_fun )
	LMAT_DEFINE_NUMERIC_BINARY_FUNMAP( divide_t,   scalar_kernel_t, divide_fun )

	LMAT_DEFINE_NUMERIC_UNARY_FUNMAP( negate_t, scalar_kernel_t, negate_fun )
	LMAT_DEFINE_NUMERIC_UNARY_FUNMAP( abs_t,    scalar_kernel_t, abs_fun )

	LMAT_DEFINE_NUMERIC_BINARY_FUNMAP( max_t, scalar_kernel_t, max_fun )
	LMAT_DEFINE_NUMERIC_BINARY_FUNMAP( min_t, scalar_kernel_t, min_fun )
}


#endif 





