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

#define LMAT_DEFINE_UNARY_ARITH_FUNCTOR( fname, expr ) \
	template<typename T> struct fname##_fun { \
		LMAT_ENSURE_INLINE fname##_fun() { } \
		LMAT_ENSURE_INLINE fname##_fun( fname##_t ) { } \
		LMAT_ENSURE_INLINE T operator() (const T& x) const { return expr; } \
	}; \
	LMAT_DEFINE_NUMERIC_UNARY_FUNMAP( fname##_t, scalar_kernel_t, fname##_fun )

#define LMAT_DEFINE_BINARY_ARITH_FUNCTOR( fname, expr ) \
	template<typename T> struct fname##_fun { \
		LMAT_ENSURE_INLINE fname##_fun() { } \
		LMAT_ENSURE_INLINE fname##_fun( fname##_t ) { } \
		LMAT_ENSURE_INLINE T operator() (const T& x, const T& y) const { return expr; } \
	}; \
	LMAT_DEFINE_NUMERIC_BINARY_FUNMAP( fname##_t, scalar_kernel_t, fname##_fun )


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

	LMAT_DEFINE_BINARY_ARITH_FUNCTOR( add, x + y )
	LMAT_DEFINE_BINARY_ARITH_FUNCTOR( subtract, x - y )
	LMAT_DEFINE_BINARY_ARITH_FUNCTOR( multiply, x * y )
	LMAT_DEFINE_BINARY_ARITH_FUNCTOR( divide, x / y )

	LMAT_DEFINE_UNARY_ARITH_FUNCTOR( negate, - x )
	LMAT_DEFINE_UNARY_ARITH_FUNCTOR( abs, math::abs(x) )

	LMAT_DEFINE_BINARY_ARITH_FUNCTOR( min, (math::min)(x, y) )
	LMAT_DEFINE_BINARY_ARITH_FUNCTOR( max, (math::max)(x, y) )
}


#endif 





