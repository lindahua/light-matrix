/*
 * @file logical_functors.h
 *
 * Logical functors for not, and, or, xor
 *
 * @author Dahua Lin
 */

#ifndef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_LOGICAL_FUNCTORS_H_
#define LIGHTMAT_LOGICAL_FUNCTORS_H_

#include <light_mat/math/functor_base.h>

#define LMAT_DEFINE_UNARY_LOGICAL_FUNCTOR( fname, expr ) \
	struct fname##_fun { \
		LMAT_ENSURE_INLINE fname##_fun() { } \
		LMAT_ENSURE_INLINE fname##_fun( fname##_t ) { } \
		LMAT_ENSURE_INLINE bool operator() (const bool& x) const { return expr; } \
	}; \
	LMAT_DEFINE_NUMERIC_UNARY_FUNMAP( fname##_t, scalar_kernel_t, fname##_fun )

#define LMAT_DEFINE_BINARY_LOGICAL_FUNCTOR( fname, expr ) \
	struct fname##_fun { \
		LMAT_ENSURE_INLINE fname##_fun() { } \
		LMAT_ENSURE_INLINE fname##_fun( fname##_t ) { } \
		LMAT_ENSURE_INLINE bool operator() (const bool& x, const bool& y) const { return expr; } \
	}; \
	LMAT_DEFINE_NUMERIC_BINARY_FUNMAP( fname##_t, scalar_kernel_t, fname##_fun )


namespace lmat
{
	LMAT_DEFINE_LOGICAL_UNARY_OP( not_t )
	LMAT_DEFINE_LOGICAL_BINARY_OP( and_t )
	LMAT_DEFINE_LOGICAL_BINARY_OP( or_t )
	LMAT_DEFINE_LOGICAL_BINARY_OP( xor_t )

	LMAT_DEFINE_UNARY_LOGICAL_FUNCTOR( not, ~x )
	LMAT_DEFINE_BINARY_LOGICAL_FUNCTOR( and, x & y )
	LMAT_DEFINE_BINARY_LOGICAL_FUNCTOR( or, x | y )
	LMAT_DEFINE_BINARY_LOGICAL_FUNCTOR( xor, x ^ y )
}

#endif 
