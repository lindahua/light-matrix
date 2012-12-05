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

	LMAT_DECLARE_NUMERIC_UNARY_OP( id_t )

	LMAT_DECLARE_NUMERIC_BINARY_OP( add_t )
	LMAT_DECLARE_NUMERIC_BINARY_OP( subtract_t )
	LMAT_DECLARE_NUMERIC_BINARY_OP( multiply_t )
	LMAT_DECLARE_NUMERIC_BINARY_OP( divide_t )

	LMAT_DECLARE_NUMERIC_UNARY_OP( negate_t )
	LMAT_DECLARE_NUMERIC_UNARY_OP( abs_t )

	LMAT_DECLARE_NUMERIC_BINARY_OP( max_t )
	LMAT_DECLARE_NUMERIC_BINARY_OP( min_t )

	// define functors

	LMAT_DEFINE_NUMERIC_UNARY_FUNCTOR( id, x )

	LMAT_DEFINE_NUMERIC_BINARY_FUNCTOR( add, x + y )
	LMAT_DEFINE_NUMERIC_BINARY_FUNCTOR( subtract, x - y )
	LMAT_DEFINE_NUMERIC_BINARY_FUNCTOR( multiply, x * y )
	LMAT_DEFINE_NUMERIC_BINARY_FUNCTOR( divide, x / y )

	LMAT_DEFINE_NUMERIC_UNARY_FUNCTOR( negate, - x )
	LMAT_DEFINE_NUMERIC_UNARY_FUNCTOR( abs, math::abs(x) )

	LMAT_DEFINE_NUMERIC_BINARY_FUNCTOR( min, (math::min)(x, y) )
	LMAT_DEFINE_NUMERIC_BINARY_FUNCTOR( max, (math::max)(x, y) )
}


#endif 





