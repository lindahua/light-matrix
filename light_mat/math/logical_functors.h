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

namespace lmat
{
	LMAT_DECLARE_LOGICAL_UNARY_OP( not_t )
	LMAT_DECLARE_LOGICAL_BINARY_OP( and_t )
	LMAT_DECLARE_LOGICAL_BINARY_OP( or_t )
	LMAT_DECLARE_LOGICAL_BINARY_OP( xor_t )

	LMAT_DEFINE_LOGICAL_UNARY_FUNCTOR( not, ~x )
	LMAT_DEFINE_LOGICAL_BINARY_FUNCTOR( and, x & y )
	LMAT_DEFINE_LOGICAL_BINARY_FUNCTOR( or, x | y )
	LMAT_DEFINE_LOGICAL_BINARY_FUNCTOR( xor, x ^ y )
}

#endif 
