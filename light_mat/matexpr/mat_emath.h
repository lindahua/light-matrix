/**
 * @file mat_emath.h
 *
 * @brief Elemental math functions on matrices
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MAT_EMATH_H_
#define LIGHTMAT_MAT_EMATH_H_

#include <light_mat/matexpr/matfun_base.h>
#include <light_mat/math/math_functors.h>

namespace lmat
{

	// power functions

	_LMAT_DEFINE_RMATFUN( pow, 2 )
	_LMAT_DEFINE_RMATFUN( cbrt, 1 )
	_LMAT_DEFINE_RMATFUN( hypot, 2 )

	// exp & log

	_LMAT_DEFINE_RMATFUN( exp, 1 )
	_LMAT_DEFINE_RMATFUN( log, 1 )
	_LMAT_DEFINE_RMATFUN( log10, 1 )

	_LMAT_DEFINE_RMATFUN( xlogx, 1 )
	_LMAT_DEFINE_RMATFUN( xlogy, 2 )

	_LMAT_DEFINE_RMATFUN( exp2, 1 )
	_LMAT_DEFINE_RMATFUN( log2, 1 )
	_LMAT_DEFINE_RMATFUN( expm1, 1 )
	_LMAT_DEFINE_RMATFUN( log1p, 1 )

	// trigonometry

	_LMAT_DEFINE_RMATFUN( sin, 1 )
	_LMAT_DEFINE_RMATFUN( cos, 1 )
	_LMAT_DEFINE_RMATFUN( tan, 1 )

	_LMAT_DEFINE_RMATFUN( asin, 1 )
	_LMAT_DEFINE_RMATFUN( acos, 1 )
	_LMAT_DEFINE_RMATFUN( atan, 1 )
	_LMAT_DEFINE_RMATFUN( atan2, 2 )

	// hyperbolic

	_LMAT_DEFINE_RMATFUN( sinh, 1 )
	_LMAT_DEFINE_RMATFUN( cosh, 1 )
	_LMAT_DEFINE_RMATFUN( tanh, 1 )

	_LMAT_DEFINE_RMATFUN( asinh, 1 )
	_LMAT_DEFINE_RMATFUN( acosh, 1 )
	_LMAT_DEFINE_RMATFUN( atanh, 1 )

}

#endif
