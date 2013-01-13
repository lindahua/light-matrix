/**
 * @file math_functors.h
 *
 * @brief Math functors
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATH_FUNCTORS_H_
#define LIGHTMAT_MATH_FUNCTORS_H_

#include <light_mat/math/basic_functors.h>
#include <light_mat/math/math.h>
#include <light_mat/math/simd_math.h>

namespace lmat
{
	// power functions

	_LMAT_DEFINE_REAL_MATH_FUN( pow, 2 )
	_LMAT_DEFINE_REAL_MATH_FUN( cbrt, 1 )
	_LMAT_DEFINE_REAL_MATH_FUN( hypot, 2 )

	// exp & log

	_LMAT_DEFINE_REAL_MATH_FUN( exp, 1 )
	_LMAT_DEFINE_REAL_MATH_FUN( log, 1 )
	_LMAT_DEFINE_REAL_MATH_FUN( log10, 1 )
	_LMAT_DEFINE_REAL_MATH_FUN( xlogx, 1 )
	_LMAT_DEFINE_REAL_MATH_FUN( xlogy, 2 )

	_LMAT_DEFINE_REAL_MATH_FUN( exp2, 1 )
	_LMAT_DEFINE_REAL_MATH_FUN( log2, 1 )
	_LMAT_DEFINE_REAL_MATH_FUN( expm1, 1 )
	_LMAT_DEFINE_REAL_MATH_FUN( log1p, 1 )

	// trigonometry

	_LMAT_DEFINE_REAL_MATH_FUN( sin, 1 )
	_LMAT_DEFINE_REAL_MATH_FUN( cos, 1 )
	_LMAT_DEFINE_REAL_MATH_FUN( tan, 1 )

	_LMAT_DEFINE_REAL_MATH_FUN( asin, 1 )
	_LMAT_DEFINE_REAL_MATH_FUN( acos, 1 )
	_LMAT_DEFINE_REAL_MATH_FUN( atan, 1 )
	_LMAT_DEFINE_REAL_MATH_FUN( atan2, 2 )

	// hyperbolic

	_LMAT_DEFINE_REAL_MATH_FUN( sinh, 1 )
	_LMAT_DEFINE_REAL_MATH_FUN( cosh, 1 )
	_LMAT_DEFINE_REAL_MATH_FUN( tanh, 1 )

	_LMAT_DEFINE_REAL_MATH_FUN( asinh, 1 )
	_LMAT_DEFINE_REAL_MATH_FUN( acosh, 1 )
	_LMAT_DEFINE_REAL_MATH_FUN( atanh, 1 )

}

#endif
