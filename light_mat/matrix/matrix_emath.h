/**
 * @file matrix_emath.h
 *
 * Element-wise math functions on matrix
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_EMATH_H_
#define LIGHTMAT_MATRIX_EMATH_H_

#include <light_mat/matrix/matrix_ewise_expr.h>
#include <light_mat/math/emath_functors.h>

namespace lmat
{

	/********************************************
	 *
	 *  Specific Math Expressions
	 *
	 ********************************************/

	LMAT_DEFINE_UNARY_MATFUNCTION( floor, floor_fun )
	LMAT_DEFINE_UNARY_MATFUNCTION( ceil,  ceil_fun )

	LMAT_DEFINE_BINARY_MATFUNCTION_EX( pow, pow_fun )

	LMAT_DEFINE_UNARY_MATFUNCTION( exp, exp_fun )
	LMAT_DEFINE_UNARY_MATFUNCTION( log, log_fun )
	LMAT_DEFINE_UNARY_MATFUNCTION( log10, log10_fun )

	LMAT_DEFINE_UNARY_MATFUNCTION( sin, sin_fun )
	LMAT_DEFINE_UNARY_MATFUNCTION( cos, cos_fun )
	LMAT_DEFINE_UNARY_MATFUNCTION( tan, tan_fun )

	LMAT_DEFINE_UNARY_MATFUNCTION( asin, asin_fun )
	LMAT_DEFINE_UNARY_MATFUNCTION( acos, acos_fun )
	LMAT_DEFINE_UNARY_MATFUNCTION( atan, atan_fun )
	LMAT_DEFINE_BINARY_MATFUNCTION( atan2, atan2_fun )

	LMAT_DEFINE_UNARY_MATFUNCTION( sinh, sinh_fun )
	LMAT_DEFINE_UNARY_MATFUNCTION( cosh, cosh_fun )
	LMAT_DEFINE_UNARY_MATFUNCTION( tanh, tanh_fun )

#ifdef LMAT_HAS_CXX11_MATH

	LMAT_DEFINE_UNARY_MATFUNCTION( round, round_fun )
	LMAT_DEFINE_UNARY_MATFUNCTION( trunc, trunc_fun )

	LMAT_DEFINE_BINARY_MATFUNCTION( hypot, hypot_fun )
	LMAT_DEFINE_UNARY_MATFUNCTION( cbrt, cbrt_fun )

	LMAT_DEFINE_UNARY_MATFUNCTION( exp2, exp2_fun )
	LMAT_DEFINE_UNARY_MATFUNCTION( log2, log2_fun )
	LMAT_DEFINE_UNARY_MATFUNCTION( expm1, expm1_fun )
	LMAT_DEFINE_UNARY_MATFUNCTION( log1p, log1p_fun )

	LMAT_DEFINE_UNARY_MATFUNCTION( asinh, asinh_fun )
	LMAT_DEFINE_UNARY_MATFUNCTION( acosh, acosh_fun )
	LMAT_DEFINE_UNARY_MATFUNCTION( atanh, atanh_fun )

#endif

}

#endif /* MATRIX_EMATH_H_ */


