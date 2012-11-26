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

#include <light_mat/matexpr/matrix_ewise_eval.h>
#include <light_mat/math/emath_functors.h>

namespace lmat
{

	/********************************************
	 *
	 *  Specific Math Expressions
	 *
	 ********************************************/

	LMAT_DEFINE_UNARY_MATFUNCTION( sqr,   sqr_t )
	LMAT_DEFINE_UNARY_MATFUNCTION( cube,  cube_t )
	LMAT_DEFINE_UNARY_MATFUNCTION( sqrt,  sqrt_t )
	LMAT_DEFINE_UNARY_MATFUNCTION( rcp,   rcp_t )
	LMAT_DEFINE_UNARY_MATFUNCTION( rsqrt, rsqrt_t )

	LMAT_DEFINE_BINARY_MATFUNCTION( pow, pow_t )
	LMAT_DEFINE_BINARY_MATFUNCTION( mod, mod_t )

	LMAT_DEFINE_UNARY_MATFUNCTION( floor, floor_t )
	LMAT_DEFINE_UNARY_MATFUNCTION( ceil,  ceil_t )

	LMAT_DEFINE_UNARY_MATFUNCTION( exp,   exp_t )
	LMAT_DEFINE_UNARY_MATFUNCTION( log,   log_t )
	LMAT_DEFINE_UNARY_MATFUNCTION( log10, log10_t )

	LMAT_DEFINE_UNARY_MATFUNCTION( sin, sin_t )
	LMAT_DEFINE_UNARY_MATFUNCTION( cos, cos_t )
	LMAT_DEFINE_UNARY_MATFUNCTION( tan, tan_t )

	LMAT_DEFINE_UNARY_MATFUNCTION( asin, asin_t )
	LMAT_DEFINE_UNARY_MATFUNCTION( acos, acos_t )
	LMAT_DEFINE_UNARY_MATFUNCTION( atan, atan_t )
	LMAT_DEFINE_BINARY_MATFUNCTION( atan2, atan2_t )

	LMAT_DEFINE_UNARY_MATFUNCTION( sinh, sinh_t )
	LMAT_DEFINE_UNARY_MATFUNCTION( cosh, cosh_t )
	LMAT_DEFINE_UNARY_MATFUNCTION( tanh, tanh_t )

#ifdef LMAT_HAS_C99_MATH

	LMAT_DEFINE_UNARY_MATFUNCTION( cbrt,  cbrt_t )
	LMAT_DEFINE_UNARY_MATFUNCTION( hypot, hypot_t )

	LMAT_DEFINE_UNARY_MATFUNCTION( round, round_t )
	LMAT_DEFINE_UNARY_MATFUNCTION( trunc, trunc_t )

	LMAT_DEFINE_UNARY_MATFUNCTION( exp2, exp2_t )
	LMAT_DEFINE_UNARY_MATFUNCTION( log2, log2_t )
	LMAT_DEFINE_UNARY_MATFUNCTION( expm1, expm1_t )
	LMAT_DEFINE_UNARY_MATFUNCTION( log1p, log1p_t )

	LMAT_DEFINE_UNARY_MATFUNCTION( asinh, asinh_t )
	LMAT_DEFINE_UNARY_MATFUNCTION( acosh, acosh_t )
	LMAT_DEFINE_UNARY_MATFUNCTION( atanh, atanh_t )

	LMAT_DEFINE_UNARY_MATFUNCTION( erf,    erf_t )
	LMAT_DEFINE_UNARY_MATFUNCTION( erfc,   erfc_t )
	LMAT_DEFINE_UNARY_MATFUNCTION( lgamma, lgamma_t )
	LMAT_DEFINE_UNARY_MATFUNCTION( tgamma, tgamma_t )

#endif

}

#endif /* MATRIX_EMATH_H_ */


