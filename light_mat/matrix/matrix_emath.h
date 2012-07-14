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
	 *  Expression Type mapping
	 *
	 ********************************************/

	LMAT_DECLARE_UNARY_TYPE_MAP( floor, floor_fun )
	LMAT_DECLARE_UNARY_TYPE_MAP( ceil,  ceil_fun )

	LMAT_DECLARE_BINARY_TYPE_MAP_EX( pow, pow_fun )

	LMAT_DECLARE_UNARY_TYPE_MAP( exp, exp_fun )
	LMAT_DECLARE_UNARY_TYPE_MAP( log, log_fun )
	LMAT_DECLARE_UNARY_TYPE_MAP( log10, log10_fun )

	LMAT_DECLARE_UNARY_TYPE_MAP( sin, sin_fun )
	LMAT_DECLARE_UNARY_TYPE_MAP( cos, cos_fun )
	LMAT_DECLARE_UNARY_TYPE_MAP( tan, tan_fun )

	LMAT_DECLARE_UNARY_TYPE_MAP( asin, asin_fun )
	LMAT_DECLARE_UNARY_TYPE_MAP( acos, acos_fun )
	LMAT_DECLARE_UNARY_TYPE_MAP( atan, atan_fun )
	LMAT_DECLARE_BINARY_TYPE_MAP_EX( atan2, atan2_fun )

	LMAT_DECLARE_UNARY_TYPE_MAP( sinh, sinh_fun )
	LMAT_DECLARE_UNARY_TYPE_MAP( cosh, cosh_fun )
	LMAT_DECLARE_UNARY_TYPE_MAP( tanh, tanh_fun )

#ifdef LMAT_HAS_CXX11_MATH

	LMAT_DECLARE_UNARY_TYPE_MAP( round, round_fun )
	LMAT_DECLARE_UNARY_TYPE_MAP( trunc, trunc_fun )

	LMAT_DECLARE_BINARY_TYPE_MAP( hypot, hypot_fun )
	LMAT_DECLARE_UNARY_TYPE_MAP( cbrt, cbrt_fun )

	LMAT_DECLARE_UNARY_TYPE_MAP( exp2, exp2_fun )
	LMAT_DECLARE_UNARY_TYPE_MAP( log2, log2_fun )
	LMAT_DECLARE_UNARY_TYPE_MAP( expm1, expm1_fun )
	LMAT_DECLARE_UNARY_TYPE_MAP( log1p, log1p_fun )

	LMAT_DECLARE_UNARY_TYPE_MAP( asinh, asinh_fun )
	LMAT_DECLARE_UNARY_TYPE_MAP( acosh, acosh_fun )
	LMAT_DECLARE_UNARY_TYPE_MAP( atanh, atanh_fun )

#endif


	/********************************************
	 *
	 *  Specific Math Expressions
	 *
	 ********************************************/

	LMAT_DEFINE_UNARY_MATFUNCTION( floor, floor, floor_fun )
	LMAT_DEFINE_UNARY_MATFUNCTION( ceil,  ceil,  ceil_fun )

	LMAT_DEFINE_BINARY_MATFUNCTION_EX( pow, pow, pow_fun )

	LMAT_DEFINE_UNARY_MATFUNCTION( exp, exp, exp_fun )
	LMAT_DEFINE_UNARY_MATFUNCTION( log, log, log_fun )
	LMAT_DEFINE_UNARY_MATFUNCTION( log10, log10, log10_fun )

	LMAT_DEFINE_UNARY_MATFUNCTION( sin, sin, sin_fun )
	LMAT_DEFINE_UNARY_MATFUNCTION( cos, cos, cos_fun )
	LMAT_DEFINE_UNARY_MATFUNCTION( tan, tan, tan_fun )

	LMAT_DEFINE_UNARY_MATFUNCTION( asin, asin, asin_fun )
	LMAT_DEFINE_UNARY_MATFUNCTION( acos, acos, acos_fun )
	LMAT_DEFINE_UNARY_MATFUNCTION( atan, atan, atan_fun )
	LMAT_DEFINE_BINARY_MATFUNCTION( atan2, atan2, atan2_fun )

	LMAT_DEFINE_UNARY_MATFUNCTION( sinh, sinh, sinh_fun )
	LMAT_DEFINE_UNARY_MATFUNCTION( cosh, cosh, cosh_fun )
	LMAT_DEFINE_UNARY_MATFUNCTION( tanh, tanh, tanh_fun )

#ifdef LMAT_HAS_CXX11_MATH

	LMAT_DEFINE_UNARY_MATFUNCTION( round, round, round_fun )
	LMAT_DEFINE_UNARY_MATFUNCTION( trunc, trunc, trunc_fun )

	LMAT_DEFINE_BINARY_MATFUNCTION( hypot, hypot, hypot_fun )
	LMAT_DEFINE_UNARY_MATFUNCTION( cbrt, cbrt, cbrt_fun )

	LMAT_DEFINE_UNARY_MATFUNCTION( exp2, exp2, exp2_fun )
	LMAT_DEFINE_UNARY_MATFUNCTION( log2, log2, log2_fun )
	LMAT_DEFINE_UNARY_MATFUNCTION( expm1, expm1, expm1_fun )
	LMAT_DEFINE_UNARY_MATFUNCTION( log1p, log1p, log1p_fun )

	LMAT_DEFINE_UNARY_MATFUNCTION( asinh, asinh, asinh_fun )
	LMAT_DEFINE_UNARY_MATFUNCTION( acosh, acosh, acosh_fun )
	LMAT_DEFINE_UNARY_MATFUNCTION( atanh, atanh, atanh_fun )

#endif

}

#endif /* MATRIX_EMATH_H_ */


