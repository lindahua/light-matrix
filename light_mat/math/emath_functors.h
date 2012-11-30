/**
 * @file emath_functors.h
 *
 * Elementary math functors
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_EMATH_FUNCTORS_H_
#define LIGHTMAT_EMATH_FUNCTORS_H_

#include <light_mat/math/functor_base.h>
#include <light_mat/math/math_base.h>


namespace lmat
{
	/********************************************
	 *
	 *  Declarations
	 *
	 ********************************************/

	// power & module

	LMAT_DECLARE_REAL_UNARY_OP( sqr_t )
	LMAT_DECLARE_REAL_UNARY_OP( cube_t )
	LMAT_DECLARE_REAL_UNARY_OP( sqrt_t )
	LMAT_DECLARE_REAL_UNARY_OP( cbrt_t )
	LMAT_DECLARE_REAL_UNARY_OP( rcp_t )
	LMAT_DECLARE_REAL_UNARY_OP( rsqrt_t )

	LMAT_DECLARE_REAL_BINARY_OP( pow_t )
	LMAT_DECLARE_REAL_BINARY_OP( fmod_t )
	LMAT_DECLARE_REAL_BINARY_OP( hypot_t )

	// rounding

	LMAT_DECLARE_REAL_UNARY_OP( floor_t )
	LMAT_DECLARE_REAL_UNARY_OP( ceil_t )
	LMAT_DECLARE_REAL_UNARY_OP( round_t )
	LMAT_DECLARE_REAL_UNARY_OP( trunc_t )

	// exponential and logarithm

	LMAT_DECLARE_REAL_UNARY_OP( exp_t )
	LMAT_DECLARE_REAL_UNARY_OP( exp2_t )
	LMAT_DECLARE_REAL_UNARY_OP( log_t )
	LMAT_DECLARE_REAL_UNARY_OP( log2_t )
	LMAT_DECLARE_REAL_UNARY_OP( log10_t )

	LMAT_DECLARE_REAL_UNARY_OP( expm1_t )
	LMAT_DECLARE_REAL_UNARY_OP( log1p_t )

	// trigonometry

	LMAT_DECLARE_REAL_UNARY_OP( sin_t )
	LMAT_DECLARE_REAL_UNARY_OP( cos_t )
	LMAT_DECLARE_REAL_UNARY_OP( tan_t )

	LMAT_DECLARE_REAL_UNARY_OP( asin_t )
	LMAT_DECLARE_REAL_UNARY_OP( acos_t )
	LMAT_DECLARE_REAL_UNARY_OP( atan_t )

	LMAT_DECLARE_REAL_BINARY_OP( atan2_t )

	// hyperbolic

	LMAT_DECLARE_REAL_UNARY_OP( sinh_t )
	LMAT_DECLARE_REAL_UNARY_OP( cosh_t )
	LMAT_DECLARE_REAL_UNARY_OP( tanh_t )

	LMAT_DECLARE_REAL_UNARY_OP( asinh_t )
	LMAT_DECLARE_REAL_UNARY_OP( acosh_t )
	LMAT_DECLARE_REAL_UNARY_OP( atanh_t )

	// error & gamma

	LMAT_DECLARE_REAL_UNARY_OP( erf_t )
	LMAT_DECLARE_REAL_UNARY_OP( erfc_t )
	LMAT_DECLARE_REAL_UNARY_OP( lgamma_t )
	LMAT_DECLARE_REAL_UNARY_OP( tgamma_t )

	/********************************************
	 *
	 *  elementary math operations
	 *
	 ********************************************/

	// power & module

	LMAT_DEFINE_REAL_UNARY_FUNCTOR( sqr, math::sqr(x) )
	LMAT_DEFINE_REAL_UNARY_FUNCTOR( cube, math::cube(x) )

	LMAT_DEFINE_REAL_UNARY_FUNCTOR( sqrt, math::sqrt(x) )
	LMAT_DEFINE_REAL_UNARY_FUNCTOR( cbrt, math::cbrt(x) )

	LMAT_DEFINE_REAL_UNARY_FUNCTOR( rcp, math::rcp(x) )
	LMAT_DEFINE_REAL_UNARY_FUNCTOR( rsqrt, math::rsqrt(x) )

	LMAT_DEFINE_REAL_BINARY_FUNCTOR( pow, math::pow(x, y) )
	LMAT_DEFINE_REAL_BINARY_FUNCTOR( fmod, math::fmod(x, y) )

#ifdef LMAT_HAS_C99_MATH
	LMAT_DEFINE_REAL_BINARY_FUNCTOR( hypot, math::hypot(x, y) )
#endif

	// rounding

	LMAT_DEFINE_REAL_UNARY_FUNCTOR( floor, math::floor(x) )
	LMAT_DEFINE_REAL_UNARY_FUNCTOR( ceil, math::ceil(x) )

#ifdef LMAT_HAS_C99_MATH
	LMAT_DEFINE_REAL_UNARY_FUNCTOR( round, math::round(x) )
	LMAT_DEFINE_REAL_UNARY_FUNCTOR( trunc, math::trunc(x) )
#endif

	// exponential and logarithm

	LMAT_DEFINE_REAL_UNARY_FUNCTOR( exp, math::exp(x) )
	LMAT_DEFINE_REAL_UNARY_FUNCTOR( log, math::log(x) )
	LMAT_DEFINE_REAL_UNARY_FUNCTOR( log10, math::log10(x) )

#ifdef LMAT_HAS_C99_MATH
	LMAT_DEFINE_REAL_UNARY_FUNCTOR( exp2, math::exp2(x) )
	LMAT_DEFINE_REAL_UNARY_FUNCTOR( log2, math::log2(x) )
	LMAT_DEFINE_REAL_UNARY_FUNCTOR( expm1, math::expm1(x) )
	LMAT_DEFINE_REAL_UNARY_FUNCTOR( log1p, math::log1p(x) )
#endif

	// trigonometry

	LMAT_DEFINE_REAL_UNARY_FUNCTOR( sin, math::sin(x) )
	LMAT_DEFINE_REAL_UNARY_FUNCTOR( cos, math::cos(x) )
	LMAT_DEFINE_REAL_UNARY_FUNCTOR( tan, math::tan(x) )

	LMAT_DEFINE_REAL_UNARY_FUNCTOR( asin, math::asin(x) )
	LMAT_DEFINE_REAL_UNARY_FUNCTOR( acos, math::acos(x) )
	LMAT_DEFINE_REAL_UNARY_FUNCTOR( atan, math::atan(x) )
	LMAT_DEFINE_REAL_BINARY_FUNCTOR( atan2, math::atan2(x, y) )

	// hyperbolic

	LMAT_DEFINE_REAL_UNARY_FUNCTOR( sinh, math::sinh(x) )
	LMAT_DEFINE_REAL_UNARY_FUNCTOR( cosh, math::cosh(x) )
	LMAT_DEFINE_REAL_UNARY_FUNCTOR( tanh, math::tanh(x) )

#ifdef LMAT_HAS_C99_MATH
	LMAT_DEFINE_REAL_UNARY_FUNCTOR( asinh, math::asinh(x) )
	LMAT_DEFINE_REAL_UNARY_FUNCTOR( acosh, math::acosh(x) )
	LMAT_DEFINE_REAL_UNARY_FUNCTOR( atanh, math::atanh(x) )
#endif LMAT_HAS_C99_MATH

	// error & gamma

#ifdef LMAT_HAS_C99_MATH
	LMAT_DEFINE_REAL_UNARY_FUNCTOR( erf, math::erf(x) )
	LMAT_DEFINE_REAL_UNARY_FUNCTOR( erfc, math::erfc(x) )
	LMAT_DEFINE_REAL_UNARY_FUNCTOR( lgamma, math::lgamma(x) )
	LMAT_DEFINE_REAL_UNARY_FUNCTOR( tgamma, math::tgamma(x) )
#endif

}

#endif /* EMATH_FUNCTORS_H_ */



