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


#define LMAT_DEFINE_UNARY_MATH_FUNCTOR( fname ) \
	template<typename T> struct fname##_fun { \
		LMAT_ENSURE_INLINE fname##_fun() { } \
		LMAT_ENSURE_INLINE fname##_fun( fname##_t ) { } \
		LMAT_ENSURE_INLINE T operator() (const T& x) const { return math::fname(x); } \
	}; \
	LMAT_DEFINE_REAL_UNARY_FUNMAP( fname##_t, scalar_kernel_t, fname##_fun )

#define LMAT_DEFINE_BINARY_MATH_FUNCTOR( fname ) \
	template<typename T> struct fname##_fun { \
		LMAT_ENSURE_INLINE fname##_fun() { } \
		LMAT_ENSURE_INLINE fname##_fun( fname##_t ) { } \
		LMAT_ENSURE_INLINE T operator() (const T& x, const T& y) const { return math::fname(x, y); } \
	}; \
	LMAT_DEFINE_REAL_BINARY_FUNMAP( fname##_t, scalar_kernel_t, fname##_fun )


namespace lmat
{
	/********************************************
	 *
	 *  elementary math operations
	 *
	 ********************************************/

	// power & module

	LMAT_DEFINE_NUMERIC_UNARY_OP( sqr_t )
	LMAT_DEFINE_NUMERIC_UNARY_OP( cube_t )

	LMAT_DEFINE_REAL_UNARY_OP( sqrt_t )
	LMAT_DEFINE_REAL_UNARY_OP( cbrt_t )

	LMAT_DEFINE_REAL_UNARY_OP( rcp_t )
	LMAT_DEFINE_REAL_UNARY_OP( rsqrt_t )

	LMAT_DEFINE_REAL_BINARY_OP( pow_t )
	LMAT_DEFINE_REAL_BINARY_OP( fmod_t )
	LMAT_DEFINE_REAL_BINARY_OP( hypot_t )

	// rounding

	LMAT_DEFINE_REAL_UNARY_OP( floor_t )
	LMAT_DEFINE_REAL_UNARY_OP( ceil_t )
	LMAT_DEFINE_REAL_UNARY_OP( round_t )
	LMAT_DEFINE_REAL_UNARY_OP( trunc_t )

	// exponential and logarithm

	LMAT_DEFINE_REAL_UNARY_OP( exp_t )
	LMAT_DEFINE_REAL_UNARY_OP( exp2_t )
	LMAT_DEFINE_REAL_UNARY_OP( log_t )
	LMAT_DEFINE_REAL_UNARY_OP( log2_t )
	LMAT_DEFINE_REAL_UNARY_OP( log10_t )

	LMAT_DEFINE_REAL_UNARY_OP( expm1_t )
	LMAT_DEFINE_REAL_UNARY_OP( log1p_t )

	// trigonometry

	LMAT_DEFINE_REAL_UNARY_OP( sin_t )
	LMAT_DEFINE_REAL_UNARY_OP( cos_t )
	LMAT_DEFINE_REAL_UNARY_OP( tan_t )

	LMAT_DEFINE_REAL_UNARY_OP( asin_t )
	LMAT_DEFINE_REAL_UNARY_OP( acos_t )
	LMAT_DEFINE_REAL_UNARY_OP( atan_t )

	LMAT_DEFINE_REAL_BINARY_OP( atan2_t )

	// hyperbolic

	LMAT_DEFINE_REAL_UNARY_OP( sinh_t )
	LMAT_DEFINE_REAL_UNARY_OP( cosh_t )
	LMAT_DEFINE_REAL_UNARY_OP( tanh_t )

	LMAT_DEFINE_REAL_UNARY_OP( asinh_t )
	LMAT_DEFINE_REAL_UNARY_OP( acosh_t )
	LMAT_DEFINE_REAL_UNARY_OP( atanh_t )

	// error & gamma

	LMAT_DEFINE_REAL_UNARY_OP( erf_t )
	LMAT_DEFINE_REAL_UNARY_OP( erfc_t )
	LMAT_DEFINE_REAL_UNARY_OP( lgamma_t )
	LMAT_DEFINE_REAL_UNARY_OP( tgamma_t )


	/********************************************
	 *
	 *  functor classes
	 *
	 ********************************************/

	LMAT_DEFINE_UNARY_MATH_FUNCTOR( sqr )
	LMAT_DEFINE_UNARY_MATH_FUNCTOR( cube )
	LMAT_DEFINE_UNARY_MATH_FUNCTOR( sqrt )
	LMAT_DEFINE_UNARY_MATH_FUNCTOR( rcp )
	LMAT_DEFINE_UNARY_MATH_FUNCTOR( rsqrt )

	LMAT_DEFINE_BINARY_MATH_FUNCTOR( pow )
	LMAT_DEFINE_BINARY_MATH_FUNCTOR( fmod )

	LMAT_DEFINE_UNARY_MATH_FUNCTOR( floor )
	LMAT_DEFINE_UNARY_MATH_FUNCTOR( ceil )

	LMAT_DEFINE_UNARY_MATH_FUNCTOR( exp )
	LMAT_DEFINE_UNARY_MATH_FUNCTOR( log )
	LMAT_DEFINE_UNARY_MATH_FUNCTOR( log10 )

	LMAT_DEFINE_UNARY_MATH_FUNCTOR( sin )
	LMAT_DEFINE_UNARY_MATH_FUNCTOR( cos )
	LMAT_DEFINE_UNARY_MATH_FUNCTOR( tan )

	LMAT_DEFINE_UNARY_MATH_FUNCTOR( asin )
	LMAT_DEFINE_UNARY_MATH_FUNCTOR( acos )
	LMAT_DEFINE_UNARY_MATH_FUNCTOR( atan )

	LMAT_DEFINE_BINARY_MATH_FUNCTOR( atan2 )

	LMAT_DEFINE_UNARY_MATH_FUNCTOR( sinh )
	LMAT_DEFINE_UNARY_MATH_FUNCTOR( cosh )
	LMAT_DEFINE_UNARY_MATH_FUNCTOR( tanh )


#ifdef LMAT_HAS_C99_MATH

	LMAT_DEFINE_UNARY_MATH_FUNCTOR( cbrt )
	LMAT_DEFINE_BINARY_MATH_FUNCTOR( hypot )

	LMAT_DEFINE_UNARY_MATH_FUNCTOR( round )
	LMAT_DEFINE_UNARY_MATH_FUNCTOR( trunc )

	LMAT_DEFINE_UNARY_MATH_FUNCTOR( exp2 )
	LMAT_DEFINE_UNARY_MATH_FUNCTOR( log2 )
	LMAT_DEFINE_UNARY_MATH_FUNCTOR( expm1 )
	LMAT_DEFINE_UNARY_MATH_FUNCTOR( log1p )

	LMAT_DEFINE_UNARY_MATH_FUNCTOR( asinh )
	LMAT_DEFINE_UNARY_MATH_FUNCTOR( acosh )
	LMAT_DEFINE_UNARY_MATH_FUNCTOR( atanh )

	LMAT_DEFINE_UNARY_MATH_FUNCTOR( erf )
	LMAT_DEFINE_UNARY_MATH_FUNCTOR( erfc )
	LMAT_DEFINE_UNARY_MATH_FUNCTOR( lgamma )
	LMAT_DEFINE_UNARY_MATH_FUNCTOR( tgamma )

#endif

}

#endif /* EMATH_FUNCTORS_H_ */



