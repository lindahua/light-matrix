/**
 * @file math_functors.h
 *
 * Math functors
 * 
 * @author Dahua Lin 
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATH_FUNCTORS_H_
#define LIGHTMAT_MATH_FUNCTORS_H_

#include <light_mat/math/math_base.h>


#define LMAT_DEFINE_GENERIC_MATH_FUNCTOR_1( Name, Expr ) \
	struct Name##_fun { \
		template<typename T> \
		LMAT_ENSURE_INLINE \
		T operator()(const T& x) { return Expr; } \
	};

#define LMAT_DEFINE_GENERIC_MATH_FUNCTOR_2( Name, Expr ) \
	struct Name##_fun { \
		template<typename T> \
		LMAT_ENSURE_INLINE \
		T operator()(const T& x, const T& y) { return Expr; } \
	};

#define LMAT_DEFINE_MATH_FUNCTOR_1( Name ) \
	LMAT_DEFINE_GENERIC_MATH_FUNCTOR_1( Name, Name(x) )

#define LMAT_DEFINE_MATH_FUNCTOR_2( Name ) \
	LMAT_DEFINE_GENERIC_MATH_FUNCTOR_2( Name, Name(x, y) )

namespace lmat { namespace math {

	// Arithmetic functors

	LMAT_DEFINE_GENERIC_MATH_FUNCTOR_2( add, x + y )
	LMAT_DEFINE_GENERIC_MATH_FUNCTOR_2( sub, x - y )
	LMAT_DEFINE_GENERIC_MATH_FUNCTOR_2( mul, x * y )
	LMAT_DEFINE_GENERIC_MATH_FUNCTOR_2( div, x / y )
	LMAT_DEFINE_GENERIC_MATH_FUNCTOR_1( negate, -x )

	LMAT_DEFINE_MATH_FUNCTOR_1( abs )
	LMAT_DEFINE_MATH_FUNCTOR_2( max )
	LMAT_DEFINE_MATH_FUNCTOR_2( min )

	// power

	LMAT_DEFINE_MATH_FUNCTOR_1( sqr )
	LMAT_DEFINE_MATH_FUNCTOR_1( cube )
	LMAT_DEFINE_MATH_FUNCTOR_1( rcp )
	LMAT_DEFINE_MATH_FUNCTOR_1( sqrt )
	LMAT_DEFINE_MATH_FUNCTOR_1( rsqrt )
	LMAT_DEFINE_MATH_FUNCTOR_2( pow )

	// floor & ceil

	LMAT_DEFINE_MATH_FUNCTOR_1( floor )
	LMAT_DEFINE_MATH_FUNCTOR_1( ceil )

	// exp & log

	LMAT_DEFINE_MATH_FUNCTOR_1( exp )
	LMAT_DEFINE_MATH_FUNCTOR_1( log )
	LMAT_DEFINE_MATH_FUNCTOR_1( log10 )

	// trigonometry

	LMAT_DEFINE_MATH_FUNCTOR_1( sin )
	LMAT_DEFINE_MATH_FUNCTOR_1( cos )
	LMAT_DEFINE_MATH_FUNCTOR_1( tan )

	LMAT_DEFINE_MATH_FUNCTOR_1( asin )
	LMAT_DEFINE_MATH_FUNCTOR_1( acos )
	LMAT_DEFINE_MATH_FUNCTOR_1( atan )
	LMAT_DEFINE_MATH_FUNCTOR_2( atan2 )

	// hyperbolic

	LMAT_DEFINE_MATH_FUNCTOR_1( sinh )
	LMAT_DEFINE_MATH_FUNCTOR_1( cosh )
	LMAT_DEFINE_MATH_FUNCTOR_1( tanh )

#ifdef LMAT_HAS_CXX11_MATH

	// cbrt and hypot

	LMAT_DEFINE_MATH_FUNCTOR_1( cbrt )
	LMAT_DEFINE_MATH_FUNCTOR_2( hypot )

	// rounding

	LMAT_DEFINE_MATH_FUNCTOR_1( round )
	LMAT_DEFINE_MATH_FUNCTOR_1( trunc )

	// exp & log

	LMAT_DEFINE_MATH_FUNCTOR_1( exp2 )
	LMAT_DEFINE_MATH_FUNCTOR_1( log2 )
	LMAT_DEFINE_MATH_FUNCTOR_1( expm1 )
	LMAT_DEFINE_MATH_FUNCTOR_1( log1p )

	// inverse hyperbolic

	LMAT_DEFINE_MATH_FUNCTOR_1( asinh )
	LMAT_DEFINE_MATH_FUNCTOR_1( acosh )
	LMAT_DEFINE_MATH_FUNCTOR_1( atanh )

	// error and gamma

	LMAT_DEFINE_MATH_FUNCTOR_1( erf )
	LMAT_DEFINE_MATH_FUNCTOR_1( erfc )
	LMAT_DEFINE_MATH_FUNCTOR_1( lgamma )
	LMAT_DEFINE_MATH_FUNCTOR_1( tgamma )

#endif


} }

#endif 


