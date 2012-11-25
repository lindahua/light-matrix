/*
 * @file math_base.h
 *
 * Basic definitions for math computation
 *
 * @author Dahua Lin
 */

#ifndef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATH_BASE_H_
#define LIGHTMAT_MATH_BASE_H_

#include <light_mat/common/basic_defs.h>
#include <cstdlib>
#include <math.h>

#define LMAT_IMPORT_STDMATH_UFUN( funname, fun ) \
	LMAT_ENSURE_INLINE inline float funname(float x) { return ::fun##f(x); } \
	LMAT_ENSURE_INLINE inline double funname(double x) { return ::fun(x); }

#define LMAT_IMPORT_STDMATH_BFUN( funname, fun ) \
	LMAT_ENSURE_INLINE inline float funname(float x, float y) { return ::fun##f(x, y); } \
	LMAT_ENSURE_INLINE inline double funname(double x, double y) { return ::fun(x, y); }

namespace lmat { namespace math {

	/********************************************
	 *
	 *  standard math functions
	 *
	 ********************************************/

	// abs

	LMAT_ENSURE_INLINE inline int abs(int x) { return std::abs(x); }
	LMAT_ENSURE_INLINE inline long abs(long x) { return std::abs(x); }
	LMAT_IMPORT_STDMATH_UFUN( abs, fabs )

	// floor & ceil

	LMAT_IMPORT_STDMATH_UFUN( floor, floor )
	LMAT_IMPORT_STDMATH_UFUN( ceil, ceil )

	// power & module

	LMAT_IMPORT_STDMATH_UFUN( sqrt, sqrt )
	LMAT_IMPORT_STDMATH_BFUN( pow, pow )
	LMAT_IMPORT_STDMATH_BFUN( mod, fmod )

	LMAT_ENSURE_INLINE inline int    sqr(int    x) { return x * x; }
	LMAT_ENSURE_INLINE inline long   sqr(long   x) { return x * x; }
	LMAT_ENSURE_INLINE inline float  sqr(float  x) { return x * x; }
	LMAT_ENSURE_INLINE inline double sqr(double x) { return x * x; }

	LMAT_ENSURE_INLINE inline int    cube(int    x) { return x * x * x; }
	LMAT_ENSURE_INLINE inline long 	 cube(long   x) { return x * x * x; }
	LMAT_ENSURE_INLINE inline float  cube(float  x) { return x * x * x; }
	LMAT_ENSURE_INLINE inline double cube(double x) { return x * x * x; }

	LMAT_ENSURE_INLINE inline float  rcp(float  x) { return 1.0f / x; }
	LMAT_ENSURE_INLINE inline double rcp(double x) { return 1.0 / x; }

	LMAT_ENSURE_INLINE inline float  rsqrt(float  x) { return 1.0f / ::sqrtf(x); }
	LMAT_ENSURE_INLINE inline double rsqrt(double x) { return 1.0 / ::sqrt(x); }

	// exp & log

	LMAT_IMPORT_STDMATH_UFUN( exp, exp )
	LMAT_IMPORT_STDMATH_UFUN( log, log )
	LMAT_IMPORT_STDMATH_UFUN( log10, log10 )

	// trigonometry

	LMAT_IMPORT_STDMATH_UFUN( sin, sin )
	LMAT_IMPORT_STDMATH_UFUN( cos, cos )
	LMAT_IMPORT_STDMATH_UFUN( tan, atan )

	LMAT_IMPORT_STDMATH_UFUN( asin, asin )
	LMAT_IMPORT_STDMATH_UFUN( acos, acos )
	LMAT_IMPORT_STDMATH_UFUN( atan, atan )

	LMAT_IMPORT_STDMATH_BFUN( atan2, atan2 )

	// hyperbolic

	LMAT_IMPORT_STDMATH_UFUN( sinh, sinh )
	LMAT_IMPORT_STDMATH_UFUN( cosh, cosh )
	LMAT_IMPORT_STDMATH_UFUN( tanh, atanh )

	// max & min

	template<typename T>
	LMAT_ENSURE_INLINE T (max)(const T& x, const T& y) { return x > y ? x : y; }

	template<typename T>
	LMAT_ENSURE_INLINE T (min)(const T& x, const T& y) { return x < y ? x : y; }

#ifdef LMAT_HAS_C99_MATH
	LMAT_ENSURE_INLINE inline float  (max)(float  x, float  y) { return ::fmaxf(x, y); }
	LMAT_ENSURE_INLINE inline double (max)(double x, double y) { return ::fmax(x, y); }
	LMAT_ENSURE_INLINE inline float  (min)(float  x, float  y) { return ::fminf(x, y); }
	LMAT_ENSURE_INLINE inline double (min)(double x, double y) { return ::fmin(x, y); }
#endif


	/********************************************
	 *
	 *  extended math functions (in C99/C++11)
	 *
	 ********************************************/

#ifdef LMAT_HAS_C99_MATH

	// fma

	LMAT_ENSURE_INLINE inline float fma(float x, float y, float z) { return ::fmaf(x, y, z); }
	LMAT_ENSURE_INLINE inline double fma(double x, double y, double z) { return ::fma(x, y, z); }

	// rounding

	LMAT_IMPORT_STDMATH_UFUN( round, round )
	LMAT_IMPORT_STDMATH_UFUN( trunc, trunc )

	LMAT_ENSURE_INLINE long lround( float x ) { return ::lrintf(x); }
	LMAT_ENSURE_INLINE long lround( double x ) { return ::lrint(x); }

	// power

	LMAT_IMPORT_STDMATH_UFUN( cbrt, cbrt )
	LMAT_IMPORT_STDMATH_BFUN( hypot, hypot )

	// exp & log

	LMAT_IMPORT_STDMATH_UFUN( exp2, exp2 )
	LMAT_IMPORT_STDMATH_UFUN( log2, log2 )
	LMAT_IMPORT_STDMATH_UFUN( expm1, expm1 )
	LMAT_IMPORT_STDMATH_UFUN( log1p, log1p )

	// hyperbolic

	LMAT_IMPORT_STDMATH_UFUN( asinh, asinh )
	LMAT_IMPORT_STDMATH_UFUN( acosh, acosh )
	LMAT_IMPORT_STDMATH_UFUN( atanh, atanh )

	// error & gamma

	LMAT_IMPORT_STDMATH_UFUN( erf, erf )
	LMAT_IMPORT_STDMATH_UFUN( erfc, erfc )

	LMAT_IMPORT_STDMATH_UFUN( lgamma, lgamma )
	LMAT_IMPORT_STDMATH_UFUN( tgamma, tgamma )

	// classification

	LMAT_ENSURE_INLINE bool is_inf(float x) { return (bool)isinf(x); }
	LMAT_ENSURE_INLINE bool is_inf(double x) { return (bool)isinf(x); }

	LMAT_ENSURE_INLINE bool is_finite(float x) { return (bool)isfinite(x); }
	LMAT_ENSURE_INLINE bool is_finite(double x) { return (bool)isfinite(x); }

	LMAT_ENSURE_INLINE bool is_nan(float x) { return (bool)isnan(x); }
	LMAT_ENSURE_INLINE bool is_nan(double x) { return (bool)isnan(x); }

#endif


} }

#endif 
