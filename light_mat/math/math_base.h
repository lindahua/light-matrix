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
#include <cmath>

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

	using std::abs;

	// floor & ceil

	using std::floor;
	using std::ceil;

	// power & module

	using std::sqrt;
	using std::pow;
	using std::fmod;

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

	LMAT_ENSURE_INLINE inline float  rsqrt(float  x) { return 1.0f / sqrt(x); }
	LMAT_ENSURE_INLINE inline double rsqrt(double x) { return 1.0 / sqrt(x); }

	// exp & log

	using std::exp;
	using std::log;
	using std::log10;

	// trigonometry

	using std::sin;
	using std::cos;
	using std::tan;

	using std::asin;
	using std::acos;
	using std::atan;
	using std::atan2;

	// hyperbolic

	using std::sinh;
	using std::cosh;
	using std::tanh;

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

	LMAT_ENSURE_INLINE inline long lround( float x ) { return ::lrintf(x); }
	LMAT_ENSURE_INLINE inline long lround( double x ) { return ::lrint(x); }

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

	LMAT_ENSURE_INLINE inline bool is_inf(float x) { return (bool)isinf(x); }
	LMAT_ENSURE_INLINE inline bool is_inf(double x) { return (bool)isinf(x); }

	// LMAT_ENSURE_INLINE bool is_finite(float x) { return (bool)isfinite(x); }
	// LMAT_ENSURE_INLINE bool is_finite(double x) { return (bool)isfinite(x); }

	LMAT_ENSURE_INLINE inline bool is_nan(float x) { return (bool)isnan(x); }
	LMAT_ENSURE_INLINE inline bool is_nan(double x) { return (bool)isnan(x); }

#endif


} }

#endif 
