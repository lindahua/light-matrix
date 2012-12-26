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
#include <light_mat/math/fun_tags.h>

#if ((LIGHTMAT_PLATFORM == LIGHTMAT_POSIX) || defined(__INTEL_COMPILER))
#define LMAT_PLATFORM_HAS_CXX11_MATH
#define LMAT_HAS_CXX11_MATH
#endif

#include <cstdlib>
#include <cmath>

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

	LMAT_ENSURE_INLINE inline float  diff_abs(float  x, float  y) { return abs(x - y); }
	LMAT_ENSURE_INLINE inline double diff_abs(double x, double y) { return abs(x - y); }

	LMAT_ENSURE_INLINE inline float  diff_sqr(float  x, float  y) { return sqr(x - y); }
	LMAT_ENSURE_INLINE inline double diff_sqr(double x, double y) { return sqr(x - y); }

	// conditional

	template<typename T>
	LMAT_ENSURE_INLINE
	inline T cond(bool tf, const T& x, const T& y)
	{
		return tf ? x : y;
	}

	template<typename T>
	LMAT_ENSURE_INLINE
	inline T cond(const mask_t<T>& m, const T& x, const T& y)
	{
		return m.bvalue ? x : y;
	}


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
	LMAT_ENSURE_INLINE inline T (max)(const T& x, const T& y) { return x > y ? x : y; }

	template<typename T>
	LMAT_ENSURE_INLINE inline T (min)(const T& x, const T& y) { return x < y ? x : y; }

#ifdef LMAT_HAS_CXX11_MATH
	LMAT_ENSURE_INLINE inline float  (max)(float  x, float  y) { return std::fmax(x, y); }
	LMAT_ENSURE_INLINE inline double (max)(double x, double y) { return std::fmax(x, y); }
	LMAT_ENSURE_INLINE inline float  (min)(float  x, float  y) { return std::fmin(x, y); }
	LMAT_ENSURE_INLINE inline double (min)(double x, double y) { return std::fmin(x, y); }
#endif

	template<typename T>
	LMAT_ENSURE_INLINE inline T clamp(const T& x, const T& lb, const T& ub)
	{
		return (min)((max)(x, lb), ub);
	}

	// xlogy

	LMAT_ENSURE_INLINE inline float xlogy(float x, float y)
	{
		return x > 0 ? x * log(y) : 0.f;
	}

	LMAT_ENSURE_INLINE inline double xlogy(double x, double y)
	{
		return x > 0 ? x * log(y) : 0.0;
	}


	/********************************************
	 *
	 *  extended math functions (in C99/C++11)
	 *
	 ********************************************/

#ifdef LMAT_HAS_CXX11_MATH

	// fma

	using std::fma;

	// rounding

	using std::round;
	using std::trunc;

	using std::lrint;
	using std::lround;

	// power

	using std::cbrt;
	using std::hypot;

	// exp & log

	using std::exp2;
	using std::log2;
	using std::expm1;
	using std::log1p;

	// hyperbolic

	using std::asinh;
	using std::acosh;
	using std::atanh;

	// error & gamma

	using std::erf;
	using std::erfc;

	using std::lgamma;
	using std::tgamma;

	// classification

	using std::signbit;
	using std::isinf;
	using std::isfinite;
	using std::isnan;

#else

	LMAT_ENSURE_INLINE inline float fma( float x, float y, float z ) { return x * y + z; }
	LMAT_ENSURE_INLINE inline double fma( double x, double y, double z ) { return x * y + z; }

#endif


} }

#endif 
