/*
 * @file math_base.h
 *
 * Basic definitions for math computation
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATH_BASE_H_
#define LIGHTMAT_MATH_BASE_H_

#include <light_mat/math/fun_tags.h>

#if ((LIGHTMAT_PLATFORM == LIGHTMAT_POSIX) || defined(__INTEL_COMPILER))
#define LMAT_PLATFORM_HAS_CXX11_MATH
#define LMAT_HAS_CXX11_MATH
#endif

#include <cstdlib>
#include <cmath>

namespace lmat { namespace math {

	// arithmetics

	using std::fma;

	template<typename T>
	LMAT_ENSURE_INLINE inline T (max)(const T& x, const T& y) { return x > y ? x : y; }

	template<typename T>
	LMAT_ENSURE_INLINE inline T (min)(const T& x, const T& y) { return x < y ? x : y; }

	LMAT_ENSURE_INLINE inline float  (max)(float  x, float  y) { return std::fmax(x, y); }
	LMAT_ENSURE_INLINE inline double (max)(double x, double y) { return std::fmax(x, y); }
	LMAT_ENSURE_INLINE inline float  (min)(float  x, float  y) { return std::fmin(x, y); }
	LMAT_ENSURE_INLINE inline double (min)(double x, double y) { return std::fmin(x, y); }

	template<typename T>
	LMAT_ENSURE_INLINE inline T clamp(const T& x, const T& lb, const T& ub)
	{
		return (min)((max)(x, lb), ub);
	}

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


	// simple power functions

	using std::abs;
	using std::sqrt;

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

	// rounding

	using std::floor;
	using std::ceil;
	using std::round;
	using std::trunc;

	using std::lrint;
	using std::lround;

	// classification

	using std::signbit;
	using std::isinf;
	using std::isfinite;
	using std::isnan;


} }

#endif 
