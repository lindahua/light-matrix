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

#include <light_mat/core/basic_defs.h>
#include <cstdlib>
#include <cmath>

namespace lmat { namespace math {

	// import of C++03 standard math functions

	using std::abs;
	using std::sqrt;
	using std::pow;
	using std::modf;

	using std::floor;
	using std::ceil;

	using std::sin;
	using std::cos;
	using std::tan;

	using std::asin;
	using std::acos;
	using std::atan;
	using std::atan2;

	using std::sinh;
	using std::cosh;
	using std::tanh;

	using std::exp;
	using std::log;
	using std::log10;

	// additional functions

	LMAT_ENSURE_INLINE inline int32_t  sqr(int32_t  x) { return x * x; }
	LMAT_ENSURE_INLINE inline uint32_t sqr(uint32_t x) { return x * x; }
	LMAT_ENSURE_INLINE inline float    sqr(float    x) { return x * x; }
	LMAT_ENSURE_INLINE inline double   sqr(double   x) { return x * x; }

	LMAT_ENSURE_INLINE inline float  rcp(float  x) { return 1.0f / x; }
	LMAT_ENSURE_INLINE inline double rcp(double x) { return 1.0 / x; }

	LMAT_ENSURE_INLINE inline float  rsqrt(float  x) { return 1.0f / std::sqrt(x); }
	LMAT_ENSURE_INLINE inline double rsqrt(double x) { return 1.0 / std::sqrt(x); }

	LMAT_ENSURE_INLINE inline int32_t  (max)(int32_t  x, int32_t  y) { return x > y ? x : y; }
	LMAT_ENSURE_INLINE inline uint32_t (max)(uint32_t x, uint32_t y) { return x > y ? x : y; }
	LMAT_ENSURE_INLINE inline int32_t  (min)(int32_t  x, int32_t  y) { return x < y ? x : y; }
	LMAT_ENSURE_INLINE inline uint32_t (min)(uint32_t x, uint32_t y) { return x < y ? x : y; }

#ifdef LMAT_HAS_CXX11_MATH
	LMAT_ENSURE_INLINE inline float  (max)(float  x, float  y) { return std::fmax(x, y); }
	LMAT_ENSURE_INLINE inline double (max)(double x, double y) { return std::fmax(x, y); }
	LMAT_ENSURE_INLINE inline float  (min)(float  x, float  y) { return std::fmin(x, y); }
	LMAT_ENSURE_INLINE inline double (min)(double x, double y) { return std::fmin(x, y); }
#else
	LMAT_ENSURE_INLINE inline float  (max)(float  x, float  y) { return x > y ? x : y; }
	LMAT_ENSURE_INLINE inline double (max)(double x, double y) { return x > y ? x : y; }
	LMAT_ENSURE_INLINE inline float  (min)(float  x, float  y) { return x < y ? x : y; }
	LMAT_ENSURE_INLINE inline double (min)(double x, double y) { return x < y ? x : y; }
#endif


	// import of C++11 standard math functions

#ifdef LMAT_HAS_CXX11_MATH

	using std::cbrt;
	using std::hypot;

	using std::exp2;
	using std::log2;
	using std::expm1;
	using std::log1p;

	using std::asinh;
	using std::acosh;
	using std::atanh;

	using std::erf;
	using std::erfc;

	using std::trunc;
	using std::round;

	using std::isinf;
	using std::isnan;
	using std::isfinite;

#endif


} }

#endif 
