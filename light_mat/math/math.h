/**
 * @file math.h
 *
 * Definition/import of C99 math functions
 * 
 * @author Dahua Lin 
 */

#ifdef _MSC_VER
#pragma once
#endif

#include <light_mat/math/math_base.h>

#ifndef LIGHTMAT_MATH_H_
#define LIGHTMAT_MATH_H_

namespace lmat { namespace math {

	// power functions

	using std::pow;
	using std::fmod;

	using std::cbrt;
	using std::hypot;

	// exp & log

	using std::exp;
	using std::log;
	using std::log10;

	using std::exp2;
	using std::log2;
	using std::expm1;
	using std::log1p;

	// xlogy

	LMAT_ENSURE_INLINE inline float xlogy(float x, float y)
	{
		return x > 0 ? x * log(y) : 0.f;
	}

	LMAT_ENSURE_INLINE inline double xlogy(double x, double y)
	{
		return x > 0 ? x * log(y) : 0.0;
	}

	LMAT_ENSURE_INLINE inline float xlogx(float x)
	{
		return x > 0 ? x * log(x) : 0.f;
	}

	LMAT_ENSURE_INLINE inline double xlogx(double x)
	{
		return x > 0 ? x * log(x) : 0.0;
	}

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

	using std::asinh;
	using std::acosh;
	using std::atanh;

} }

#endif 
