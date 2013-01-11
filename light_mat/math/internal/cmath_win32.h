/**
 * @file cmath_win32.h
 *
 * Declare C99 functions that are lacking for WIN32
 * 
 * @author Dahua Lin 
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_CMATH_WIN32_H_
#define LIGHTMAT_CMATH_WIN32_H_

#include <light_mat/common/prim_types.h>

#if (LIGHTMAT_PLATFORM != LIGHTMAT_WIN32)
#error "Only include this file in WIN32 platform"
#endif

#define _LMAT_DECLARE_C99_MATH_1(Name) \
	float Name##f(float); \
	double Name(double);

#define _LMAT_DECLARE_C99_MATH_2(Name) \
	float Name##f(float, float); \
	double Name(double, double);

#define _LMAT_IMPORT_C99_MATH_1(Name) \
	LMAT_ENSURE_INLINE \
	inline float Name(float x) { return ::Name##f(x); } \
	LMAT_ENSURE_INLINE \
	inline double Name(double x) { return ::Name(x); }

#define _LMAT_IMPORT_C99_MATH_2(Name) \
	LMAT_ENSURE_INLINE \
	inline float Name(float x, float y) { return ::Name##f(x, y); } \
	LMAT_ENSURE_INLINE \
	inline double Name(double x, double y) { return ::Name(x, y); }

extern "C"
{
	_LMAT_DECLARE_C99_MATH_1( cbrt )
	_LMAT_DECLARE_C99_MATH_2( hypot )
}

// Import into std namespace

namespace std
{
	_LMAT_IMPORT_C99_MATH_1( cbrt )
	_LMAT_IMPORT_C99_MATH_2( hypot )

}

#endif 
