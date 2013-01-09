/**
 * @file mat_emath.h
 *
 * @brief Elemental math functions on matrices
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MAT_EMATH_H_
#define LIGHTMAT_MAT_EMATH_H_

#include <light_mat/mateval/map_expr.h>

#define _LMAT_DEFINE_MAT_MATH_FUN_S1( FunName, T, FTag ) \
	template<class X> \
	LMAT_ENSURE_INLINE \
	inline map_expr<ftags::FTag, X> \
	FunName (const IEWiseMatrix<X, T>& x) \
	{ return make_map_expr(ftags::FTag(), x); }

#define _LMAT_DEFINE_MAT_MATH_FUN_S2( FunName, T, FTag ) \
	template<class X, class Y> \
	LMAT_ENSURE_INLINE \
	inline map_expr<ftags::FTag, X, Y> \
	FunName (const IEWiseMatrix<X, T>& x, const IEWiseMatrix<Y, T>& y) \
	{ return make_map_expr(ftags::FTag(), x, y); } \
	template<class X> \
	LMAT_ENSURE_INLINE \
	inline map_expr<ftags::FTag, X, T> \
	FunName (const IEWiseMatrix<X, T>& x, const T& y) \
	{ return make_map_expr_fix2(ftags::FTag(), x, y); } \
	template<class Y> \
	LMAT_ENSURE_INLINE \
	inline map_expr<ftags::FTag, T, Y> \
	FunName (const T& x, const IEWiseMatrix<Y, T>& y) \
	{ return make_map_expr_fix1(ftags::FTag(), x, y); }

#define _LMAT_DEFINE_MAT_MATH_FUN_S3( FunName, T, FTag ) \
	template<class X, class Y, class Z> \
	LMAT_ENSURE_INLINE \
	inline map_expr<ftags::FTag, X, Y, Z> \
	FunName (const IEWiseMatrix<X, T>& x, const IEWiseMatrix<Y, T>& y, const IEWiseMatrix<Z, T>& z) \
	{ return make_map_expr(ftags::FTag(), x, y, z); } \
	template<class X, class Y> \
	LMAT_ENSURE_INLINE \
	inline map_expr<ftags::FTag, X, Y, T> \
	FunName (const IEWiseMatrix<X, T>& x, const IEWiseMatrix<Y, T>& y, const T& z) \
	{ return make_map_expr_fix3(ftags::FTag(), x, y, z); } \
	template<class X, class Z> \
	LMAT_ENSURE_INLINE \
	inline map_expr<ftags::FTag, X, T, Z> \
	FunName (const IEWiseMatrix<X, T>& x, const T& y, const IEWiseMatrix<Z, T>& z) \
	{ return make_map_expr_fix2(ftags::FTag(), x, y, z); } \
	template<class Y, class Z> \
	LMAT_ENSURE_INLINE \
	inline map_expr<ftags::FTag, T, Y, Z> \
	FunName (const T& x, const IEWiseMatrix<Y, T>& y, const IEWiseMatrix<Z, T>& z) \
	{ return make_map_expr_fix1(ftags::FTag(), x, y, z); } \
	template<class X> \
	LMAT_ENSURE_INLINE \
	inline map_expr<ftags::FTag, X, T, T> \
	FunName (const IEWiseMatrix<X, T>& x, const T& y, const T& z) \
	{ return make_map_expr_fix23(ftags::FTag(), x, y, z); } \
	template<class Y> \
	LMAT_ENSURE_INLINE \
	inline map_expr<ftags::FTag, T, Y, T> \
	FunName (const T& x, const IEWiseMatrix<Y, T>& y, const T& z) \
	{ return make_map_expr_fix13(ftags::FTag(), x, y, z); } \
	template<class Z> \
	LMAT_ENSURE_INLINE \
	inline map_expr<ftags::FTag, T, T, Z> \
	FunName (const T& x, const T& y, const IEWiseMatrix<Z, T>& z) \
	{ return make_map_expr_fix12(ftags::FTag(), x, y, z); }



#define _LMAT_DEFINE_MAT_MATH_FUN_R1( FunName ) \
		_LMAT_DEFINE_MAT_MATH_FUN_S1( FunName, float, FunName##_ ) \
		_LMAT_DEFINE_MAT_MATH_FUN_S1( FunName, double, FunName##_ )

#define _LMAT_DEFINE_MAT_MATH_FUN_R2( FunName ) \
		_LMAT_DEFINE_MAT_MATH_FUN_S2( FunName, float, FunName##_ ) \
		_LMAT_DEFINE_MAT_MATH_FUN_S2( FunName, double, FunName##_ )

#define _LMAT_DEFINE_MAT_MATH_FUN_R3( FunName ) \
		_LMAT_DEFINE_MAT_MATH_FUN_S3( FunName, float, FunName##_ ) \
		_LMAT_DEFINE_MAT_MATH_FUN_S3( FunName, double, FunName##_ )

namespace lmat
{
	_LMAT_DEFINE_MAT_MATH_FUN_R3( fma )
	_LMAT_DEFINE_MAT_MATH_FUN_R3( clamp )

	_LMAT_DEFINE_MAT_MATH_FUN_R1( rcp )
	_LMAT_DEFINE_MAT_MATH_FUN_R1( sqrt )
	_LMAT_DEFINE_MAT_MATH_FUN_R1( rsqrt )
	_LMAT_DEFINE_MAT_MATH_FUN_R2( pow )

	_LMAT_DEFINE_MAT_MATH_FUN_R1( floor )
	_LMAT_DEFINE_MAT_MATH_FUN_R1( ceil )

	_LMAT_DEFINE_MAT_MATH_FUN_R1( exp )
	_LMAT_DEFINE_MAT_MATH_FUN_R1( log )
	_LMAT_DEFINE_MAT_MATH_FUN_R1( log10 )
	_LMAT_DEFINE_MAT_MATH_FUN_R2( xlogy )
	_LMAT_DEFINE_MAT_MATH_FUN_R1( xlogx )

	_LMAT_DEFINE_MAT_MATH_FUN_R1( sin )
	_LMAT_DEFINE_MAT_MATH_FUN_R1( cos )
	_LMAT_DEFINE_MAT_MATH_FUN_R1( tan )

	_LMAT_DEFINE_MAT_MATH_FUN_R1( asin )
	_LMAT_DEFINE_MAT_MATH_FUN_R1( acos )
	_LMAT_DEFINE_MAT_MATH_FUN_R1( atan )
	_LMAT_DEFINE_MAT_MATH_FUN_R2( atan2 )

	_LMAT_DEFINE_MAT_MATH_FUN_R1( sinh )
	_LMAT_DEFINE_MAT_MATH_FUN_R1( cosh )
	_LMAT_DEFINE_MAT_MATH_FUN_R1( tanh )

#ifdef LMAT_HAS_CXX11_MATH

	_LMAT_DEFINE_MAT_MATH_FUN_R1( cbrt )
	_LMAT_DEFINE_MAT_MATH_FUN_R2( hypot )

	_LMAT_DEFINE_MAT_MATH_FUN_R1( round )
	_LMAT_DEFINE_MAT_MATH_FUN_R1( trunc )

	_LMAT_DEFINE_MAT_MATH_FUN_R1( exp2 )
	_LMAT_DEFINE_MAT_MATH_FUN_R1( log2 )
	_LMAT_DEFINE_MAT_MATH_FUN_R1( expm1 )
	_LMAT_DEFINE_MAT_MATH_FUN_R1( log1p )

	_LMAT_DEFINE_MAT_MATH_FUN_R1( asinh )
	_LMAT_DEFINE_MAT_MATH_FUN_R1( acosh )
	_LMAT_DEFINE_MAT_MATH_FUN_R1( atanh )

	_LMAT_DEFINE_MAT_MATH_FUN_R1( erf )
	_LMAT_DEFINE_MAT_MATH_FUN_R1( erfc )
	_LMAT_DEFINE_MAT_MATH_FUN_R1( lgamma )
	_LMAT_DEFINE_MAT_MATH_FUN_R1( tgamma )

#endif
}


#endif /* MAT_EMATH_H_ */
