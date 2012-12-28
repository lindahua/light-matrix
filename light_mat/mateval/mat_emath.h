/**
 * @file mat_emath.h
 *
 * @brief Elemental math functions on matrices
 *
 * @author Dahua Lin
 */

#ifndef LIGHTMAT_MAT_EMATH_H_
#define LIGHTMAT_MAT_EMATH_H_

#include <light_mat/mateval/map_expr.h>

#define LMAT_DEFINE_MAT_MATH_FUN_S1( FunName, T, FTag ) \
	template<class X> \
	LMAT_ENSURE_INLINE \
	inline map_expr<FTag, X> \
	FunName (const IEWiseMatrix<X, T>& x) \
	{ return make_map_expr(FTag(), x); }

#define LMAT_DEFINE_MAT_MATH_FUN_S2( FunName, T, FTag ) \
	template<class X, class Y> \
	LMAT_ENSURE_INLINE \
	inline map_expr<FTag, X, Y> \
	FunName (const IEWiseMatrix<X, T>& x, const IEWiseMatrix<Y, T>& y) \
	{ return make_map_expr(FTag(), x, y); } \
	template<class X> \
	LMAT_ENSURE_INLINE \
	inline map_expr<FTag, X, T> \
	FunName (const IEWiseMatrix<X, T>& x, const T& y) \
	{ return make_map_expr_fix2(FTag(), x, y); } \
	template<class Y> \
	LMAT_ENSURE_INLINE \
	inline map_expr<FTag, T, Y> \
	FunName (const T& x, const IEWiseMatrix<Y, T>& y) \
	{ return make_map_expr_fix1(FTag(), x, y); }

#define LMAT_DEFINE_MAT_MATH_FUN_S3( FunName, T, FTag ) \
	template<class X, class Y, class Z> \
	LMAT_ENSURE_INLINE \
	inline map_expr<FTag, X, Y, Z> \
	FunName (const IEWiseMatrix<X, T>& x, const IEWiseMatrix<Y, T>& y, const IEWiseMatrix<Z, T>& z) \
	{ return make_map_expr(FTag(), x, y, z); } \
	template<class X, class Y> \
	LMAT_ENSURE_INLINE \
	inline map_expr<FTag, X, Y, T> \
	FunName (const IEWiseMatrix<X, T>& x, const IEWiseMatrix<Y, T>& y, const T& z) \
	{ return make_map_expr_fix3(FTag(), x, y, z); } \
	template<class X, class Z> \
	LMAT_ENSURE_INLINE \
	inline map_expr<FTag, X, T, Z> \
	FunName (const IEWiseMatrix<X, T>& x, const T& y, const IEWiseMatrix<Z, T>& z) \
	{ return make_map_expr_fix2(FTag(), x, y, z); } \
	template<class Y, class Z> \
	LMAT_ENSURE_INLINE \
	inline map_expr<FTag, T, Y, Z> \
	FunName (const T& x, const IEWiseMatrix<Y, T>& y, const IEWiseMatrix<Z, T>& z) \
	{ return make_map_expr_fix1(FTag(), x, y, z); } \
	template<class X> \
	LMAT_ENSURE_INLINE \
	inline map_expr<FTag, X, T, T> \
	FunName (const IEWiseMatrix<X, T>& x, const T& y, const T& z) \
	{ return make_map_expr_fix23(FTag(), x, y, z); } \
	template<class Y> \
	LMAT_ENSURE_INLINE \
	inline map_expr<FTag, T, Y, T> \
	FunName (const T& x, const IEWiseMatrix<Y, T>& y, const T& z) \
	{ return make_map_expr_fix13(FTag(), x, y, z); } \
	template<class Z> \
	LMAT_ENSURE_INLINE \
	inline map_expr<FTag, T, T, Z> \
	FunName (const T& x, const T& y, const IEWiseMatrix<Z, T>& z) \
	{ return make_map_expr_fix12(FTag(), x, y, z); }



#define LMAT_DEFINE_MAT_MATH_FUN_R1( FunName ) \
		LMAT_DEFINE_MAT_MATH_FUN_S1( FunName, float, FunName##_ ) \
		LMAT_DEFINE_MAT_MATH_FUN_S1( FunName, double, FunName##_ )

#define LMAT_DEFINE_MAT_MATH_FUN_R2( FunName ) \
		LMAT_DEFINE_MAT_MATH_FUN_S2( FunName, float, FunName##_ ) \
		LMAT_DEFINE_MAT_MATH_FUN_S2( FunName, double, FunName##_ )

#define LMAT_DEFINE_MAT_MATH_FUN_R3( FunName ) \
		LMAT_DEFINE_MAT_MATH_FUN_S3( FunName, float, FunName##_ ) \
		LMAT_DEFINE_MAT_MATH_FUN_S3( FunName, double, FunName##_ )

namespace lmat
{
	LMAT_DEFINE_MAT_MATH_FUN_R3( fma )
	LMAT_DEFINE_MAT_MATH_FUN_R3( clamp )

	LMAT_DEFINE_MAT_MATH_FUN_R1( rcp )
	LMAT_DEFINE_MAT_MATH_FUN_R1( sqrt )
	LMAT_DEFINE_MAT_MATH_FUN_R1( rsqrt )
	LMAT_DEFINE_MAT_MATH_FUN_R2( pow )

	LMAT_DEFINE_MAT_MATH_FUN_R1( floor )
	LMAT_DEFINE_MAT_MATH_FUN_R1( ceil )

	LMAT_DEFINE_MAT_MATH_FUN_R1( exp )
	LMAT_DEFINE_MAT_MATH_FUN_R1( log )
	LMAT_DEFINE_MAT_MATH_FUN_R1( log10 )
	LMAT_DEFINE_MAT_MATH_FUN_R2( xlogy )
	LMAT_DEFINE_MAT_MATH_FUN_R1( xlogx )

	LMAT_DEFINE_MAT_MATH_FUN_R1( sin )
	LMAT_DEFINE_MAT_MATH_FUN_R1( cos )
	LMAT_DEFINE_MAT_MATH_FUN_R1( tan )

	LMAT_DEFINE_MAT_MATH_FUN_R1( asin )
	LMAT_DEFINE_MAT_MATH_FUN_R1( acos )
	LMAT_DEFINE_MAT_MATH_FUN_R1( atan )
	LMAT_DEFINE_MAT_MATH_FUN_R2( atan2 )

	LMAT_DEFINE_MAT_MATH_FUN_R1( sinh )
	LMAT_DEFINE_MAT_MATH_FUN_R1( cosh )
	LMAT_DEFINE_MAT_MATH_FUN_R1( tanh )

#ifdef LMAT_HAS_CXX11_MATH

	LMAT_DEFINE_MAT_MATH_FUN_R1( cbrt )
	LMAT_DEFINE_MAT_MATH_FUN_R2( hypot )

	LMAT_DEFINE_MAT_MATH_FUN_R1( round )
	LMAT_DEFINE_MAT_MATH_FUN_R1( trunc )

	LMAT_DEFINE_MAT_MATH_FUN_R1( exp2 )
	LMAT_DEFINE_MAT_MATH_FUN_R1( log2 )
	LMAT_DEFINE_MAT_MATH_FUN_R1( expm1 )
	LMAT_DEFINE_MAT_MATH_FUN_R1( log1p )

	LMAT_DEFINE_MAT_MATH_FUN_R1( asinh )
	LMAT_DEFINE_MAT_MATH_FUN_R1( acosh )
	LMAT_DEFINE_MAT_MATH_FUN_R1( atanh )

	LMAT_DEFINE_MAT_MATH_FUN_R1( erf )
	LMAT_DEFINE_MAT_MATH_FUN_R1( erfc )
	LMAT_DEFINE_MAT_MATH_FUN_R1( lgamma )
	LMAT_DEFINE_MAT_MATH_FUN_R1( tgamma )

#endif
}


#endif /* MAT_EMATH_H_ */
