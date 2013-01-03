/**
 * @file mat_arith.h
 *
 * @brief Matrix arithmetics
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MAT_ARITH_H_
#define LIGHTMAT_MAT_ARITH_H_

#include <light_mat/mateval/map_expr.h>


#define LMAT_DEFINE_MAT_ARITH_FUN_1( FunName, FTag ) \
	template<typename T, class X> \
	LMAT_ENSURE_INLINE \
	inline map_expr<FTag, X> \
	FunName (const IEWiseMatrix<X, T>& x) \
	{ return make_map_expr(FTag(), x); }

#define LMAT_DEFINE_MAT_INPLACE_ARITH( FunName, FTag ) \
	template<typename T, class X, class Y> \
	LMAT_ENSURE_INLINE \
	inline X& FunName (IEWiseMatrix<X, T>& x, const IEWiseMatrix<Y, T>& y) \
	{ x.derived() = make_map_expr(FTag(), x, y); return x.derived(); } \
	template<typename T, class X> \
	LMAT_ENSURE_INLINE \
	inline X& FunName (IEWiseMatrix<X, T>& x, const T& y) \
	{ x.derived() = make_map_expr_fix2(FTag(), x, y); return x.derived(); }

#define LMAT_DEFINE_MAT_ARITH_FUN_2( FunName, FTag ) \
	template<typename T, class X, class Y> \
	LMAT_ENSURE_INLINE \
	inline map_expr<FTag, X, Y> \
	FunName (const IEWiseMatrix<X, T>& x, const IEWiseMatrix<Y, T>& y) \
	{ return make_map_expr(FTag(), x, y); } \
	template<typename T, class X> \
	LMAT_ENSURE_INLINE \
	inline map_expr<FTag, X, T> \
	FunName (const IEWiseMatrix<X, T>& x, const T& y) \
	{ return make_map_expr_fix2(FTag(), x, y); } \
	template<typename T, class Y> \
	LMAT_ENSURE_INLINE \
	inline map_expr<FTag, T, Y> \
	FunName (const T& x, const IEWiseMatrix<Y, T>& y) \
	{ return make_map_expr_fix1(FTag(), x, y); }


namespace lmat
{

	LMAT_DEFINE_MAT_ARITH_FUN_2( operator +, add_ )
	LMAT_DEFINE_MAT_ARITH_FUN_2( operator -, sub_ )
	LMAT_DEFINE_MAT_ARITH_FUN_2( operator *, mul_ )
	LMAT_DEFINE_MAT_ARITH_FUN_2( operator /, div_ )
	LMAT_DEFINE_MAT_ARITH_FUN_1( operator -, neg_ )

	LMAT_DEFINE_MAT_INPLACE_ARITH( operator +=, add_ )
	LMAT_DEFINE_MAT_INPLACE_ARITH( operator -=, sub_ )
	LMAT_DEFINE_MAT_INPLACE_ARITH( operator *=, mul_ )
	LMAT_DEFINE_MAT_INPLACE_ARITH( operator /=, div_ )

	LMAT_DEFINE_MAT_ARITH_FUN_1( abs, abs_ )
	LMAT_DEFINE_MAT_ARITH_FUN_1( sqr, sqr_ )
	LMAT_DEFINE_MAT_ARITH_FUN_1( cube, cube_ )

	LMAT_DEFINE_MAT_ARITH_FUN_2( max, max_ )
	LMAT_DEFINE_MAT_ARITH_FUN_2( min, min_ )

}

#endif /* MAT_ARITH_H_ */
