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


#define _LMAT_DEFINE_MAT_ARITH_FUN_1( FunName, FTag ) \
	template<typename T, class X> \
	LMAT_ENSURE_INLINE \
	inline map_expr<ftags::FTag, X> \
	FunName (const IEWiseMatrix<X, T>& x) \
	{ return make_map_expr(ftags::FTag(), x); }

#define _LMAT_DEFINE_MAT_INPLACE_ARITH( FunName, FTag ) \
	template<typename T, class X, class Y> \
	LMAT_ENSURE_INLINE \
	inline X& FunName (IEWiseMatrix<X, T>& x, const IEWiseMatrix<Y, T>& y) \
	{ x.derived() = make_map_expr(ftags::FTag(), x, y); return x.derived(); } \
	template<typename T, class X> \
	LMAT_ENSURE_INLINE \
	inline X& FunName (IEWiseMatrix<X, T>& x, const T& y) \
	{ x.derived() = make_map_expr_fix2(ftags::FTag(), x, y); return x.derived(); }

#define _LMAT_DEFINE_MAT_ARITH_FUN_2( FunName, FTag ) \
	template<typename T, class X, class Y> \
	LMAT_ENSURE_INLINE \
	inline map_expr<ftags::FTag, X, Y> \
	FunName (const IEWiseMatrix<X, T>& x, const IEWiseMatrix<Y, T>& y) \
	{ return make_map_expr(ftags::FTag(), x, y); } \
	template<typename T, class X> \
	LMAT_ENSURE_INLINE \
	inline map_expr<ftags::FTag, X, T> \
	FunName (const IEWiseMatrix<X, T>& x, const T& y) \
	{ return make_map_expr_fix2(ftags::FTag(), x, y); } \
	template<typename T, class Y> \
	LMAT_ENSURE_INLINE \
	inline map_expr<ftags::FTag, T, Y> \
	FunName (const T& x, const IEWiseMatrix<Y, T>& y) \
	{ return make_map_expr_fix1(ftags::FTag(), x, y); }


namespace lmat
{

	_LMAT_DEFINE_MAT_ARITH_FUN_2( operator +, add_ )
	_LMAT_DEFINE_MAT_ARITH_FUN_2( operator -, sub_ )
	_LMAT_DEFINE_MAT_ARITH_FUN_2( operator *, mul_ )
	_LMAT_DEFINE_MAT_ARITH_FUN_2( operator /, div_ )
	_LMAT_DEFINE_MAT_ARITH_FUN_1( operator -, neg_ )

	_LMAT_DEFINE_MAT_INPLACE_ARITH( operator +=, add_ )
	_LMAT_DEFINE_MAT_INPLACE_ARITH( operator -=, sub_ )
	_LMAT_DEFINE_MAT_INPLACE_ARITH( operator *=, mul_ )
	_LMAT_DEFINE_MAT_INPLACE_ARITH( operator /=, div_ )

	_LMAT_DEFINE_MAT_ARITH_FUN_1( abs, abs_ )
	_LMAT_DEFINE_MAT_ARITH_FUN_1( sqr, sqr_ )
	_LMAT_DEFINE_MAT_ARITH_FUN_1( cube, cube_ )

	_LMAT_DEFINE_MAT_ARITH_FUN_2( max, max_ )
	_LMAT_DEFINE_MAT_ARITH_FUN_2( min, min_ )

}

#endif /* MAT_ARITH_H_ */
