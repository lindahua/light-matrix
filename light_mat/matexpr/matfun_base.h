/**
 * @file matfun_base.h
 *
 * @brief The basic facilities to help defining functions on matrices
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATFUN_BASE_H_
#define LIGHTMAT_MATFUN_BASE_H_

#include <light_mat/matexpr/map_expr.h>

/************************************************
 *
 *  Internal macros
 *
 ************************************************/

// generic functions

#define _LMAT_DEFINE_GMATFUN_1( FunName, FTag ) \
	template<typename T, class X> \
	LMAT_ENSURE_INLINE \
	inline map_expr<FTag, X> \
	FunName (const IEWiseMatrix<X, T>& x) \
	{ return make_map_expr(FTag(), x); }

#define _LMAT_DEFINE_GMATFUN_2( FunName, FTag ) \
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

#define _LMAT_DEFINE_GMATFUN_3( FunName, FTag ) \
	template<typename T, class X, class Y, class Z> \
	LMAT_ENSURE_INLINE \
	inline map_expr<FTag, X, Y, Z> \
	FunName (const IEWiseMatrix<X, T>& x, const IEWiseMatrix<Y, T>& y, const IEWiseMatrix<Z, T>& z) \
	{ return make_map_expr(FTag(), x, y, z); } \
	template<typename T, class X, class Y> \
	LMAT_ENSURE_INLINE \
	inline map_expr<FTag, X, Y, T> \
	FunName (const IEWiseMatrix<X, T>& x, const IEWiseMatrix<Y, T>& y, const T& z) \
	{ return make_map_expr_fix3(FTag(), x, y, z); } \
	template<typename T, class X, class Z> \
	LMAT_ENSURE_INLINE \
	inline map_expr<FTag, X, T, Z> \
	FunName (const IEWiseMatrix<X, T>& x, const T& y, const IEWiseMatrix<Z, T>& z) \
	{ return make_map_expr_fix2(FTag(), x, y, z); } \
	template<typename T, class Y, class Z> \
	LMAT_ENSURE_INLINE \
	inline map_expr<FTag, T, Y, Z> \
	FunName (const T& x, const IEWiseMatrix<Y, T>& y, const IEWiseMatrix<Z, T>& z) \
	{ return make_map_expr_fix1(FTag(), x, y, z); } \
	template<typename T, class X> \
	LMAT_ENSURE_INLINE \
	inline map_expr<FTag, X, T, T> \
	FunName (const IEWiseMatrix<X, T>& x, const T& y, const T& z) \
	{ return make_map_expr_fix23(FTag(), x, y, z); } \
	template<typename T, class Y> \
	LMAT_ENSURE_INLINE \
	inline map_expr<FTag, T, Y, T> \
	FunName (const T& x, const IEWiseMatrix<Y, T>& y, const T& z) \
	{ return make_map_expr_fix13(FTag(), x, y, z); } \
	template<typename T, class Z> \
	LMAT_ENSURE_INLINE \
	inline map_expr<FTag, T, T, Z> \
	FunName (const T& x, const T& y, const IEWiseMatrix<Z, T>& z) \
	{ return make_map_expr_fix12(FTag(), x, y, z); }

// specific functions

#define _LMAT_DEFINE_SMATFUN_1( FunName, FTag, T ) \
	template<class X> \
	LMAT_ENSURE_INLINE \
	inline map_expr<FTag, X> \
	FunName (const IEWiseMatrix<X, T>& x) \
	{ return make_map_expr(FTag(), x); }

#define _LMAT_DEFINE_SMATFUN_2( FunName, FTag, T ) \
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

#define _LMAT_DEFINE_SMATFUN_3( FunName, FTag, T ) \
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


/************************************************
 *
 *  Macros to be used by developers
 *
 ************************************************/

#define LMAT_DEF_GMATFUN( FunName, FTag, NA ) _LMAT_DEFINE_GMATFUN_##NA( FunName, FTag )

#define LMAT_DEF_SMATFUN( FunName, FTag, T, NA ) _LMAT_DEFINE_SMATFUN_##NA( FunName, FTag, T )

#define LMAT_DEF_RMATFUN( FunName, FTag, NA ) \
	_LMAT_DEFINE_SMATFUN_##NA( FunName, FTag, float ) \
	_LMAT_DEFINE_SMATFUN_##NA( FunName, FTag, double )


/************************************************
 *
 *  Macros specialized for LightMatrix
 *
 ************************************************/

#define _LMAT_DEFINE_GMATOP( Name, Op, NA ) LMAT_DEF_GMATFUN( operator Op, ftags::Name##_, NA )

#define _LMAT_DEFINE_GMATOP2( Name, Op ) \
	LMAT_DEF_GMATFUN( operator Op, ftags::Name##_, 2 ) \
	template<typename T, class X, class Y> \
	LMAT_ENSURE_INLINE \
	inline X& operator Op##= (IRegularMatrix<X, T>& x, const IEWiseMatrix<Y, T>& y) { \
		x.derived() = x Op y; \
		return x.derived(); } \
	template<typename T, class X> \
	LMAT_ENSURE_INLINE \
	inline X& operator Op##= (IRegularMatrix<X, T>& x, const T& y) { \
		x.derived() = x Op y; \
		return x.derived(); } \

#define _LMAT_DEFINE_GMATFUN( Name, NA ) LMAT_DEF_GMATFUN( Name, ftags::Name##_, NA )
#define _LMAT_DEFINE_RMATFUN( Name, NA ) LMAT_DEF_RMATFUN( Name, ftags::Name##_, NA )

#endif
