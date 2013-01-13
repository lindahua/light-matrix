/**
 * @file mat_pred.h
 *
 * @brief Predicate functions on matrices
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MAT_PRED_H_
#define LIGHTMAT_MAT_PRED_H_

#include <light_mat/matexpr/mat_cast.h>
#include <light_mat/math/basic_functors.h>

#define _LMAT_DEFINE_MAT_LOGICAL_FUN_1( FunName, FTag ) \
	template<typename T, class X> \
	LMAT_ENSURE_INLINE \
	inline map_expr<ftags::FTag, X> \
	FunName (const IEWiseMatrix<X, mask_t<T> >& x) \
	{ return make_map_expr(ftags::FTag(), x); } \
	template<class X> \
	LMAT_ENSURE_INLINE \
	inline map_expr<ftags::FTag, X> \
	FunName (const IEWiseMatrix<X, bool >& x) \
	{ return make_map_expr(ftags::FTag(), x); }

#define _LMAT_DEFINE_MAT_LOGICAL_FUN_2( FunName, FTag ) \
	template<typename T, class X, class Y> \
	LMAT_ENSURE_INLINE \
	inline map_expr<ftags::FTag, X, Y> \
	FunName (const IEWiseMatrix<X, mask_t<T> >& x, const IEWiseMatrix<Y, mask_t<T> >& y) \
	{ return make_map_expr(ftags::FTag(), x, y); } \
	template<typename T, class X, class Y> \
	LMAT_ENSURE_INLINE \
	inline map_expr<ftags::FTag, X, Y> \
	FunName (const IEWiseMatrix<X, mask_t<T> >& x, const IEWiseMatrix<Y, bool>& y) \
	{ return make_map_expr(ftags::FTag(), x, y); } \
	template<typename T, class X, class Y> \
	LMAT_ENSURE_INLINE \
	inline map_expr<ftags::FTag, X, Y> \
	FunName (const IEWiseMatrix<X, bool>& x, const IEWiseMatrix<Y, mask_t<T> >& y) \
	{ return make_map_expr(ftags::FTag(), x, y); } \
	template<class X, class Y> \
	LMAT_ENSURE_INLINE \
	inline map_expr<ftags::FTag, X, Y> \
	FunName (const IEWiseMatrix<X, bool>& x, const IEWiseMatrix<Y, bool>& y) \
	{ return make_map_expr(ftags::FTag(), x, y); } \


namespace lmat
{
	// comparison

	_LMAT_DEFINE_GMATOP( eq, ==, 2 )
	_LMAT_DEFINE_GMATOP( ne, !=, 2 )
	_LMAT_DEFINE_GMATOP( lt, <,  2 )
	_LMAT_DEFINE_GMATOP( le, <=, 2 )
	_LMAT_DEFINE_GMATOP( gt, >,  2 )
	_LMAT_DEFINE_GMATOP( ge, >=, 2 )

	_LMAT_DEFINE_RMATFUN( signbit, 1 )
	_LMAT_DEFINE_RMATFUN( isfinite, 1 )
	_LMAT_DEFINE_RMATFUN( isinf, 1 )
	_LMAT_DEFINE_RMATFUN( isnan, 1 )


	// logical

	_LMAT_DEFINE_MAT_LOGICAL_FUN_1( operator ~,  logical_not_ )
	_LMAT_DEFINE_MAT_LOGICAL_FUN_2( operator &,  logical_and_ )
	_LMAT_DEFINE_MAT_LOGICAL_FUN_2( operator |,  logical_or_ )
	_LMAT_DEFINE_MAT_LOGICAL_FUN_2( operator ==, logical_eq_ )
	_LMAT_DEFINE_MAT_LOGICAL_FUN_2( operator !=, logical_ne_ )

}

#endif /* MAT_PRED_H_ */
