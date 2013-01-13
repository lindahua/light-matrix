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

#include <light_mat/matexpr/matfun_base.h>
#include <light_mat/math/basic_functors.h>

namespace lmat
{

	// arithmetics

	_LMAT_DEFINE_GMATOP2( add, + )
	_LMAT_DEFINE_GMATOP2( sub, - )
	_LMAT_DEFINE_GMATOP2( mul, * )
	_LMAT_DEFINE_GMATOP2( div, / )

	_LMAT_DEFINE_GMATOP( neg, -, 1 )

	_LMAT_DEFINE_RMATFUN( fma, 3 )

	// min & max

	_LMAT_DEFINE_GMATFUN( max, 2 )
	_LMAT_DEFINE_GMATFUN( min, 2 )
	_LMAT_DEFINE_GMATFUN( clamp, 3 )

	// simple power functions

	_LMAT_DEFINE_GMATFUN( abs, 1 )
	_LMAT_DEFINE_GMATFUN( sqr, 1 )
	_LMAT_DEFINE_RMATFUN( cube, 1 )
	_LMAT_DEFINE_RMATFUN( rcp, 1 )
	_LMAT_DEFINE_RMATFUN( sqrt, 1 )
	_LMAT_DEFINE_RMATFUN( rsqrt, 1 )

	// rounding

	_LMAT_DEFINE_RMATFUN( floor, 1 )
	_LMAT_DEFINE_RMATFUN( ceil, 1 )
	_LMAT_DEFINE_RMATFUN( round, 1 )
	_LMAT_DEFINE_RMATFUN( trunc, 1 )

	// conditional operator

	template<typename T, class C, class X, class Y>
	LMAT_ENSURE_INLINE
	inline map_expr<ftags::cond_, C, X, Y>
	cond(const IEWiseMatrix<C, bool>& c, const IEWiseMatrix<X, T>& x, const IEWiseMatrix<Y, T>& y)
	{
		return make_map_expr(ftags::cond_(), c, x, y);
	}

	template<typename T, class C, class X, class Y>
	LMAT_ENSURE_INLINE
	inline map_expr<ftags::cond_, C, X, Y>
	cond(const IEWiseMatrix<C, mask_t<T> >& c, const IEWiseMatrix<X, T>& x, const IEWiseMatrix<Y, T>& y)
	{
		return make_map_expr(ftags::cond_(), c, x, y);
	}

	template<typename T, class C, class X>
	LMAT_ENSURE_INLINE
	inline map_expr<ftags::cond_, C, X, T>
	cond(const IEWiseMatrix<C, bool>& c, const IEWiseMatrix<X, T>& x, const T& y)
	{
		return make_map_expr_fix3(ftags::cond_(), c, x, y);
	}

	template<typename T, class C, class X>
	LMAT_ENSURE_INLINE
	inline map_expr<ftags::cond_, C, X, T>
	cond(const IEWiseMatrix<C, mask_t<T> >& c, const IEWiseMatrix<X, T>& x, const T& y)
	{
		return make_map_expr_fix3(ftags::cond_(), c, x, y);
	}

	template<typename T, class C, class Y>
	LMAT_ENSURE_INLINE
	inline map_expr<ftags::cond_, C, T, Y>
	cond(const IEWiseMatrix<C, bool>& c, const T& x, const IEWiseMatrix<Y, T>& y)
	{
		return make_map_expr_fix2(ftags::cond_(), c, x, y);
	}

	template<typename T, class C, class Y>
	LMAT_ENSURE_INLINE
	inline map_expr<ftags::cond_, C, T, Y>
	cond(const IEWiseMatrix<C, mask_t<T> >& c, const T& x, const IEWiseMatrix<Y, T>& y)
	{
		return make_map_expr_fix2(ftags::cond_(), c, x, y);
	}

}

#endif /* MAT_ARITH_H_ */
